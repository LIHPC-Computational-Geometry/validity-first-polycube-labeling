#pragma once

#include <geogram/mesh/mesh.h>          // for Mesh
#include <geogram/mesh/mesh_geometry.h> // for mesh_facet_normal()
#include <geogram/basic/numeric.h>      // for index_t
#include <geogram/basic/vecg.h>         // for normalize()
#include <geogram/basic/assert.h>       // for geo_assert()

#include <vector>
#include <set>
#include <utility> // for std::pair
#include <iterator>
#include <algorithm> // for std::swap()

#include "geometry_halfedges.h" // for MeshHalfedgesExt

using namespace GEO;

#define FEATURE_EDGES_MIN_ANGLE	0.5

class MeshExtFacetNormals {

public:
    MeshExtFacetNormals(const Mesh& mesh) : mesh_(mesh) {
        recompute();
    }

    void clear() {
        normals_.clear();
    }

    void recompute() {
        normals_.resize(mesh_.facets.nb());
        FOR(f,mesh_.facets.nb()) {
            normals_[f] = normalize(Geom::mesh_facet_normal(mesh_,f));
        }
    }

    std::vector<vec3>& as_vector() {
        return normals_;
    }

    const std::vector<vec3>& as_vector() const {
        return normals_;
    }

    vec3& operator[](index_t facet_index) {
        return normals_[facet_index];
    }

    const vec3& operator[](index_t facet_index) const {
        return normals_[facet_index];
    }

protected:

    std::vector<vec3> normals_; // facet normals
    const Mesh& mesh_;
};

struct FacetCorner {
    index_t facet_index;
    index_t local_vertex;
};

class MeshExtAdjFacetCorners {
public:
    MeshExtAdjFacetCorners(const Mesh& mesh) : mesh_(mesh) {
        recompute();
    };

    void clear() {
        adj_facets_corners_.clear();
    }

    void recompute() {
        adj_facets_corners_.resize(mesh_.vertices.nb());
        // link each facet to adjacent vertices
        FOR(f,mesh_.facets.nb()) {
            FOR(lv,mesh_.facets.nb_vertices(f)) { // for each local vertices of the facet f
                adj_facets_corners_[mesh_.facets.vertex(f,lv)].push_back(FacetCorner{f,lv});
            }
        }
    }

    std::vector<std::vector<FacetCorner>>& as_vector() {
        return adj_facets_corners_;
    }

    const std::vector<std::vector<FacetCorner>>& as_vector() const {
        return adj_facets_corners_;
    }

    std::vector<FacetCorner>& of_vertex(index_t vertex_index) {
        return adj_facets_corners_[vertex_index];
    }

    const std::vector<FacetCorner>& of_vertex(index_t vertex_index) const {
        return adj_facets_corners_[vertex_index];
    }

    bool size_matches_nb_vertices() const {
        return adj_facets_corners_.size() == (size_t) mesh_.vertices.nb();
    }

protected:

    std::vector<std::vector<FacetCorner>> adj_facets_corners_; // for each vertex, store adjacent facets. no ordering. has mesh_.vertices.nb() elements, each of them is a list of adjacent facet, with 2 components : the facet index, and the local vertex index for this facet
    const Mesh& mesh_;
};

class MeshExtFeatureEdges {

public:

    MeshExtFeatureEdges(const Mesh& mesh, const MeshExtAdjFacetCorners& adj_facet_corners) : mesh_(mesh), adj_facet_corners_(adj_facet_corners) {
        recompute();
    };

    void clear() {
        feature_edges_.clear();
    }

    // parse mesh_.edges
    // filter edges where angle >= FEATURE_EDGES_MIN_ANGLE
    // store them in a set of pair of vertices
    // do not edit mesh_.edges
    void recompute() {
        geo_assert(adj_facet_corners_.size_matches_nb_vertices()); // assert the adjacency has been computed

        index_t nb_edges = mesh_.edges.nb(),
                v0_index = index_t(-1),
                v1_index = index_t(-1),
                facet_corner_on_v0 = index_t(-1),
                facet_at_left = index_t(-1),
                facet_at_right = index_t(-1);
        unsigned int nb_edges_removed = 0;
        MeshHalfedgesExt::Halfedge outgoing_halfedge;

        FOR(e,nb_edges) {
            // get indices of the 2 vertices of the current edge
            v0_index = mesh_.edges.vertex(e,0);
            v1_index = mesh_.edges.vertex(e,1);
            // parse adjacent facets of v0 until we found an halfedge
            // whose origin is on v0 and tip on v1
            for(const FacetCorner& fc : adj_facet_corners_.of_vertex(v0_index)) {
                FOR(lv,mesh_.facets.nb_vertices(fc.facet_index)) {
                    facet_corner_on_v0 = mesh_.facets.corner(fc.facet_index,lv);
                    if(mesh_.facet_corners.vertex(facet_corner_on_v0) == v0_index) {
                        break; // we found which facet corner is on v0
                    }
                }
                // now that we have a facet and a facet cornet, we can have an halfedge going out of v0
                outgoing_halfedge.facet = fc.facet_index;
                outgoing_halfedge.corner = facet_corner_on_v0;
                geo_assert(halfedge_vertex_index_from(mesh_,outgoing_halfedge) == v0_index)
                // check of the vertex at the tip is v1
                if(halfedge_vertex_index_to(mesh_,outgoing_halfedge) == v1_index) {
                    break;
                }
                // else : continue and check another adjacent facet
            }
            facet_at_left = halfedge_facet_left(mesh_,outgoing_halfedge);
            facet_at_right = halfedge_facet_right(mesh_,outgoing_halfedge);
            if(angle(mesh_facet_normal(mesh_,facet_at_left),mesh_facet_normal(mesh_,facet_at_right)) < FEATURE_EDGES_MIN_ANGLE) {
                nb_edges_removed++;
            }
            else {
                // Keep this feature edge because significant angle

                // Use a pair of vertices and not GEO::Mesh.edges, because when working on the labeling graph
                // (charts, boundaries & corners) we use oriented edges,
                // and we need to know if a given oriented edge is on a feature edge
                // When feature edges are stored in GEO::Mesh.edges, see https://github.com/BrunoLevy/geogram/wiki/Mesh#mesh-edges
                // finding if they contains (v0,v1) is expensive, we have to check all edges
                // Instead we can use a set of pair of vertices, the pair being sorted by ascending index,
                // to quickly check if (v0,v1) is a feature edge with `feature_edges_.contains(min(v0,v1),max(v0,v1))`
                feature_edges_.insert(std::make_pair(
                    std::min(v0_index,v1_index),
                    std::max(v0_index,v1_index)
                ));
            }
        }

        fmt::println(
            Logger::out("feature edges"),
            "{} feature edge(s) removed {} because dihedral angle < {}",
            nb_edges_removed,
            (nb_edges==0 ? "" : fmt::format("({:.1f}%)",(nb_edges_removed * 100.0) / (double) nb_edges)),
            FEATURE_EDGES_MIN_ANGLE);
        Logger::out("feature edges").flush();
    }

    index_t nb() const {
        return (index_t) feature_edges_.size();
    }

    bool contain(std::pair<index_t,index_t> pair_of_vertices) const {
        // ensure the first vertex has smaller index that the second
        if(pair_of_vertices.first > pair_of_vertices.second) {
            std::swap(pair_of_vertices.first,pair_of_vertices.second);
        }
        return feature_edges_.contains(pair_of_vertices);
    }

    bool contain_halfedge(MeshHalfedges::Halfedge halfedge) const {
        index_t v0 = halfedge_vertex_index_from(mesh_,halfedge),
                v1 = halfedge_vertex_index_to(mesh_,halfedge);
        return feature_edges_.contains(std::make_pair(
            std::min(v0,v1),
            std::max(v0,v1)
        ));
    }

    bool contain_facet_edge(index_t facet_index, index_t local_edge) const {
        // local edge k is the one between local vertices k and (k+1)%3
        // see https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
        geo_assert(facet_index < mesh_.facets.nb());
        geo_assert(local_edge < 3);
        index_t v0 = mesh_.facets.vertex(facet_index,local_edge),
                v1 = mesh_.facets.vertex(facet_index,(local_edge+1)%3);
        return feature_edges_.contains(std::make_pair(
            std::min(v0,v1),
            std::max(v0,v1)
        ));
    }

    auto begin()   { return feature_edges_.begin();   }
    auto end()     { return feature_edges_.end();     }
    auto rbegin()  { return feature_edges_.rbegin();  }
    auto rend()    { return feature_edges_.rend();    }
    auto cbegin() const  { return feature_edges_.cbegin();  }
    auto cend() const    { return feature_edges_.cend();    }
    auto crbegin() const { return feature_edges_.crbegin(); }
    auto crend() const   { return feature_edges_.crend();   }

protected:

    std::set<std::pair<index_t,index_t>> feature_edges_; // significant feature edges (with FEATURE_EDGES_MIN_ANGLE as threshold)
    const Mesh& mesh_;
    const MeshExtAdjFacetCorners& adj_facet_corners_;
};

class MeshExt {

public:
    MeshExt(Mesh& mesh) : 
        halfedges(mesh),
        facet_normals(mesh),
        adj_facet_corners(mesh),
        feature_edges(mesh,adj_facet_corners),
        vertices(mesh.vertices),
        edges(mesh.edges),
        facets(mesh.facets),
        facet_corners(mesh.facet_corners),
        cells(mesh.cells),
        cell_corners(mesh.cell_corners),
        cell_facets(mesh.cell_facets),
        mesh_(mesh)
    {};

    void clear() {
        // also clear the underlying GEO::Mesh?
        facet_normals.clear();
        adj_facet_corners.clear();
        feature_edges.clear();
    }

    // Auto-convert to GEO::Mesh

    operator Mesh&() {
        return mesh_;
    }

    operator const Mesh&() const {
        return mesh_;
    }

public:
    MeshHalfedgesExt halfedges; // half-edges interface of the mesh
    MeshExtFacetNormals facet_normals; // per-facet normals
    MeshExtAdjFacetCorners adj_facet_corners; // vertex to facets adjacency
    MeshExtFeatureEdges feature_edges; // significant feature edges
    // Expose the GEO::Mesh interface
    MeshVertices& vertices;
    MeshEdges& edges;
    MeshFacets& facets;
    MeshFacetCornersStore& facet_corners;
    MeshCells& cells;
    MeshCellCornersStore& cell_corners;
    MeshCellFacetsStore& cell_facets;

protected:
    Mesh& mesh_;
};