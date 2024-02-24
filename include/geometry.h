#pragma once

#include <geogram/basic/numeric.h>          // for GEO::index_t
#include <geogram/basic/vecg.h>             // for GEO::vecng
#include <geogram/mesh/mesh_halfedges.h>    // for MeshHalfedges::Halfedge

#include <utility>  // for std::pair

#include "CustomMeshHalfedges.h"

#define FEATURE_EDGES_MIN_ANGLE	0.5

inline int label2axis(GEO::index_t label) {
    // GEO::index_t is unsigned, so no need to check (0 <= label)
    geo_assert(label <= 5);
    return (int) (label / 2); // integer division
}

const int OTHER_AXIS[3][3] = {{-1, 2, 1},
                              { 2,-1, 0},
                              { 1, 0,-1}};

inline int other_axis(int axis1, int axis2) {
    // axis1, axis2 in {0 = X, 1 = Y, 2 = Z}
    // other_axis in {-1 = undefined, 0 = X, 1 = Y, 2 = Z}
    geo_assert((0 <= axis1) && (axis1 <= 2));
    geo_assert((0 <= axis2) && (axis2 <= 2));
    return OTHER_AXIS[axis1][axis2];
}

// in order to have halfedges as key of std::map, we have to define comparison between halfedge

struct HalfedgeCompare {
    bool operator()(const MeshHalfedges::Halfedge& a, const MeshHalfedges::Halfedge& b) const {
        // compare facet, then corner 
        if( a.facet < b.facet ) return true;
        else if( a.facet == b.facet ) return (a.corner < b.corner);
        else return false;
    }
};

// in order to have vecng as key of std::map, we have to define comparison between vecng

template <index_t DIM, typename T>
struct VecngCompare {
    bool operator()(const GEO::vecng<DIM,T>& a, const GEO::vecng<DIM,T>& b) const {
        for(index_t d = 0; d < DIM; ++d) {
            if(a[d] < b[d]) {
                return true;
            }
            else if (a[d] > b[d]) {
                return false;
            }
            // else : compare other dimensions
        }
        return false; // case a == b
    }
};

// to compare 2 vecng in tests, we have to define equality between vecng

// operator must be in the same namespace as the operand type
// https://itecnote.com/tecnote/c-how-to-get-a-custom-operator-to-work-with-google-test/ 
namespace GEO {
    template <index_t DIM, typename T>
    bool operator==(const vecng<DIM,T>& a, const vecng<DIM,T>& b) {
        for(index_t d = 0; d < DIM; ++d) {
            if(a[d] != b[d]) {
                return false;
            }
            // else : compare other dimensions
        }
        return true;
    }
}

// map keys = surface to explore with all the facet indices
// map values = 0 if distance of 0, unsigned int(-1) if to be computed in accordance
void per_facet_distance(const Mesh& mesh, std::map<index_t,unsigned int>& distance);

bool facet_normals_are_inward(Mesh& mesh); // a mutable Mesh is required by MeshHalfedges

// TODO use geogram's invert_normals() ? in mesh/mesh_preprocessing.h
void flip_facet_normals(Mesh& mesh);

// shift a mesh so that 0,0,0 is at the center of the bounding box
// if normalize==true, scale the mesh so that the bounding box is in [-1,1]^3
void center_mesh(Mesh& mesh, bool normalize);

void compute_adjacent_facets_of_vertices(const Mesh& mesh, std::vector<std::vector<index_t>>& adj_facets);

struct AdjacentFacetOfVertex {
    index_t facet_index;
    index_t local_vertex;
};

// adj_facet_corners has mesh.vertices.nb() elements, each of them is a list of adjacent facet, with 2 components : the facet index, and the local vertex index for this facet
void compute_adjacent_facets_of_vertices(const Mesh& mesh, std::vector<std::vector<AdjacentFacetOfVertex>>& adj_facet_corners);

// /!\ if adj_facets not empty, adjacency will not be recomputed
void remove_feature_edges_with_low_dihedral_angle(Mesh& mesh, std::vector<std::vector<index_t>>& adj_facets);

// Move feature edges from mesh.edges to set of pair of vertices
// Why? 
// Because when working on the labeling graph (charts, boundaries & corners), we use oriented edges
// and we need to know if a given oriented edge is on a feature edge
// When feature edges are stored in mesh edges, see https://github.com/BrunoLevy/geogram/wiki/Mesh#mesh-edges
// finding if they contains (v0,v1) is expensive, we have to check all edges
// Instead we can use a set of pair of vertices, the pair being sorted by ascending index,
// to quickly check if (v0,v1) is a feature edge with `feature_edges.contains(min(v0,v1),max(v0,v1))`
void transfer_feature_edges(Mesh& mesh, std::set<std::pair<index_t,index_t>>& feature_edges);

bool halfedge_is_on_feature_edge(const Mesh& mesh, const MeshHalfedges::Halfedge& H, const std::set<std::pair<index_t,index_t>>& feature_edges);

void rotate_mesh_according_to_principal_axes(Mesh& mesh);

// Return an outgoing halfedge,
// specifically the one for which H.facet == adj_facets[vertex_index],
// but you should use this function when you don't care about the specific halfedge returned,
// as long as `vertex_index` is at its origin
MeshHalfedges::Halfedge get_an_outgoing_halfedge_of_vertex(const Mesh& mesh, const std::vector<std::vector<index_t>>& adj_facets, index_t vertex_index);

MeshHalfedges::Halfedge get_halfedge_between_two_vertices(const CustomMeshHalfedges& mesh_he, const std::vector<std::vector<index_t>>& adj_facets, index_t origin_vertex, index_t extremity_vertex);

// the vertex around which halfedges will be tested is the origin vertex of `init_halfedge` 
MeshHalfedges::Halfedge get_most_aligned_halfedge_around_vertex(const MeshHalfedges::Halfedge& init_halfedge, const CustomMeshHalfedges& mesh_he, const vec3& reference);

void get_adjacent_facets_conditional(const Mesh& mesh, index_t facet_index, index_t which_chart, const std::vector<index_t>& facet2chart, std::set<index_t>& out);

// for given 3D coodinates and a start vertex not too far, find the nearest vertex
index_t get_nearest_vertex_of_coordinates(const CustomMeshHalfedges& mesh_he, const std::vector<std::vector<index_t>>& adj_facets, vec3 target_coordinates, index_t start_vertex, size_t max_dist);