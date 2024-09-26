#pragma once

#include <geogram/basic/numeric.h>          // for GEO::index_t
#include <geogram/basic/vecg.h>             // for GEO::vecng
#include <geogram/mesh/mesh_halfedges.h>    // for MeshHalfedges::Halfedge

#include <DisjointSet.hpp>

#include <utility>  // for std::pair
#include <optional>

#include "geometry_halfedges.h"
#include "geometry_mesh_ext.h"

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

bool move_to_next_halfedge_on_feature_edge(const MeshExt& mesh, MeshHalfedges::Halfedge& H);

void rotate_mesh_according_to_principal_axes(Mesh& mesh);

// Return an outgoing halfedge,
// specifically the one for which H.facet == adj_facets[vertex_index],
// but you should use this function when you don't care about the specific halfedge returned,
// as long as `vertex_index` is at its origin
MeshHalfedges::Halfedge get_an_outgoing_halfedge_of_vertex(const MeshExt& mesh, index_t vertex_index);

MeshHalfedges::Halfedge get_halfedge_between_two_vertices(const MeshExt& mesh, index_t origin_vertex, index_t extremity_vertex);

// the vertex around which halfedges will be tested is the origin vertex of `init_halfedge` 
MeshHalfedges::Halfedge get_most_aligned_halfedge_around_vertex(const MeshExt& mesh, const MeshHalfedges::Halfedge& init_halfedge, const vec3& reference);

void get_adjacent_facets_conditional(const Mesh& mesh, index_t facet_index, index_t which_chart, const std::vector<index_t>& facet2chart, std::set<index_t>& out);

// for given 3D coordinates and a start vertex not too far, find the nearest vertex
index_t get_nearest_vertex_of_coordinates(const MeshExt& mesh, vec3 target_coordinates, index_t start_vertex, size_t max_dist);

index_t nearest_axis_of_edges(const Mesh& mesh, std::initializer_list<MeshHalfedges::Halfedge> edges, std::initializer_list<index_t> forbidden_axes);

// return value between 0 and 1
double dot_product_between_halfedge_and_axis(const Mesh& mesh, MeshHalfedges::Halfedge halfedge, index_t axis);

// return value between 0 and pi
double angle_between_halfedge_and_axis(const Mesh& mesh, MeshHalfedges::Halfedge halfedge, index_t axis);

void trace_path_on_chart(const MeshExt& mesh, const std::vector<index_t>& facet2chart, const std::map<index_t,std::vector<std::pair<index_t,index_t>>>& turning_point_vertices, index_t start_vertex, vec3 direction, std::set<index_t>& facets_at_left, std::set<index_t>& facets_at_right, std::vector<MeshHalfedges::Halfedge>& halfedges);

// TODO move to geometry_mesh_ext.h > MeshExtFacetNormals::average_over(const std::set<index_t>& facets)
vec3 average_facets_normal(const std::vector<vec3>& normals, const std::set<index_t>& facets);

// TODO move to geometry_mesh_ext.h > MeshExtFacetNormals::average_dot_product_over(const std::set<index_t>& facets)
double average_dot_product(const std::vector<vec3>& normals, const std::set<index_t>& facets, vec3 reference);

// TODO move to geometry_mesh_ext.h > MeshExtFacetNormals::average_angle_over(const std::set<index_t>& facets)
double average_angle(const std::vector<vec3>& normals, const std::set<index_t>& facets, vec3 reference);

bool vertex_has_lost_feature_edge_in_neighborhood(const MeshExt& mesh, index_t vertex, MeshHalfedges::Halfedge& outgoing_halfedge_on_feature_edges);

// return last halfedge on the feature edge that is still on the same chart, or the last halfedge before unsuccessful move
MeshHalfedges::Halfedge follow_feature_edge_on_chart(const MeshExt& mesh, MeshHalfedges::Halfedge halfedge, const std::vector<index_t>& facet2chart, std::set<index_t>& facets_at_left, std::set<index_t>& facets_at_right);

bool is_a_facet_to_tilt(const vec3& facet_normal, double sensitivity);

// return nb facets to tilt (size of set)
size_t get_facets_to_tilt(const MeshExt& mesh_ext, std::set<index_t>& facets_to_tilt, double sensitivity);

// return nb groups
template <template<class> class C> // C a container -> usable with std::vector<> or GEO::Attribute<>
index_t group_facets(const Mesh& mesh, const std::set<index_t>& facets_to_tilt, C<index_t>& per_facet_group_index) {
    geo_assert(per_facet_group_index.size() == mesh.facets.nb());
    DisjointSet<index_t> ds(mesh.facets.nb());
    index_t adjacent_facet = index_t(-1);
    std::optional<index_t> a_facet_not_to_tilt = std::nullopt;
    FOR(f,mesh.facets.nb()) {
        FOR(le,3) { // for each local edge
            adjacent_facet = mesh.facets.adjacent(f,le);
            if(facets_to_tilt.contains(f) == facets_to_tilt.contains(adjacent_facet)) {
                // both are in the set or both are out
                ds.mergeSets(f,adjacent_facet);
            }
            if(!facets_to_tilt.contains(f)) {
                a_facet_not_to_tilt = f;
            }
        }
    }
    return ds.getSetsMap(per_facet_group_index.data(),a_facet_not_to_tilt); // return association between facet and group index, and impose facets outside `facets_to_tilt` to have 0 as group index
}

// compute standard deviation of the adjacent facets area
double sd_adjacent_facets_area(const Mesh& mesh, const std::vector<std::vector<index_t>>& adj_facets, index_t vertex_index);

void triangulate_facets(Mesh& M, std::vector<index_t>& triangle_index_to_old_facet_index, std::vector<index_t>& corner_index_to_old_corner_index);

mat3 rotation_matrix(double OX_OY_OZ_angle);

bool vertex_has_a_feature_edge_in_its_ring(const MeshExt& M, index_t vertex_index, MeshHalfedges::Halfedge& found_halfedge);