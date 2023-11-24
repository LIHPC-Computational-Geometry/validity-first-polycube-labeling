// Extended halfedges interface of Geogram
// 
// CustomMeshHalfedges is the interface itself
// CustomMeshHalfedges::Halfedge is an oriented edge, which stores a facet index and a facet corner index
// CustomMeshHalfedges::Halfedge.facet is the facet at the right of the halfedge
// CustomMeshHalfedges::Halfedge.corner is the corner of the facet that is along the halfedge origin
// 
//      halfedge
//   o ==========> o
//    \/corner    /
//     \         /
//      \ facet /
//       \     /
//        \   /
//         \ /
//          o
// 
// Around a vertex:
// - move_to_next_around_vertex() moves counterclockwise
// - move_to_prev_around_vertex() moves clockwise
// so counterclockwise, the halfedge is ahead of the facet
// and clockwise, the facet is ahead of the halfedge
//
// see https://github.com/BrunoLevy/geogram/pull/116

#pragma once

#include <geogram/mesh/mesh_halfedges.h>
#include <geogram/mesh/mesh_geometry.h> // for mesh_facet_normal()

#include <fmt/core.h>
#include <fmt/ostream.h>

using namespace GEO;
using namespace GEO::Geom;

namespace GEO {

    class CustomMeshHalfedges : public MeshHalfedges {
    public:

        CustomMeshHalfedges(Mesh& mesh) : MeshHalfedges(mesh) {}

        // re-expose this method from MeshHalfedges
        void set_use_facet_region(std::string attribute_name) {
            MeshHalfedges::set_use_facet_region(attribute_name);
        }

        // enforce the right function when passing a char array
        // fixed since https://github.com/BrunoLevy/geogram/commit/ddacc7013124fbb1b8e4514f8d9613bf48d2d95e
        void set_use_facet_region(const char* attribute_name) {
            set_use_facet_region(std::string(attribute_name));
        }

        // re-expose this method from MeshHalfedges
        void move_to_next_around_facet(Halfedge& H) const {
            MeshHalfedges::move_to_next_around_facet(H);
        }

        // re-expose this method from MeshHalfedges
        void move_to_prev_around_facet(Halfedge& H) const {
            MeshHalfedges::move_to_prev_around_facet(H);
        }

        inline void move_clockwise_around_facet(Halfedge& H) const {
            // around facets, next is clockwise and prev is counterclockwise
            // -> call move_to_next_around_facet()
            move_to_next_around_facet(H);
        }

        inline void move_counterclockwise_around_facet(Halfedge& H) const {
            // around facets, next is clockwise and prev is counterclockwise
            // -> call move_to_prev_around_facet()
            geo_debug_assert(halfedge_is_valid(H));
            move_to_prev_around_facet(H);
        }

        // re-expose this method from MeshHalfedges and add 2nd argument
        bool move_to_next_around_vertex(Halfedge& H, bool ignore_borders = false) const {
            geo_debug_assert(halfedge_is_valid(H));
            index_t v = mesh_.facet_corners.vertex(H.corner);
            index_t f = mesh_.facet_corners.adjacent_facet(H.corner);
            if(f == NO_FACET) {
                return false;
            }
            if(
                ignore_borders == false &&
                facet_region_.is_bound() &&
                facet_region_[H.facet] != facet_region_[f]
            ) {
                return false;
            }
            for(index_t c: mesh_.facets.corners(f)) {
                index_t pc = mesh_.facets.prev_corner_around_facet(f, c);
                if(
                    mesh_.facet_corners.vertex(c) == v &&
                    mesh_.facet_corners.adjacent_facet(pc) == H.facet
                ) {
                    H.corner = c;
                    H.facet = f;
                    return true;
                }
            }
            geo_assert_not_reached;
        }

        // re-expose this method from MeshHalfedges and add 2nd argument
        bool move_to_prev_around_vertex(Halfedge& H, bool ignore_borders = false) const {
            geo_debug_assert(halfedge_is_valid(H));
            index_t v = mesh_.facet_corners.vertex(H.corner);
            index_t pc = mesh_.facets.prev_corner_around_facet(H.facet, H.corner);
            index_t f = mesh_.facet_corners.adjacent_facet(pc);
            if(f == NO_FACET) {
                return false;
            }
            if(
                ignore_borders == false &&
                facet_region_.is_bound() &&
                facet_region_[H.facet] != facet_region_[f]
            ) {
                return false;
            }
            for(index_t c: mesh_.facets.corners(f)) {
                if(
                    mesh_.facet_corners.vertex(c) == v &&
                    mesh_.facet_corners.adjacent_facet(c) == H.facet
                ) {
                    H.corner = c;
                    H.facet = f;
                    return true;
                }
            }
            geo_assert_not_reached;
        }

        inline bool move_clockwise_around_vertex(Halfedge& H, bool ignore_borders = false) const {
            // around vertices, next is counterclockwise and prev is clockwise
            // -> call move_to_prev_around_vertex()
            return move_to_prev_around_vertex(H,ignore_borders);
        }

        inline bool move_counterclockwise_around_vertex(Halfedge& H, bool ignore_borders = false) const {
            // around vertices, next is counterclockwise and prev is clockwise
            // -> call move_to_next_around_vertex()
            return move_to_next_around_vertex(H,ignore_borders);
        }

        // re-expose this method from MeshHalfedges
        void move_to_next_around_border(Halfedge& H) const {
            MeshHalfedges::move_to_next_around_border(H);
        }

        // re-expose this method from MeshHalfedges
        void move_to_prev_around_border(Halfedge& H) const {
            MeshHalfedges::move_to_prev_around_border(H);
        }

        // TODO remove. the right one is the one from MeshHalfedges (vanilla Geogram)
        // -> equivalent to move_to_opposite() !!!
        void custom_move_to_next_around_border(Halfedge& H) const {
            geo_debug_assert(halfedge_is_valid(H));
            geo_debug_assert(halfedge_is_border(H));
            index_t count = 0;
            do {
                ++count;
                geo_assert(count < 10000);
                move_to_next_around_vertex(H,true);
            } while(!halfedge_is_border(H));
        }

        // TODO remove. the right one is the one from MeshHalfedges (vanilla Geogram)
        // -> the new halfedge is not even on a border !!!
        void custom_move_to_prev_around_border(Halfedge& H) const {
            geo_debug_assert(halfedge_is_valid(H));
            geo_debug_assert(halfedge_is_border(H));
            index_t count = 0;
            while(move_to_prev_around_vertex(H)) {
                ++count;
                geo_assert(count < 10000);
            }
        }

        bool is_on_lower_than_180_degrees_edge(Halfedge& H) const {
            // define the plane passing through three points of H.facet
            index_t vertex0 = mesh_.facet_corners.vertex(mesh_.facets.corner(H.facet,0));
            index_t vertex1 = mesh_.facet_corners.vertex(mesh_.facets.corner(H.facet,1));
            index_t vertex2 = mesh_.facet_corners.vertex(mesh_.facets.corner(H.facet,2));
            Plane plane(mesh_.vertices.point(vertex0),mesh_.vertices.point(vertex1),mesh_.vertices.point(vertex2));
            // change direction, so that H.facet is the facet at the other side of the edge
            move_to_opposite(H);
            // move so that the 3rd vertex of this facet (the one that is not on the original halfedge) is on the tip of the halhedge
            move_to_next_around_facet(H);
            vec3 tip = halfedge_vertex_to(mesh_,H); // get the coordinates
            // move back to the initial halfedge
            move_to_next_around_facet(H);
            move_to_next_around_facet(H);
            move_to_opposite(H);
            // evaluate the plane equation at tip, return the sign
            return ((plane.a * tip.x) + (plane.b * tip.y) + (plane.c * tip.z) + plane.d < 0);
        }

        vec3 normal(Halfedge& H) const {
            vec3 normal = Geom::mesh_facet_normal(mesh_,H.facet);
            // change direction, so that H.facet is the facet at the other side of the edge
            move_to_opposite(H);
            normal += Geom::mesh_facet_normal(mesh_,H.facet);
            // move back to the initial halfedge
            move_to_opposite(H);
            return normal;
        }
    };

    namespace Geom {

        inline index_t halfedge_vertex_index_from(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            return M.facet_corners.vertex(H.corner);
        }

        inline index_t halfedge_vertex_index_to(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            index_t c = M.facets.next_corner_around_facet(H.facet, H.corner);
            return M.facet_corners.vertex(c);
        }

        inline index_t halfedge_facet_left(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            return M.facet_corners.adjacent_facet(H.corner); // see H.facet new value in move_to_opposite()
        }

        inline index_t halfedge_facet_right(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            geo_argused(M);
            return H.facet;
        }

        inline index_t halfedge_bottom_left_corner(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            index_t left_facet = halfedge_facet_left(M,H);
            index_t corner = index_t(-1);
            if(left_facet == NO_FACET) {
                return index_t(-1);
            }
            FOR(lv,3) {
                corner = M.facets.corner(left_facet,lv);
                if(M.facet_corners.vertex(corner) == halfedge_vertex_index_from(M,H)) {
                    return corner;
                }
            }
            geo_assert_not_reached;
        }

        inline index_t halfedge_bottom_right_corner(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            geo_argused(M);
            return (H.corner == NO_CORNER) ? index_t(-1) : H.corner;
            // // Should be equivalent to:
            // index_t right_facet = halfedge_facet_right(M,H);
            // index_t corner = index_t(-1);
            // if(right_facet == NO_FACET) {
            //     return index_t(-1);
            // }
            // FOR(lv,3) {
            //     corner = M.facets.corner(right_facet,lv);
            //     if(M.facet_corners.vertex(corner) == halfedge_vertex_index_from(M,H)) {
            //         return corner;
            //     }
            // }
            // geo_assert_not_reached;
        }

        inline index_t halfedge_top_left_corner(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            index_t left_facet = halfedge_facet_left(M,H);
            index_t corner = index_t(-1);
            if(left_facet == NO_FACET) {
                return index_t(-1);
            }
            FOR(lv,3) {
                corner = M.facets.corner(left_facet,lv);
                if(M.facet_corners.vertex(corner) == halfedge_vertex_index_to(M,H)) {
                    return corner;
                }
            }
            geo_assert_not_reached;
        }

        inline index_t halfedge_top_right_corner(
            const Mesh& M, const MeshHalfedges::Halfedge& H
        ) {
            index_t right_facet = halfedge_facet_right(M,H);
            index_t corner = index_t(-1);
            if(right_facet == NO_FACET) {
                return index_t(-1);
            }
            FOR(lv,3) {
                corner = M.facets.corner(right_facet,lv);
                if(M.facet_corners.vertex(corner) == halfedge_vertex_index_to(M,H)) {
                    return corner;
                }
            }
            geo_assert_not_reached;
        }
    }

}

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<CustomMeshHalfedges::Halfedge> : ostream_formatter {};