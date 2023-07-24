#pragma once

#include <geogram/mesh/mesh_halfedges.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

using namespace GEO;
using namespace GEO::Geom;

class CustomMeshHalfedges : public MeshHalfedges {
public:

    CustomMeshHalfedges(Mesh& mesh) : MeshHalfedges(mesh) {}

    // re-expose this method from MeshHalfedges
    void set_use_facet_region(std::string attribute_name) {
        MeshHalfedges::set_use_facet_region(attribute_name);
    }

    // enforce the right function when passing a char array
    void set_use_facet_region(const char* attribute_name) {
        set_use_facet_region(std::string(attribute_name));
    }

    bool move_to_next_around_vertex(Halfedge& H, bool ignore_borders = false) const { // 2nd argument added
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

    bool move_to_prev_around_vertex(Halfedge& H, bool ignore_borders = false) const { // 2nd argument added
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

    void custom_move_to_prev_around_border(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        geo_debug_assert(halfedge_is_border(H));
        index_t count = 0;
        while(move_to_prev_around_vertex(H)) {
            ++count;
            geo_assert(count < 10000);
        }
        move_to_prev_around_facet(H);
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

    index_t opposite_corner(Halfedge& H) {
        move_to_opposite(H);
        index_t corner_index = H.corner;
        move_to_opposite(H);
        return corner_index;
    }
};

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<CustomMeshHalfedges::Halfedge> : ostream_formatter {};