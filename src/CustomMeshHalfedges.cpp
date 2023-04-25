// based on ext/geogram/src/lib/geogram/mesh/mesh_halfedges.h
//
// # Changelog - https://keepachangelog.com
//
// ## Added
//
// - using namespace GEO
// - inclusion of CustomMeshHalfedges.h
// - definition of custom_move_to_next_around_border()
//
// ## Changed
//
// - move_to_next_around_vertex() and move_to_prev_around_vertex() are independant of facet_region_ if ignore_borders==true
//
// ## Removed
//
// - namespace GEO wrapping the code
// - inclusion of <geogram/mesh/mesh_halfedges.h>

#include "CustomMeshHalfedges.h"

using namespace GEO;

    bool CustomMeshHalfedges::move_to_next_around_vertex(Halfedge& H, bool ignore_borders) const {
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

    bool CustomMeshHalfedges::move_to_prev_around_vertex(Halfedge& H, bool ignore_borders) const {
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

    void CustomMeshHalfedges::move_to_next_around_border(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        geo_debug_assert(halfedge_is_border(H));
        move_to_next_around_facet(H);
        index_t count = 0;
        while(move_to_next_around_vertex(H)) {
            ++count;
            geo_assert(count < 10000);
        }
    }

    void CustomMeshHalfedges::custom_move_to_next_around_border(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        geo_debug_assert(halfedge_is_border(H));
        index_t count = 0;
        do {
            ++count;
            geo_assert(count < 10000);
            move_to_next_around_vertex(H,true);
        } while(!halfedge_is_border(H));
    }

    void CustomMeshHalfedges::move_to_prev_around_border(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        geo_debug_assert(halfedge_is_border(H));
        index_t count = 0;
        while(move_to_prev_around_vertex(H)) {
            ++count;
            geo_assert(count < 10000);
        }
        move_to_prev_around_facet(H);        
    }
    
    void CustomMeshHalfedges::move_to_opposite(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        index_t v = mesh_.facet_corners.vertex(
            mesh_.facets.next_corner_around_facet(H.facet, H.corner)
        );
        index_t f = mesh_.facet_corners.adjacent_facet(H.corner);
        geo_assert(f != NO_FACET);
        for(index_t c: mesh_.facets.corners(f)) {
            if(mesh_.facet_corners.vertex(c) == v) {
                H.facet = f;
                H.corner = c;
                return;
            }
        }
        geo_assert_not_reached;
    }
