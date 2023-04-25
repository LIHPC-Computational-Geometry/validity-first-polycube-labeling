#pragma once

#include <geogram/basic/numeric.h>  // for GEO::index_t
#include <geogram/basic/vecg.h>     // for GEO::vecng

#include "CustomMeshHalfedges.h"

typedef vecng<3, Numeric::float64> vec3d; // Geogram only has a float version of vec3

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
    bool operator()(const CustomMeshHalfedges::Halfedge& a, const CustomMeshHalfedges::Halfedge& b) const {
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