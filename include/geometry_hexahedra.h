#pragma once

#include <geogram/basic/numeric.h> // for GEO::index_t
#include <geogram/mesh/mesh.h>     // for GEO::Mesh

#include "stats.h"

// Ordering of hexahedron corners
//
//   UltiMaille convention
//         6-------7
//        /|      /|
//       / |     / |
//      4-------5  |
//      |  2----|--3
//      | /     | /
//      |/      |/
//      0-------1
//
//   MEDIT convention (.mesh files)
//        5-------6
//       /|      /|
//      / |     / |
//     1-------2  |
//     |  4----|--7
//     | /     | /
//     |/      |/
//     0-------3
//
//   Geogram convention
//        4-------6
//       /|      /|
//      / |     / |
//     0-------2  |
//     |  5----|--7
//     | /     | /
//     |/      |/
//     1-------3
constexpr GEO::index_t HEX_CORNER_SPLITTING[8][4] = {
    // corner index, then its 3 ordered neighbors
    {0,4,2,1},
    {1,3,5,0},
    {2,0,6,3},
    {3,7,1,2},
    {4,6,0,5},
    {5,1,7,4},
    {6,2,4,7},
    {7,5,3,6},
};

// compute the Scaled Jacobian of each cell
// in a cell attribute 'SJ' of M
void compute_scaled_jacobian(GEO::Mesh& M, IncrementalStats& stats);