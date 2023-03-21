#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/basic/vecg.h>         // for GEO::vec3.normalize()
#include <geogram/basic/logger.h>       // for GEO::Logger
#include <geogram/mesh/mesh.h>          // for GEO::Mesh

#include <algorithm>    // for std::min()
#include <ostream>      // for std::endl

#include "hex_mesh.h"

double compute_scaled_jacobian(GEO::Mesh& M) {
    GEO::Attribute<double> SJ(M.cells.attributes(), "SJ");
    double minSJ = 1.0;
    for (GEO::index_t hex_index : M.cells) { // for each hexahedron in mesh M
        
        if (M.cells.type(hex_index) != GEO::MESH_HEX) { // check if the cell is not an hexahedron
            SJ[hex_index] = nanf64("");
            GEO::Logger::out("") << "Found a cell that is not an hexahedron. type = " << M.cells.type(hex_index) << std::endl;
            continue;
        }
        double scaled_jacobian = 1.0;
        for (int hex_corner = 0; hex_corner < 8; hex_corner++) { // for each corner of the current hexahedron
            std::array<GEO::vec3, 4> v;
            for (std::array<GEO::vec3, 4>::size_type i = 0; i < 4; i++) { // for each neighboring corner of the current corner
                
                GEO::index_t vertex_index = M.cells.vertex(hex_index,HEX_CORNER_SPLITTING[hex_corner][i]); // get vertex index of this corner
                // or
                //  GEO::index_t corner_index = M.cells.corner(hex_index,HEX_CORNER_SPLITTING[hex_corner][i]);
                //  GEO::index_t vertex_index = M.cell_corners.vertex(corner_index);
                v[i] = M.vertices.point(vertex_index); // get coordinates of vertex_index
            }
            GEO::vec3 n1 = GEO::normalize(v[1] - v[0]);
            GEO::vec3 n2 = GEO::normalize(v[2] - v[0]);
            GEO::vec3 n3 = GEO::normalize(v[3] - v[0]);
            scaled_jacobian = std::min(scaled_jacobian, dot(n3,cross(n1, n2)));
        }
        SJ[hex_index] = scaled_jacobian; // store value in mesh attribute
        minSJ = std::min(minSJ, scaled_jacobian); // update min Scaled Jacobian
    }
    return minSJ;
}