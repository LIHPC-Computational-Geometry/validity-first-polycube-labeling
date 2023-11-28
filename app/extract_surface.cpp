// Rewriting of https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/tet_boundary.cpp
// using Geogram instead of libigl.
// Idea : Because polycube labeling is only involving the surface of a tetrahedral mesh,
// we extract the surface of the volume mesh, working on the lighter, triangle mesh.
// The (surface) labeling is a text file of nb_triangles lines, of values in {0 = +X, 1 = -X, 2 = +Y, 3 = -Y, 4 = +Z, 5 = -Z}.
// The volume mesh is used later in the pipeline to extract a hexahedral mesh from a polycube (parametrization & quantization).
// 
// But : parametrization algorithms https://github.com/fprotais/polycube_withHexEx and https://github.com/fprotais/robustPolycube
// want a volume labeling/flagging, that is, instead of associating labels to surface triangles, they need labels associated to all cell facets (4 per tetrahedra).
// For inner facets, -1 value is used, for "no label".
// 
// So we need a map between surface triangle and cell facets, to generate, when we need it, the volume labeling from the surface one:
// 1. Read the surface labeling and the map between surface and volume meshes
// 2. Fill the volume labeling (a vector of 4 * nb_cells) with -1 = no label
// 3. For each surface triangle, overwrite label of cell facet (4 * cell + local_facet) with surface_labeling[surface_triangle]
// 
// The map between surface and volume mesh must give, for each surface triangle (1) which cell (2) which of the 4 facets of this cell
// We can write (4*cell + cell_facet) for each line (= surface triangle), so we can get both with integer division or modulo by 4.
// For preallocation, it is convenient to have the number of triangles and the number of tetrehedra, so they are written in a header (2 lines)
//
// Format of surface map ("surface_map.txt" in the help message):
//   line 1      (header)                  | 158110 triangles
//   line 2      (header)                  | 2015758 tetrahedra
//   line 3      = surface triangle 0      | 631
//   line 4      = surface triangle 1      | 2815
//   line 5      = surface triangle 2      | 2935
//         ...          ...                |  ...
//   line 158112 = surface triangle 158109 | 8062990
// 
// In this example,
//   surface triangle 0 is the (631 % 4 =) 3rd facet of cell (631 // 4 =) 157
//   surface triangle 1 is the (2815 % 4 =) 3rd facet of cell (2815 // 4 =) 703
//   surface triangle 2 is the (2935 % 4 =) 3rd facet of cell (2935 // 4 =) 733
//   surface triangle 158109 is the (8062990 % 4 =) 2nd facet of cell (8062990 // 4 =) 2015747
//
// Issue : Geogram can compute the surface mesh with mesh.cells.compute_borders(),
// but only gives us the corresponding cell, without the cell facet.
// So we need to recover which cell facet it was, by comparing surface triangle vertices with vertices of each cell facet.

#include <geogram/mesh/mesh.h>  // for GEO::Mesh
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/attributes.h>

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>

#include <vector>
#include <string>
#include <set>

#include "geometry.h"

using namespace GEO;

int main(int argc, const char** argv) {

    if (argc<3) {
        fmt::println("Missing arguments");
        fmt::println("./extract_surface input_tetra.mesh output_surface.obj surface_map.txt");
        return 1;
    }

    GEO::initialize();

    std::string input = argv[1];

    Mesh tetra_mesh, triangle_mesh;
    MeshIOFlags flags;
    flags.set_element(MESH_CELLS);
    geo_assert(FileSystem::is_file(input));

    if(!GEO::mesh_load(input,tetra_mesh,flags)) {
        fmt::println("Unable to open {}",input);
        return 1;
    }
    geo_assert(tetra_mesh.cells.nb() != 0); // must be a volumetric mesh
    geo_assert(tetra_mesh.cells.are_simplices()); // must be a tetrahedral mesh
    index_t nb_tetra = tetra_mesh.cells.nb();

    triangle_mesh.copy(tetra_mesh);
    Attribute<index_t> surface_map(triangle_mesh.facets.attributes(),"cell");
    triangle_mesh.facets.clear();
    triangle_mesh.cells.compute_borders(surface_map);
    triangle_mesh.cells.clear();
    triangle_mesh.vertices.remove_isolated();

    index_t nb_facets = triangle_mesh.facets.nb();
    geo_assert(surface_map.size() == nb_facets);

    // ensure facet normals are outwards
    if(facet_normals_are_inwards(triangle_mesh)) {
        flip_facet_normals(triangle_mesh);
        fmt::println(Logger::out("normals dir."),"Facet normals were inwards and are now outwards"); Logger::out("normals dir.").flush();
    }

    mesh_save(triangle_mesh,argv[2]);

    // Write map between surface facets & cells
    // For each surface facet, we want to store on which cell it was, and which of the 4 facets of this cell it is
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/tet_boundary.cpp#L33
    // https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/tet_boundary.cpp#L43
    // Same tetrahedra facets ordering as https://github.com/fprotais/polycube_withHexEx and https://github.com/fprotais/robustPolycube
    //   cell facet 0 -> vertices {1,2,3}
    //   cell facet 1 -> vertices {0,3,2}
    //   cell facet 2 -> vertices {0,1,3}
    //   cell facet 3 -> vertices {0,2,1}
    // which is the facets ordering of https://github.com/ssloy/ultimaille
    // Geogram seems to follow the same ordering for tetra_mesh.cells.corner(cell,_)...
    fmt::println(Logger::out("I/O"),"Writing {}",argv[3]); Logger::out("I/O").flush();
    auto out = fmt::output_file(argv[3]);
    // write 2 lines of header : number of triangles and number of tetrahedra
    out.print("{} triangles\n",nb_facets);
    out.print("{} tetrahedra\n",nb_tetra);
    index_t cell = 0, v0 = 0, v1 = 0, v2 = 0, v3 = 0;
    std::set<index_t> facet_vertices;
    FOR(facet,nb_facets) { // for each surface facet
        // get the 3 vertices around this facet
        facet_vertices = {
            triangle_mesh.facets.vertex(facet,0),
            triangle_mesh.facets.vertex(facet,1),
            triangle_mesh.facets.vertex(facet,2)
        };
        // retreive on which cell was this facet
        cell = surface_map[facet];
        // get the 4 vertices of this cell
        v0 = tetra_mesh.cell_corners.vertex(tetra_mesh.cells.corner(cell,0));
        v1 = tetra_mesh.cell_corners.vertex(tetra_mesh.cells.corner(cell,1));
        v2 = tetra_mesh.cell_corners.vertex(tetra_mesh.cells.corner(cell,2));
        v3 = tetra_mesh.cell_corners.vertex(tetra_mesh.cells.corner(cell,3));
        // now we need to know which of the 4 cell facets is the current facet
        if (facet_vertices == std::set<index_t>({v1,v2,v3})) {
            out.print("{}\n",4*cell+0); // current facet is the cell facet n°0
        }
        else if (facet_vertices == std::set<index_t>({v0,v3,v2})) {
            out.print("{}\n",4*cell+1); // current facet is the cell facet n°1
        }
        else if (facet_vertices == std::set<index_t>({v0,v1,v3})) {
            out.print("{}\n",4*cell+2); // current facet is the cell facet n°2
        }
        else if (facet_vertices == std::set<index_t>({v0,v2,v1})) {
            out.print("{}\n",4*cell+3); // current facet is the cell facet n°3
        }
        else {
            geo_assert_not_reached;
        }
    }

    return 0;
}