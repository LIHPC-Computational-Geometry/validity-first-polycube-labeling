#include <geogram/mesh/mesh.h>  // for GEO::Mesh
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/attributes.h>

#include <fmt/core.h>

#include <vector>
#include <string>

using namespace GEO;

int main(int argc, const char** argv) {

    if (argc<3) {
        fmt::println("Missing arguments");
        fmt::println("./extract_surface input_tetra.mesh output_surface.obj");
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

    mesh_save(triangle_mesh,argv[2]);

    // TODO write map between surface facets & cells

    return 0;
}