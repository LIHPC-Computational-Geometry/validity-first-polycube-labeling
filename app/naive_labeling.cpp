// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/naive_labeling.cpp
// using Geogram instead of libigl.
// Write only the surface labeling. The volume one can be generated later with volume_labeling.cpp

#include <geogram/mesh/mesh.h>          // for Mesh
#include <geogram/mesh/mesh_io.h>       // for mesh_load()
#include <geogram/basic/file_system.h>  // for FileSystem::is_file()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "labeling.h"
#include "labeling_generators.h"
#include "labeling_io.h"

using namespace GEO;

int main(int argc, const char** argv) {

    if (argc<2) {
        fmt::println("Missing arguments");
        fmt::println("./volume_labeling triangle_mesh.obj labeling.txt");
        return 1;
    }

    initialize();

    geo_assert(FileSystem::is_file(argv[1]));

    Mesh triangle_mesh;
    if(!mesh_load(argv[1],triangle_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to open {}",argv[1]); Logger::err("I/O").flush();
        return 1;
    }
    geo_assert(triangle_mesh.cells.nb()==0); // must be a surface only mesh
    geo_assert(triangle_mesh.facets.are_simplices()); // must be a triangle mesh

    MeshExt mesh_ext(triangle_mesh); // will compute facet normals

    Attribute<index_t> labeling(triangle_mesh.facets.attributes(),LABELING_ATTRIBUTE_NAME);
    naive_labeling(mesh_ext,labeling);

    fmt::println(Logger::out("I/O"),"Writing {}",argv[2]); Logger::out("I/O").flush();
    save_labeling(argv[2],triangle_mesh,labeling);

    return 0;
}