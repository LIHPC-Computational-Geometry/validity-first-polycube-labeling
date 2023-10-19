#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <string>

#include "dump_mesh.h"

bool dump_vertex(std::string filename, const Mesh& mesh, index_t vertex_index) {
    Mesh out;
    out.vertices.create_vertices(1);
    out.vertices.point(0) = mesh.vertices.point(vertex_index);
    return mesh_save(out,filename + ".geogram");
}