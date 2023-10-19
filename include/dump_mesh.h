// inspired by https://github.com/fprotais/robustPolycube/blob/main/lib/utils/trace.h

#pragma once

#include <geogram/mesh/mesh.h>

#include <string>

using namespace GEO;

bool dump_vertex(std::string filename, const Mesh& mesh, index_t vertex_index);