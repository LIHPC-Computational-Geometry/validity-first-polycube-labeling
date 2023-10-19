// inspired by https://github.com/fprotais/robustPolycube/blob/main/lib/utils/trace.h

#pragma once

#include <geogram/mesh/mesh.h>

#include <string>

#include "LabelingGraph.h"          // for Boundary
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

using namespace GEO;

bool dump_vertex(std::string filename, const Mesh& mesh, index_t vertex_index);

bool dump_all_boundaries(std::string filename, const Mesh& mesh, const CustomMeshHalfedges& mesh_he, const std::vector<Boundary>& boundaries);