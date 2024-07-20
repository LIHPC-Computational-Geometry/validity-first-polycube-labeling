#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include <string>

using namespace GEO;

bool load_labeling(const std::string& filename, const Mesh& mesh, Attribute<index_t>& labeling);

bool save_labeling(const std::string& filename, const Mesh& mesh, Attribute<index_t>& labeling);
