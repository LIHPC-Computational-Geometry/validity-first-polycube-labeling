#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include <string>

using namespace GEO;

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name);

bool save_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name);
