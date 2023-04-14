#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

using namespace GEO;

/**
 * \brief Compute the naive labeling of a given mesh
 * \details Compute the per-facet nearest-to-normal signed direction +/-{X,Y,Z} of a surface mesh
 * and store it in a facet attribute
 * \param[in,out] mesh A surface triangle mesh
 * \param[in] attribute_name The name of the facet attribute in which the labeling will be stored
 */
void naive_labeling(GEO::Mesh& mesh, const char* attribute_name);