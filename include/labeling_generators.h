#pragma once

#include <geogram/basic/geometry.h> // for GEO::mat3

#include <cmath> // for std::pow(), std::abs()

#include "geometry_mesh_ext.h" // for MeshExt
#include "geometry.h" // for rotation_matrix()

// PolyCut [1, Section 3.1, page 5] and Evocube [2, section 4, p. 7] use a fidelity/compactness ratio of 3 for the first try
// [1] Livesu, Vining, Sheffer, Gregson, Scateni, "Polycut: Monotone graph-cuts for polycube base-complex construction", ACM Trans. on Graphics, 2013
// [2] Dumery, Protais, Mestrallet, Bourcier, Ledoux, "Evocube: a Genetic Labeling Framework for Polycube-Maps", Computer Graphics Forum, 2022
#define DEFAULT_COMPACTNESS 1
#define DEFAULT_FIDELITY	3

// normals pre-processing
#define DEFAULT_SENSITIVITY		  10e-10 // max difference between the 2 closest label weights for a facet normal rotation
#define DEFAULT_ANGLE_OF_ROTATION 0.05   // amount of the rotation, around the OX, OY and OZ axes

using namespace GEO;

// assign randomly one of the six labels to each surface facet (useless, unless to show that not all labelings lead to polycubes)
void random_labeling(const Mesh& mesh, Attribute<index_t>& labeling);

/**
 * \brief Compute the naive labeling of a given mesh
 * \details Compute the per-facet nearest-to-normal signed direction +/-{X,Y,Z} of a surface mesh
 * and store it in a facet attribute. `sensitivity` controls a facet normals pre-processing.
 * \param[in] mesh_ext A surface triangle mesh
 * \param[in] labeling The facet attribute in which the labeling will be stored
 * \param[in] sensitivity Max difference between the 2 closest label weights for a facet normal rotation. DEFAULT IS OFF (= zero).
 * \param[in] angle_of_rotation Amount of the rotation, around the OX, OY and OZ axes
 */
void naive_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling, double sensitivity = 0.0, double angle_of_rotation = DEFAULT_ANGLE_OF_ROTATION);
// if sensitivity!=0.0, behave like what was tweaked_naive_labeling()

void graphcut_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling, int compactness_coeff = DEFAULT_COMPACTNESS, int fidelity_coeff = DEFAULT_FIDELITY, double sensitivity = DEFAULT_SENSITIVITY, double angle_of_rotation = DEFAULT_ANGLE_OF_ROTATION);
