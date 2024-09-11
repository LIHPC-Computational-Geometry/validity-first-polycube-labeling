#pragma once

#include <geogram/basic/geometry.h> // for GEO::mat3

#include <cmath> // for std::pow(), std::abs()

#include "geometry_mesh_ext.h" // for MeshExt
#include "geometry.h" // for rotation_matrix()

#define NAIVE_LABELING_TWEAK_SENSITIVITY 0.1 // min difference between the 2 closest labels before the rotation tweak
#define NAIVE_LABELING_TWEAK_ANGLE 0.1 // angle of rotation of the normal when we cannot choose between 2 or 3 labels

// PolyCut [1, Section 3.1, page 5] and Evocube [2, section 4, p. 7] use a fidelity/compactness ratio of 3 for the first try
// [1] Livesu, Vining, Sheffer, Gregson, Scateni, "Polycut: Monotone graph-cuts for polycube base-complex construction", ACM Trans. on Graphics, 2013
// [2] Dumery, Protais, Mestrallet, Bourcier, Ledoux, "Evocube: a Genetic Labeling Framework for Polycube-Maps", Computer Graphics Forum, 2022
#define DEFAULT_COMPACTNESS 1
#define DEFAULT_FIDELITY	3

// normals pre-processing
#define DEFAULT_SENSITIVITY		  10e-10
#define DEFAULT_ANGLE_OF_ROTATION 0.1

using namespace GEO;

// rotation of NAIVE_LABELING_TWEAK_ANGLE around the x, y and z axes
const mat3 ROTATION_MATRIX = rotation_matrix(NAIVE_LABELING_TWEAK_ANGLE);

// assign randomly one of the six labels to each surface facet (useless, unless to show that not all labelings lead to polycubes)
void random_labeling(const Mesh& mesh, Attribute<index_t>& labeling);

/**
 * \brief Compute the naive labeling of a given mesh
 * \details Compute the per-facet nearest-to-normal signed direction +/-{X,Y,Z} of a surface mesh
 * and store it in a facet attribute
 * \param[in,out] mesh_ext A surface triangle mesh
 * \param[in] attribute_name The name of the facet attribute in which the labeling will be stored
 */
void naive_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling);

/*
 * Like naive_labeling() but triangle normals are considered slightly rotated
 * when 2 or 3 labels have almost the same cost (like areas at 45°)
 * -> avoid labeling fragmentation
 */
void tweaked_naive_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling);

// apply either naive_labeling() or tweaked_naive_labeling() depending on the mesh
void smart_init_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling);

void graphcut_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling, int compactness_coeff = DEFAULT_COMPACTNESS, int fidelity_coeff = DEFAULT_FIDELITY, double sensitivity = DEFAULT_SENSITIVITY, double angle_of_rotation = DEFAULT_ANGLE_OF_ROTATION);
