#pragma once

#include <geogram/basic/geometry.h> // for GEO::mat3

#include <cmath> // for std::pow(), std::abs()

#include <geometry_mesh_ext.h> // for MeshExt

#define NAIVE_LABELING_TWEAK_SENSITIVITY 0.1 // min difference between the 2 closest labels before the rotation tweak
#define NAIVE_LABELING_TWEAK_ANGLE 0.1 // angle of rotation of the normal when we cannot choose between 2 or 3 labels

using namespace GEO;

// https://en.wikipedia.org/wiki/Rotation_matrix#Basic_3D_rotations
// rotation of NAIVE_LABELING_TWEAK_ANGLE around the x, y and z axes
const double COS_TILT_ANGLE = cos(NAIVE_LABELING_TWEAK_ANGLE);
const double SIN_TILT_ANGLE = sin(NAIVE_LABELING_TWEAK_ANGLE);
const double COS_SQUARED_TILT_ANGLE = COS_TILT_ANGLE*COS_TILT_ANGLE;
const double SIN_SQUARED_TILT_ANGLE = SIN_TILT_ANGLE*SIN_TILT_ANGLE;
const double SIN_BY_COS_TILT_ANGLE = SIN_TILT_ANGLE*COS_TILT_ANGLE;
const mat3 rotation({
    {
        COS_SQUARED_TILT_ANGLE,                                 // <=> cos*cos              @ (0,0)
        SIN_BY_COS_TILT_ANGLE*(SIN_TILT_ANGLE-1),               // <=> sin*sin*cos-cos*sin  @ (0,1)
        SIN_TILT_ANGLE*(COS_SQUARED_TILT_ANGLE+SIN_TILT_ANGLE)  // <=> cos*sin*cos+sin*sin  @ (0,2)
    },
    {
        SIN_BY_COS_TILT_ANGLE,                                          // <=> cos*sin              @ (1,0)
        SIN_SQUARED_TILT_ANGLE*SIN_TILT_ANGLE+COS_SQUARED_TILT_ANGLE,   // <=> sin*sin*sin+cos*cos  @ (1,1)
        SIN_BY_COS_TILT_ANGLE*(SIN_TILT_ANGLE-1)                        // <=> cos*sin*sin-sin*cos  @ (1,2)
    },
    {
        -SIN_TILT_ANGLE,        // <=> -sin     @ (2,0)
        SIN_BY_COS_TILT_ANGLE,  // <=> sin*cos  @ (2,1)
        COS_SQUARED_TILT_ANGLE  // <=> cos*cos  @ (2,2)
    }
});

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

void graphcut_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling, int compactness_coeff = 1, int fidelity_coeff = 1);
