#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include "LabelingGraph.h"

#define LABEL2STR(label) ((label==0) ? "+X" : ( \
                          (label==1) ? "-X" : ( \
                          (label==2) ? "+Y" : ( \
                          (label==3) ? "-Y" : ( \
                          (label==4) ? "+Z" : ( \
                          (label==5) ? "-Z" :   \
                          "undef."))))))

using namespace GEO;

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name);

bool save_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name);

index_t nearest_label(const vec3& normal);

// similar to https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/flagging_utils.cpp#LL5C17-L5C27 axesMatrix()
const vec3 label2vector[6] = {
    { 1.0, 0.0, 0.0}, // +X
    {-1.0, 0.0, 0.0}, // -X
    { 0.0, 1.0, 0.0}, // +Y
    { 0.0,-1.0, 0.0}, // -Y
    { 0.0, 0.0, 1.0}, // +Z
    { 0.0, 0.0,-1.0}  // -Z
};

/**
 * \brief Compute the naive labeling of a given mesh
 * \details Compute the per-facet nearest-to-normal signed direction +/-{X,Y,Z} of a surface mesh
 * and store it in a facet attribute
 * \param[in,out] mesh A surface triangle mesh
 * \param[in] attribute_name The name of the facet attribute in which the labeling will be stored
 */
void naive_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name);

void graphcut_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, int compactness_coeff = 1, int fidelity_coeff = 1);

// returns {min,max,avg}
std::tuple<double,double,double> compute_per_facet_fidelity(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* labeling_attribute_name, const char* fidelity_attribute_name);

unsigned int remove_surrounded_charts(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg);

unsigned int fix_invalid_boundaries(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg);

unsigned int fix_invalid_corners(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg);

void remove_invalid_charts(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg);

unsigned int move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg);