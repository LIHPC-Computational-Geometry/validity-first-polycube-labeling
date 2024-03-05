#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include "LabelingGraph.h"
#include "basic_stats.h"

#define LABEL2STR(label) ((label==0) ? "+X" : ( \
                          (label==1) ? "-X" : ( \
                          (label==2) ? "+Y" : ( \
                          (label==3) ? "-Y" : ( \
                          (label==4) ? "+Z" : ( \
                          (label==5) ? "-Z" :   \
                          "undef."))))))

#define NAIVE_LABELING_TWEAK_SENSITIVITY 0.1 // min difference between the 2 closest labels before the rotation tweak
#define NAIVE_LABELING_TWEAK_ANGLE 0.1 // angle of rotation of the normal when we cannot choose between 2 or 3 labels
const double COS_TILT_ANGLE = cos(NAIVE_LABELING_TWEAK_ANGLE);
const double SIN_TILT_ANGLE = sin(NAIVE_LABELING_TWEAK_ANGLE);
const double COS_SQUARED_TILT_ANGLE = COS_TILT_ANGLE*COS_TILT_ANGLE;
const double SIN_SQUARED_TILT_ANGLE = SIN_TILT_ANGLE*SIN_TILT_ANGLE;
const double SIN_BY_COS_TILT_ANGLE = SIN_TILT_ANGLE*COS_TILT_ANGLE;

using namespace GEO;

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name);

bool save_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name);

index_t nearest_label(const vec3& normal);

index_t nearest_axis(const vec3& vector);

bool is_better_label(const vec3& facet_normal, index_t current_label, index_t new_label);

// similar to https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/flagging_utils.cpp#LL5C17-L5C27 axesMatrix()
const vec3 label2vector[6] = {
    { 1.0, 0.0, 0.0}, // +X
    {-1.0, 0.0, 0.0}, // -X
    { 0.0, 1.0, 0.0}, // +Y
    { 0.0,-1.0, 0.0}, // -Y
    { 0.0, 0.0, 1.0}, // +Z
    { 0.0, 0.0,-1.0}  // -Z
};

bool are_orthogonal_labels(index_t label1, index_t label2);

/**
 * \brief Compute the naive labeling of a given mesh
 * \details Compute the per-facet nearest-to-normal signed direction +/-{X,Y,Z} of a surface mesh
 * and store it in a facet attribute
 * \param[in,out] mesh A surface triangle mesh
 * \param[in] attribute_name The name of the facet attribute in which the labeling will be stored
 */
void naive_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name);

/*
 * Like naive_labeling() but triangle normals are considered slightly rotated
 * when 2 or 3 labels have almost the same cost (like areas at 45°)
 * -> avoid labeling fragmentation
 */
void tweaked_naive_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name);

void graphcut_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, int compactness_coeff = 1, int fidelity_coeff = 1);

void compute_per_facet_fidelity(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* labeling_attribute_name, const char* fidelity_attribute_name, BasicStats& stats);

unsigned int remove_surrounded_charts(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg);

unsigned int fix_invalid_boundaries(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, const std::vector<vec3>& facet_normals);

unsigned int fix_invalid_corners(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg);

// return true if all invalid charts are locked, ie we cannot remove them because they are surrounded by feature edges
bool remove_invalid_charts(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg);

void remove_charts_around_invalid_boundaries(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg);

// return true if successfully found a valid labeling
bool auto_fix_validity(Mesh& mesh, std::vector<vec3>& normals, const char* attribute_name, StaticLabelingGraph& slg, unsigned int max_nb_loop, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<vec3>& facet_normals);

unsigned int move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg);

void straighten_boundary_with_GCO(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg, index_t boundary_index);

// return true if successful
bool straighten_boundary(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, index_t boundary_index, const std::vector<std::vector<index_t>>& adj_facets);

void straighten_boundaries(GEO::Mesh& mesh, const char* attribute_name, StaticLabelingGraph& slg, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges);

// Typically called after straighten_boundary()
// Adjust postion of corners so that they are more aligned with adjacent boundaries
void move_corners(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets);

// only for turning-points on feature edges
// on a non-monotone boundary, find turning-points on feature edges, and for each of them, move the closest corner to the same vertex
// return bool if a turning point was merge with a corner, else false
bool merge_a_turning_point_and_a_corner_on_non_monotone_boundary(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets, index_t non_monotone_boundary_index);

// only for turning-points on feature edges
void merge_turning_points_and_corners(GEO::Mesh& mesh, const char* attribute_name, StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets);

// true if we reached 0 turning-point
bool auto_fix_monotonicity(Mesh& mesh, const char* attribute_name, StaticLabelingGraph& slg, std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges);

void increase_chart_valence(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets, index_t invalid_chart_index);

unsigned int count_lost_feature_edges(const CustomMeshHalfedges& mesh_he, const std::set<std::pair<index_t,index_t>>& feature_edges);