#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include <initializer_list>

#include "LabelingGraph.h"  // for Boundary
#include "stats.h"          // for BasicStats

// get the string representation of a label
#define LABEL2STR(label) ((label==0) ? "+X" : ( \
                          (label==1) ? "-X" : ( \
                          (label==2) ? "+Y" : ( \
                          (label==3) ? "-Y" : ( \
                          (label==4) ? "+Z" : ( \
                          (label==5) ? "-Z" :   \
                          "undef."))))))

using namespace GEO;

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

index_t opposite_label(index_t label);

void flip_labeling(Mesh& mesh, const char* attribute_name);

index_t find_optimal_label(std::initializer_list<index_t> forbidden_axes = {}, std::initializer_list<index_t> forbidden_labels = {}, std::initializer_list<index_t> orthogonal_labels = {}, vec3 close_vector = vec3(0.0,0.0,0.0));

void propagate_label(const Mesh& mesh, const char* attribute_name, index_t new_label, const std::set<index_t>& facets_in, const std::set<index_t> facets_out, const std::vector<index_t>& facet2chart, index_t chart_index);

void compute_per_facet_fidelity(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* labeling_attribute_name, const char* fidelity_attribute_name, BasicStats& stats);

unsigned int count_lost_feature_edges(const MeshHalfedgesExt& mesh_he, const std::set<std::pair<index_t,index_t>>& feature_edges);

// TODO move to LabelingGraph.h
index_t adjacent_chart_in_common(const Boundary& b0, const Boundary& b1);