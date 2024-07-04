#pragma once

#include <geogram/mesh/mesh.h>      // for GEO::Mesh
#include "labeling_graph.h"

using namespace GEO;

// return nb turning-points moved
size_t move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges);

void straighten_boundary_with_GCO(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const LabelingGraph& lg, index_t boundary_index);

// return true if successful
bool straighten_boundary(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, index_t boundary_index, const std::vector<std::vector<index_t>>& adj_facets);

void straighten_boundaries(GEO::Mesh& mesh, const char* attribute_name, LabelingGraph& lg, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges);

// Typically called after straighten_boundary()
// Adjust postion of corners so that they are more aligned with adjacent boundaries
void move_corners(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets);

// only for turning-points on feature edges
// on a non-monotone boundary, find turning-points on feature edges, and for each of them, move the closest corner to the same vertex
// return true if a turning point was merge with a corner
// -> if returned false, no need to call the function again
bool merge_a_turning_point_and_its_closest_corner(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets);

// return true if a pair of turning-point was processed
// -> if returned false, no need to call the function again
bool join_turning_points_pair_with_new_chart(GEO::Mesh& mesh, const char* attribute_name, LabelingGraph& lg, const std::vector<vec3>& normals, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets);

// return true if we reached 0 turning-point
bool auto_fix_monotonicity(Mesh& mesh, const char* attribute_name, LabelingGraph& lg, std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<vec3>& normals);