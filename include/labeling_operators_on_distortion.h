#pragma once

#include "labeling_graph.h"

#include "geometry_mesh_ext.h" // for MeshExt

using namespace GEO;

// return nb turning-points moved
size_t move_boundaries_near_turning_points(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg);

void straighten_boundary_with_GCO(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg, index_t boundary_index);

// return true if successful
bool straighten_boundary(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg, index_t boundary_index);

void straighten_boundaries(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);

// Typically called after straighten_boundary()
// Adjust postion of corners so that they are more aligned with adjacent boundaries
void move_corners(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg);

// only for turning-points on feature edges
// on a non-monotone boundary, find turning-points on feature edges, and for each of them, move the closest corner to the same vertex
// return true if a turning point was merge with a corner
// -> if returned false, no need to call the function again
bool merge_a_turning_point_and_its_closest_corner(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);

// return true if a pair of turning-point was processed
// -> if returned false, no need to call the function again
bool join_turning_points_pair_with_new_chart(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);

// return true if we reached 0 turning-point
bool auto_fix_monotonicity(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);