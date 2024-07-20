#pragma once

#include <geogram/mesh/mesh.h>      // for GEO::Mesh

#include <vector>
#include <set>
#include <utility>  // for std::pair

#include "labeling_graph.h"  // for LabelingGraph
#include "geometry_mesh_ext.h" // for MeshExt

using namespace GEO;

size_t remove_surrounded_charts(Attribute<index_t>& labeling, const LabelingGraph& lg);

// return true if an invalid boundary was processed
// -> if returned false, no need to call the function again
bool fix_an_invalid_boundary(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);

// fix as much invalid corners as possible until the labeling graph needs to be recomputed
// return the number of invalid corners processed
// -> if returned 0, no need to call the function again
size_t fix_as_much_invalid_corners_as_possible(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);

// return true if all invalid charts are locked, ie we cannot remove them because they are surrounded by feature edges
bool remove_invalid_charts(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg);

void remove_charts_around_invalid_boundaries(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg);

// return true if a chart has been processed
// -> if returned false, no need to call the function again
bool increase_chart_valence(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg);

// return true if successfully found a valid labeling
bool auto_fix_validity(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg, unsigned int max_nb_loop);