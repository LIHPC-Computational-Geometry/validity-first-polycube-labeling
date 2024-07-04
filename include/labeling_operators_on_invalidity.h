#pragma once

#include <geogram/mesh/mesh.h>      // for GEO::Mesh

#include <vector>
#include <set>
#include <utility>  // for std::pair

#include "labeling_graph.h"  // for LabelingGraph

using namespace GEO;

size_t remove_surrounded_charts(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg);

// return true if an invalid boundary was processed
// -> if returned false, no need to call the function again
bool fix_an_invalid_boundary(GEO::Mesh& mesh, const char* attribute_name, LabelingGraph& lg, const std::vector<vec3>& facet_normals, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets);

// fix as much invalid corners as possible until the labeling graph needs to be recomputed
// return the number of invalid corners processed
// -> if returned 0, no need to call the function again
size_t fix_as_much_invalid_corners_as_possible(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<index_t>& facet2chart, const std::vector<std::vector<index_t>>& adj_facets);

// return true if all invalid charts are locked, ie we cannot remove them because they are surrounded by feature edges
bool remove_invalid_charts(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const LabelingGraph& lg);

void remove_charts_around_invalid_boundaries(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const LabelingGraph& lg);

// return true if a chart has been processed
// -> if returned false, no need to call the function again
bool increase_chart_valence(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, LabelingGraph& lg, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges);

// return true if successfully found a valid labeling
bool auto_fix_validity(Mesh& mesh, std::vector<vec3>& normals, const char* attribute_name, LabelingGraph& lg, unsigned int max_nb_loop, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<vec3>& facet_normals, const std::vector<std::vector<index_t>>& adj_facets);