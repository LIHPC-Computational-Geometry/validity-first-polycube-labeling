/* 
 * Where a labeling is a simple data structure that maps a surface triangle to a label in [0,5] = +/- {X,Y,Z},
 * a LabelingGraph is a data structure that stores 
 *  - charts (a chart is set of adjacent triangles having the same label)
 *  - boundaries (a boundary is a set of edges between 2 charts)
 *  - corners (a corner is a vertex where several boundaries are meeting)
 * 
 * Warning : in the context of meshes in Geogram, a corner is a vertex seen from a facet or a cell https://github.com/BrunoLevy/geogram/wiki/Mesh
 * 
 * StaticLabelingGraph is a LabelingGraph computed from a labeling, that is rather read-only. If you update the labeling, you need to re-compute the whole LabelingGraph.
 * Maybe one day there will be a DynamicLabelingGraph with operators on charts, locally updating charts/boundaries/corners.
 */

#pragma once

#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/os.h>

#include <iostream>
#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <utility> // for std::pair
#include <limits>
#include <optional>

#include "geometry.h" // for HalfedgeCompare
#include "CustomMeshHalfedges.h"

#define FALLOFF_BINARY 1.0 // smaller value -> more turning points. https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/graphcut_labeling.cpp#L183

using namespace GEO;

// to print an index_t that can be undefined
#define OPTIONAL_TO_STRING(value) ((value == index_t(-1)) ? "undefined" : std::to_string(value).c_str())

// we must declare that 'Boundary' will be defined later
struct Boundary;

// A set of adjacent triangles having the same label
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L61
struct Chart {

    //// Properties //////////////////

    index_t label = index_t(-1);
    // validity is not stored because it's inexpensive to compute (boundaries.size() >= 4), and not needed for the visualization (facet attribute instead)

    //// Underlying geometry //////////////////

    std::set<index_t> facets;

    //// Adjacency //////////////////

    std::set<index_t> boundaries; // indices of adjacent boundaries
    // adjacent charts can be retrieved with Boundary::chart_at_the_other_side()

    //// Methods //////////////////

    bool is_surrounded_by_feature_edges(const std::vector<Boundary>& all_boundaries) const;
    void counterclockwise_boundaries_order(const CustomMeshHalfedges& mesh_he, const std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare> halfedge2boundary, const std::vector<Boundary>& all_boundaries, std::vector<std::pair<index_t,bool>>& counterclockwise_order) const; // `counterclockwise_order` list boundary indices & whether or not the boundary is in the same direction when going counterclockwise
};

std::ostream& operator<< (std::ostream &out, const Chart& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<Chart> : ostream_formatter {};

struct VertexRingWithBoundaries {

    //// Underlying geometry //////////////////

    std::vector<MeshHalfedges::Halfedge> boundary_edges; // outgoing boundary edges

    //// Methods //////////////////

    std::size_t valence() const; // number of adjacent boundary edges

    bool halfedge_is_in_boundary_edges(const MeshHalfedges::Halfedge& halfedge) const;

    void explore(const MeshHalfedges::Halfedge& initial_halfedge,
                 const CustomMeshHalfedges& mesh_halfedges);

    // Is it possible to get a facet corner from a vertex ? If so, use a vertex index as first argument
    void explore(index_t init_facet_corner, const CustomMeshHalfedges& mesh_halfedges);

    void check_boundary_edges(const CustomMeshHalfedges& mesh_halfedges) const;

    std::size_t circular_previous(std::size_t index) const;

    std::size_t circular_next(std::size_t index) const;
};

std::ostream& operator<< (std::ostream &out, const VertexRingWithBoundaries& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<VertexRingWithBoundaries> : ostream_formatter {};

struct StaticLabelingGraph; // defined later

// A vertex where 3 or more boundaries are meeting
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L126
struct Corner {

    //// Properties //////////////////

    bool is_valid = false;

    //// Underlying geometry //////////////////

    index_t vertex = index_t(-1);

    //// Adjacency //////////////////

    std::vector<VertexRingWithBoundaries> vertex_rings_with_boundaries; // a corner can have multiple groups of boundary edges if several solids are adjacent by a vertex. On the intersecting vertex, there is a vertex ring per adjacent solid

    //// Methods //////////////////

    std::size_t valence() const; // number of adjacent boundary edges

    bool halfedge_is_in_boundary_edges(const MeshHalfedges::Halfedge& halfedge) const;

    bool compute_validity(bool allow_boundaries_between_opposite_labels, const std::vector<Boundary>& boundaries, const std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary);

    bool all_adjacent_boundary_edges_are_on_feature_edges(const Mesh& mesh, const std::set<std::pair<index_t,index_t>>& feature_edges) const;

    // retreive near vertices along adjacent boundaries and compute the average coordinate
    vec3 average_coordinates_of_neighborhood(const Mesh& mesh, const StaticLabelingGraph& slg, bool include_itself, size_t max_dist) const;

    MeshHalfedges::Halfedge get_most_aligned_boundary_halfedge(const Mesh& mesh, const vec3& reference) const;

    void get_outgoing_halfedges_on_feature_edge(const Mesh& mesh, const std::set<std::pair<index_t,index_t>>& all_feature_edges, std::vector<MeshHalfedges::Halfedge>& outgoing_halfedges_on_feature_edge) const;
};

std::ostream& operator<< (std::ostream &out, const Corner& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<Corner> : ostream_formatter {};

struct TurningPoint {

    index_t outgoing_local_halfedge_index_; // index in Boundary::halfedges
    bool is_towards_right_;
    
    void fill_from(index_t outgoing_local_halfedge_index, const Boundary& boundary, const CustomMeshHalfedges& mesh_he); // important : mesh_he must use facet region (on the labeling)
    bool is_towards_left() const  { return is_towards_right_ == false; }
    bool is_towards_right() const { return is_towards_right_; }

    // `boundary` must be the boundary on which the turning point is
    index_t get_closest_corner(const Boundary& boundary, const CustomMeshHalfedges& mesh_he) const;
    index_t vertex(const Boundary& boundary, const Mesh& mesh) const; // get the vertex on which the turning point is
};

std::ostream& operator<< (std::ostream &out, const TurningPoint& data);

bool operator==(const TurningPoint& a, const TurningPoint& b);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<TurningPoint> : ostream_formatter {};

enum BoundarySide {
    LeftAndRight,
    OnlyLeft,
    OnlyRight
};

// A boundary is a list of halfedges between 2 charts
struct Boundary {

    //// Properties //////////////////

    int axis = -1; // in {0 = X, 1 = Y, 2 = Z, -1 = boundary between opposite labels }
    bool is_valid = false;
    std::vector<TurningPoint> turning_points;
    bool on_feature_edge = true;

    //// Underlying geometry //////////////////

    std::vector<MeshHalfedges::Halfedge> halfedges;
    vec3 average_normal = {0,0,0};

    //// Adjacency //////////////////

    index_t left_chart = index_t(-1);
    index_t right_chart = index_t(-1);
    index_t start_corner = index_t(-1);
    index_t end_corner = index_t(-1);

    //// Methods //////////////////

    bool empty() const;

    void explore(const MeshHalfedges::Halfedge& initial_halfedge,
                 const CustomMeshHalfedges& mesh_halfedges,
                 index_t index_of_self,
                 const std::set<std::pair<index_t,index_t>>& feature_edges,
                 const std::vector<index_t>& facet2chart,
                 std::vector<index_t>& vertex2corner,
                 std::vector<Chart>& charts,
                 std::vector<Corner>& corners,
                 std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary,
                 std::vector<MeshHalfedges::Halfedge>& boundary_edges_to_explore);

    bool contains_lower_than_180_degrees_angles(const CustomMeshHalfedges& mesh_halfedges);

    bool compute_validity(bool allow_boundaries_between_opposite_labels, const CustomMeshHalfedges& mesh_halfedges);

    bool find_turning_points(const CustomMeshHalfedges& mesh_halfedges); // return true if the boundary contains turning-points

    index_t chart_at_other_side(index_t origin_chart) const;

    // if `label` is the label of the `left_chart`, returns the label of the `right_chart`
    // if `label` is the label of the `right_chart`, returns the label of the `left_chart`
    index_t other_label(const std::vector<Chart>& charts, index_t label);

    void print_successive_halfedges(fmt::v9::ostream& out, Mesh& mesh);

    index_t turning_point_vertex(index_t turning_point_index, const Mesh& mesh) const;

    bool halfedge_has_turning_point_at_base(index_t local_halfedge_index) const;

    index_t get_closest_boundary_of_turning_point(const TurningPoint& turning_point, index_t closest_corner, const CustomMeshHalfedges& mesh_he, const std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary, const std::vector<Corner>& corners) const;

    vec3 vector_between_corners(const Mesh& mesh, const std::vector<Corner>& corners) const;

    vec3 average_vector_between_corners(const Mesh& mesh, const std::vector<Corner>& corners) const;

    void get_adjacent_facets(const Mesh& mesh, std::set<index_t>& adjacent_facets, BoundarySide boundary_side_to_explore, const std::vector<index_t>& facet2chart, size_t max_dist = 0) const;

    void get_flipped(const MeshHalfedges& mesh_he, Boundary& flipped_boundary) const;

    void split_at_turning_point(const MeshHalfedges& mesh_he, Boundary& downward_boundary, Boundary& upward_boundary, index_t local_turning_point_index = 0) const;

    void per_edges_axis_assignement_cost(const Mesh& mesh, index_t axis, std::vector<double>& costs) const;

    void per_edges_cumulative_axis_assignement_cost(const Mesh& mesh, index_t axis, std::vector<double>& costs, bool accumulation_from_start_corner) const;

    double length(const Mesh& mesh) const;
};

std::ostream& operator<< (std::ostream &out, const Boundary& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<Boundary> : ostream_formatter {};

// A labeling stored with charts, boundaries and corners
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L135
struct StaticLabelingGraph {

    //// LabelingGraph features //////////////////

    std::vector<Chart> charts;
    std::vector<Boundary> boundaries;
    std::vector<Corner> corners;

    //// Mapping from geometry to LabelingGraph features //////////////////

    std::vector<index_t> facet2chart;
    std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare> halfedge2boundary; // for each halfedge, store 1/ the associated boundary (may be index_t(-1)) and 2/ if the boundary is in the same orientation
    std::vector<index_t> vertex2corner;

    //// problematic LabelingGraph features //////////////////

    vector<index_t> invalid_charts;
    vector<index_t> invalid_boundaries;
    vector<index_t> invalid_corners;
    vector<index_t> non_monotone_boundaries;
    std::map<index_t,std::vector<std::pair<index_t,index_t>>> turning_point_vertices; // map a vertex index on which there is a turning-point to an invalid boundary index and a local turning-point index

    //// Parameter given to last fill_from() //////////////////
    bool allow_boundaries_between_opposite_labels_ = false;

    //// Fill & clear //////////////////

    void fill_from(Mesh& mesh, std::string facet_attribute, const std::set<std::pair<index_t,index_t>>& feature_edges, std::optional<bool> allow_boundaries_between_opposite_labels = std::nullopt);
    void clear();

    //// Validity //////////////////

    bool is_valid();

    //// Getters for sizes //////////////////

    std::size_t nb_charts() const;
    std::size_t nb_boundaries() const;
    std::size_t nb_corners() const;
    std::size_t nb_facets() const; // to iterate over facet2chart
    std::size_t nb_vertices() const; // to iterate over vertex2corner
    std::size_t nb_invalid_charts() const;
    std::size_t nb_invalid_boundaries() const;
    std::size_t nb_invalid_corners() const;
    std::size_t nb_non_monotone_boundaries() const;
    std::size_t nb_turning_points() const;

    //// Getters for parameter(s) //////////////////

    bool is_allowing_boundaries_between_opposite_labels() const;

    //// Other getters //////////////////

    bool vertex_is_only_surrounded_by(index_t vertex_index, std::vector<index_t> expected_charts, const std::vector<std::vector<index_t>>& vertex_to_adj_facets) const;
    void get_adjacent_charts_of_vertex(index_t vertex_index, const std::vector<std::vector<index_t>>& vertex_to_adj_facets, std::set<index_t>& adjacent_charts) const;

    //// Export //////////////////

    void dump_to_text_file(const char* filename, Mesh& mesh);
    void dump_to_D3_graph(const char* filename); // not tested
};

std::ostream& operator<< (std::ostream &out, const StaticLabelingGraph& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<StaticLabelingGraph> : ostream_formatter {};