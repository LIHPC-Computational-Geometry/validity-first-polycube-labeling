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

#include <iostream>
#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <utility> // for std::pair
#include <limits>

#include "geometry.h" // for HalfedgeCompare
#include "CustomMeshHalfedges.h"

using namespace GEO;

namespace LabelingGraph {
    // For some associations, we need an "undefined" value
    // If we use -1, we need to convert index_t to signed_index_t, loosing half of the range
    // So instead, we use the maximal value for an index_t as "undefined"
    // Note that std::numeric_limits<index_t>::max() should be replaced by GEO::max_index_t(), but the latter loses constexpr !
    constexpr index_t UNDEFINED = std::numeric_limits<index_t>::max();
}

// to print an index_t that can be undefined
#define OPTIONAL_TO_STRING(value) ((value == LabelingGraph::UNDEFINED) ? "undefined" : std::to_string(value).c_str())

// A set of adjacent triangles having the same label
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L61
struct Chart {

    //// Properties //////////////////

    index_t label = LabelingGraph::UNDEFINED;

    //// Underlying geometry //////////////////

    std::set<index_t> facets;

    //// Adjacency //////////////////

    std::set<index_t> boundaries; // indices of adjacent boundaries
    // adjacent charts can be retrieved by "looking at the other side" of each boundary
};

std::ostream& operator<< (std::ostream &out, const Chart& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<Chart> : ostream_formatter {};

struct VertexRingWithBoundaries {

    //// Underlying geometry //////////////////

    std::vector<CustomMeshHalfedges::Halfedge> boundary_edges;

    //// Methods //////////////////

    std::size_t valence() const; // number of adjacent boundary edges

    bool halfedge_is_in_boundary_edges(const CustomMeshHalfedges::Halfedge& halfedge) const;

    void explore(const CustomMeshHalfedges::Halfedge& initial_halfedge,
                 const CustomMeshHalfedges& mesh_halfedges);

    void check_boundary_edges(const CustomMeshHalfedges& mesh_halfedges) const;
};

std::ostream& operator<< (std::ostream &out, const VertexRingWithBoundaries& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<VertexRingWithBoundaries> : ostream_formatter {};

// A vertex where 3 or more boundaries are meeting
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L126
struct Corner {

    //// Underlying geometry //////////////////

    index_t vertex;

    //// Adjacency //////////////////

    std::vector<VertexRingWithBoundaries> vertex_rings_with_boundaries; // a corner can have multiple groups of boundary edges if several solids are adjacent by a vertex. On the intersecting vertex, there is a vertex ring per adjacent solid

    //// Methods //////////////////

    std::size_t valence() const; // number of adjacent boundary edges

    bool halfedge_is_in_boundary_edges(const CustomMeshHalfedges::Halfedge& halfedge) const;
};

std::ostream& operator<< (std::ostream &out, const Corner& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<Corner> : ostream_formatter {};

// A boundary is a list of halfedges between 2 charts
struct Boundary {

    //// Properties //////////////////

    int axis = -1; // in {0 = X, 1 = Y, 2 = Z, -1 = boundary between opposite labels }

    //// Underlying geometry //////////////////

    std::vector<CustomMeshHalfedges::Halfedge> halfedges;

    //// Adjacency //////////////////

    index_t left_chart = LabelingGraph::UNDEFINED;
    index_t right_chart = LabelingGraph::UNDEFINED;
    index_t start_corner = LabelingGraph::UNDEFINED;
    index_t end_corner = LabelingGraph::UNDEFINED;

    //// Methods //////////////////

    bool empty() const;

    void explore(const CustomMeshHalfedges::Halfedge& initial_halfedge,
                 const CustomMeshHalfedges& mesh_halfedges,
                 index_t index_of_self,
                 const std::vector<index_t>& facet2chart,
                 std::vector<index_t>& vertex2corner,
                 std::vector<Chart>& charts,
                 std::vector<Corner>& corners,
                 std::map<CustomMeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary,
                 std::vector<CustomMeshHalfedges::Halfedge>& boundary_edges_to_explore);
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
    std::map<CustomMeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare> halfedge2boundary; // for each halfedge, store 1/ the associated boundary (may be LabelingGraph::UNDEFINED) and 2/ if the boundary is in the same orientation
    std::vector<index_t> vertex2corner;

    //// Quick access to problematic LabelingGraph features //////////////////

    vector<index_t> invalid_charts;
    // TODO vector<int> boundaries_between_opposite_labels
    // TODO vector<int> invalid_corners

    //// Fill & clear //////////////////

    void fill_from(Mesh& mesh, std::string facet_attribute);
    void clear();

    //// Getters for sizes //////////////////

    std::size_t nb_charts() const;
    std::size_t nb_boundaries() const;
    std::size_t nb_corners() const;
    std::size_t nb_facets() const; // to iterate over facet2chart
    std::size_t nb_vertices() const; // to iterate over vertex2corner
    std::size_t nb_invalid_charts() const;

    //// Export //////////////////

    void dump_to_file(const char* filename) const;
};

std::ostream& operator<< (std::ostream &out, const StaticLabelingGraph& data);

// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<StaticLabelingGraph> : ostream_formatter {};