#pragma once

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram/basic/attributes.h>

#include <iostream>
#include <ostream>
#include <vector>
#include <set>
#include <limits>

using namespace GEO;

namespace LabelingGraph {
    // For some associations, we need an "undefined" value
    // If we use -1, we need to convert index_t to signed_index_t, loosing half of the range
    // So instead, we use the maximal value for an index_t as "undefined"
    // Note that std::numeric_limits<index_t>::max() should be replaced by GEo::max_index_t(), but the latter loses constexpr !
    constexpr index_t UNDEFINED = std::numeric_limits<index_t>::max();
}

// A set of adjacent triangles having the same label
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L61
struct Chart {

    //// Properties //////////////////

    int label = -1;

    //// Underlying geometry //////////////////

    std::set<index_t> facets;

    //// Adjacency //////////////////

    std::set<int> neighbors; // indices of adjacent charts
    // TODO store adjacent boundaries
};

// A vertex where 3 or more boundaries are meeting
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L126
struct Corner {

    //// Underlying geometry //////////////////

    index_t vertex;

    //// Adjacency //////////////////

    std::vector<std::vector<MeshHalfedges::Halfedge>> boundary_edges; // a corner can have multiple groups of boundary edges if several solids are adjacent by a vertex. On the intersecting vertex, there is a vertex ring per adjacent solid

    //// Methods //////////////////

    std::size_t valence() const ; // number of adjacent boundary edges
    bool halfedge_is_in_boundary_edges(const MeshHalfedges::Halfedge& halfedge) const ;
    void find_boundary_edges_in_vertex_ring(const MeshHalfedges::Halfedge& initial_halfedge, const MeshHalfedges& mesh_halfedges, const Attribute<int>& labeling); // explore a vertex ring with halfedges
};

// A labeling stored with charts, boundaries and corners
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L135
class StaticLabelingGraph {

public:
    StaticLabelingGraph(Mesh& mesh, const char* facet_attribute);

    //// Getters for sizes //////////////////

    std::size_t nb_charts() const;
    // TODO GEO::index_t nb_boundaries() const;
    std::size_t nb_corners() const;
    std::size_t nb_facets() const; // to iterate over facet2chart
    std::size_t nb_vertices() const; // to iterate over vertex2corner

    //// Getters for objects //////////////////

    const Chart& chart(std::size_t index) const;
    const Corner& corner(std::size_t index) const;
    index_t facet2chart(std::size_t index) const;
    index_t vertex2corner(std::size_t index) const;

    //// Special getters //////////////////

    const index_t* get_facet2chart_ptr() const; // useful to copy the whole chunk of data

private:

    //// LabelingGraph features //////////////////

    std::vector<Chart> charts_;
    // TODO std::vector<Boundary> boundaries;
    std::vector<Corner> corners_;

    //// Mapping from geometry to LabelingGraph features //////////////////

    std::vector<index_t> facet2chart_;
    std::vector<index_t> edge2boundary_;
    std::vector<index_t> vertex2corner_;

    //// Quick access to problematic LabelingGraph features //////////////////

    // TODO vector<int> invalid_charts
    // TODO vector<int> boundaries_between_opposite_labels
    // TODO vector<int> invalid_corners

    //// Half edges API //////////////////
    MeshHalfedges mesh_half_edges_;
};

// output stream overloading

inline std::ostream& operator<< (std::ostream &out, const Chart& data) {
    out << "\t label " << data.label;
    out << "\n\t facets : "; for (auto f: data.facets) out << f << " ";
    out << "\n\t neighbors : "; for (auto n: data.neighbors) out << n << " ";
    out << '\n';
    return out;
}

inline std::ostream& operator<< (std::ostream &out, const Corner& data) {
    out << "\t vertex " << data.vertex;
    out << "\n\t boundary_edges";
    for (std::size_t i = 0; i < data.boundary_edges.size(); ++i) { // iterate over boundary edges, grouped by vertex ring
        out << "\n\t\t[" << i << "] "; for (auto be: data.boundary_edges[i]) out << be << " ";
    }
    out << '\n';
    return out;
}

inline std::ostream& operator<< (std::ostream &out, const StaticLabelingGraph& data) {
    for(std::size_t chart_index = 0; chart_index < data.nb_charts(); ++chart_index) {
        out << "charts[" << chart_index << "]\n";
        out << data.chart(chart_index);
    }
    for(std::size_t corner_index = 0; corner_index < data.nb_corners(); ++corner_index) {
        out << "corners[" << corner_index << "]\n";
        out << data.corner(corner_index);
    }
    out << "facet2chart\n";
    for(std::size_t f = 0; f < data.nb_facets(); ++f) {
        out << '[' << f << "] " << data.facet2chart(f) << '\n';
    }
    out << "vertex2corner\n";
    index_t corner;
    for(std::size_t v = 0; v < data.nb_vertices(); ++v) {
        corner = data.vertex2corner(v);
        out << '[' << v << "] ";
        if (corner == LabelingGraph::UNDEFINED) out << "undefined";
        else out << corner;
        out << '\n';
    }
    return out;
}