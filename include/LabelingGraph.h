#pragma once

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram/basic/attributes.h>

#include <iostream>
#include <ostream>
#include <vector>
#include <set>

using namespace GEO;

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

    std::vector<MeshHalfedges::Halfedge> boundary_edges;

    //// Methods //////////////////

    std::size_t valence(); // number of adjacent boundary edges
    void autofill_from_vertex(const MeshHalfedges::Halfedge&, const MeshHalfedges& mesh_he, const Attribute<int>& labeling); // explore a vertex with halfedges. /!\ the result is not a corner if valence() < 3
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

    //// Getters for objects //////////////////

    const Chart& chart(std::size_t index) const;
    const Corner& corner(std::size_t index) const;
    index_t facet2chart(std::size_t index) const;

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
    out << "\n\t boundary_edges : "; for (auto be: data.boundary_edges) out << be << " ";
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
    return out;
}