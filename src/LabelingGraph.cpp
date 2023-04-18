#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <DisjointSet.hpp>

#include <algorithm>

#include "LabelingGraph.h"

StaticLabelingGraph::StaticLabelingGraph(const Mesh& mesh, const char* facet_attribute) {

    // based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L795
    // but here we use a disjoint-set for step 1

    // STEP 1 : Aggregete adjacent triangles of same label as chart

    Attribute<int> labeling(mesh.facets.attributes(),facet_attribute); // read and store the facets attribute corresponding to the labeling

    DisjointSet<index_t> ds(labeling.size()); // create a disjoint-set data structure to group facets by labels
    for(index_t f: mesh.facets) { // for each facet of the mesh
        for(index_t le = 0; le < 3; ++le) { // for each local edge of facet f in {0,1,2}
            index_t adjacent_facet = mesh.facets.adjacent(f,le); // get the facet id of the facet beyond le
            if(labeling[f] == labeling[adjacent_facet]) { // if facets f and adjacent_facet have the same label
                ds.mergeSets(f,adjacent_facet); // merge facets f and adjacent_facet
            }
        }
    }
    facet2chart_.resize(labeling.size()); // important: memory allocation allowing to call ds.getSetsMap() on the underlying array
    std::size_t nb_charts = ds.getSetsMap(facet2chart_.data()); // get the map (facet id -> chart id) and the number of charts
    charts_.resize(nb_charts);

    // fill the Chart objects
    for(index_t f: mesh.facets) { // for each facet of the mesh
        Chart& current_chart = charts_[facet2chart_[f]]; // get the chart associated to this facet
        current_chart.label = labeling[f]; // (re)define its label
        current_chart.facets.emplace(f); // register the facet
    }

    // STEP 2 : Walk along boundary edges to assemble boundaries + Identify corners and turning-points
    // TODO

    // STEP 3 : Find boundaries with no corners
    // TODO
}

std::size_t StaticLabelingGraph::nb_charts() const {
    return charts_.size();
}

std::size_t StaticLabelingGraph::nb_facets() const {
    return facet2chart_.size();
}

const Chart& StaticLabelingGraph::chart(std::size_t index) const {
    return charts_.at(index);
}

index_t StaticLabelingGraph::facet2chart(std::size_t index) const {
    return facet2chart_.at(index);
}

const index_t* StaticLabelingGraph::get_facet2chart_ptr() const {
    return facet2chart_.data();
}