#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <DisjointSet.hpp>

#include <algorithm>

#include "LabelingGraph.h"

StaticLabelingGraph::StaticLabelingGraph(Mesh& mesh, const char* facet_attribute) : mesh_half_edges_(mesh) {

    // based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L795
    // but here we use a disjoint-set for step 1

    facet2chart_.resize(mesh.facets.nb()); // important: memory allocation allowing to call ds.getSetsMap() on the underlying array
    edge2boundary_.resize(mesh.facets.nb()*3);
    vertex2corner_.resize(mesh.vertices.nb());
    std::vector<index_t> vertex_explored(mesh.vertices.nb(),false);

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
    std::size_t nb_charts = ds.getSetsMap(facet2chart_.data()); // get the map (facet id -> chart id) and the number of charts
    charts_.resize(nb_charts);

    // fill the Chart objects
    for(index_t f: mesh.facets) { // for each facet of the mesh
        Chart& current_chart = charts_[facet2chart_[f]]; // get the chart associated to this facet
        current_chart.label = labeling[f]; // (re)define its label
        current_chart.facets.emplace(f); // register the facet
    }

    // STEP 2 : Identify corners
    // Ideally, we have to iterate over all vertices, and,
    // for each of them, explore adjacent facets to compute valence (number of boundary edges).
    // But I don't know how to get the adjacent facets of a given vertex.
    // What I can do, with half-edges, is to iterate over all facet corners
    // check if the vertex at this corner was not already explored (because several facet corners are incident to the same vertex),
    // and go around this vertex with MeshHalfedges::move_to_next_around_vertex()

    mesh_half_edges_.set_use_facet_region(facet_attribute); // could be useful when exploring boundaries (= borders for Geogram)

    int previous_label = -1;
    int current_label = -1;
    int valence = -1;
    for(index_t f: mesh.facets) { // for each facet

        for(index_t c: mesh.facets.corners(f)) { // for each corner of the current facet

            MeshHalfedges::Halfedge initial_half_edge(f,c); // halfedge on facet f having corner c as base
            index_t current_vertex = mesh.facet_corners.vertex(initial_half_edge.corner); // get the vertex at the base of initial_half_edge

            if(vertex_explored[current_vertex]) {
                continue;
            }
            
            // Go around the vertex to compute its valence

            previous_label = labeling[initial_half_edge.facet];
            valence = 0; // number of boundary edges around this vertex
            MeshHalfedges::Halfedge current_half_edge = initial_half_edge; // current_half_edge will be our iterator around current_vertex
            do {
                mesh_half_edges_.move_to_next_around_vertex(current_half_edge); // moving to the next halfedge allows us to move to the next facet
                current_label = labeling[current_half_edge.facet];
                if (previous_label != current_label) { // check if current_half_edge is a boundary, ie between facets of different label/chart
                    valence++;
                }
                previous_label = current_label; // prepare next iteration
            } while (current_half_edge != initial_half_edge); // until we have not gone all the way around current_vertex
            
            if (valence >= 3) {
                corners_.push_back(Corner()); // add a corner
                corners_.back().vertex = current_vertex; // link this corner to current_vertex
                // TODO fill the boundary_edges attribute during the exploration
                vertex2corner_[current_vertex] = (index_t) corners_.size()-1; // link this vertex to the new (last) corner
                vertex_explored[current_vertex] = true; // mark this vertex
            }   
        }
    }

    // STEP 3 : Explore boundaries to connect corners
    // TODO

    // STEP 4 : Find boundaries with no corners
    // TODO
}

std::size_t StaticLabelingGraph::nb_charts() const {
    return charts_.size();
}

std::size_t StaticLabelingGraph::nb_corners() const {
    return corners_.size();
}

std::size_t StaticLabelingGraph::nb_facets() const {
    return facet2chart_.size();
}

const Chart& StaticLabelingGraph::chart(std::size_t index) const {
    return charts_.at(index);
}

const Corner& StaticLabelingGraph::corner(std::size_t index) const {
    return corners_.at(index);
}

index_t StaticLabelingGraph::facet2chart(std::size_t index) const {
    return facet2chart_.at(index);
}

const index_t* StaticLabelingGraph::get_facet2chart_ptr() const {
    return facet2chart_.data();
}