#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <DisjointSet.hpp>

#include <algorithm>

#include "LabelingGraph.h"

std::size_t Corner::valence() {
    return boundary_edges.size();
}

/**
 * \brief Explore the vertex at the base of \c initial_halfedge and fill the struct with a potential corner
 * \param[in] initial_halfedge An half-edge going outwards the vertex to explore
 * \param[in] mesh_halfedges The half-edges interface of the mesh
 * \param[in] labeling The facet-to-label association, to detect boundary edges (= edges between 2 labels)
 */
void Corner::autofill_from_vertex(const MeshHalfedges::Halfedge& initial_halfedge, const MeshHalfedges& mesh_halfedges, const Attribute<int>& labeling) {

    // init the attributes
    vertex = mesh_halfedges.mesh().facet_corners.vertex(initial_halfedge.corner); // get the index of the vertex at the base of initial_halfedge
    boundary_edges.clear();

    // prepare the vertex exploration
    int previous_label = labeling[initial_halfedge.facet];
    int current_label = -1;
    MeshHalfedges::Halfedge current_halfedge = initial_halfedge; // create another Halfedge that we can modify (we need to keep the initial one)
    
    // go around the vertex
    do {
        mesh_halfedges.move_to_next_around_vertex(current_halfedge); // moving to the next halfedge allows us to move to the next facet
        current_label = labeling[current_halfedge.facet]; // get the label of the facet associated to this current_halfedge
        if (previous_label != current_label) { // check if current_halfedge is a boundary, ie between facets of different label
            boundary_edges.push_back(current_halfedge); // register the current halfedge as outgoing boundary edge
        }
        previous_label = current_label; // prepare next iteration
    } while (current_halfedge != initial_halfedge);
}

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
    // and go around this vertex with MeshHalfedges::move_to_next_around_vertex(). See Corner::autofill_from_vertex()
    //
    // This does not work for shapes having parts connected by a vertex only
    // TODO (maybe optionnal) re-exploration of vertices

    mesh_half_edges_.set_use_facet_region(facet_attribute); // could be useful when exploring boundaries (= borders for Geogram)

    for(index_t f: mesh.facets) { for(index_t c: mesh.facets.corners(f)) { // for each facet corner (f,c)
        
        MeshHalfedges::Halfedge initial_half_edge(f,c); // halfedge on facet f having corner c as base
        index_t current_vertex = mesh.facet_corners.vertex(initial_half_edge.corner); // get the vertex at the base of initial_half_edge

        if(vertex_explored[current_vertex]) {
            continue;
        }
        
        // Go around the vertex to compute its valence
        Corner potential_corner;
        potential_corner.autofill_from_vertex(initial_half_edge,mesh_half_edges_,labeling);
        
        if (potential_corner.valence() >= 3) {
            corners_.push_back(potential_corner); // add a corner
            vertex2corner_[current_vertex] = (index_t) corners_.size()-1; // link this vertex to the new (last) corner
            vertex_explored[current_vertex] = true; // mark this vertex
        }
    }}

    // STEP 3 : Explore boundaries to connect corners, fill Chart::neighbors
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