#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/assert.h>

#include <DisjointSet.hpp>

#include <algorithm> // for std::find()
#define CONTAINS(container,value) (std::find(container.cbegin(),container.cend(),value) != container.cend())

#include "LabelingGraph.h"

std::size_t Corner::valence() const {
    std::size_t size = 0;
    for (const std::vector<MeshHalfedges::Halfedge>& vertex_ring_boundary_edges: boundary_edges) {
        size += vertex_ring_boundary_edges.size();
    }
    return size;
}

bool Corner::halfedge_is_in_boundary_edges(const MeshHalfedges::Halfedge& halfedge) const {
    for(const std::vector<MeshHalfedges::Halfedge>& vertex_ring_boundary_edges: boundary_edges) { // parse boundary edges grouped by vertex ring
        if(CONTAINS(vertex_ring_boundary_edges,halfedge)) { // if the boundary edges of the current vertex ring contain the input halfedge
            return true;
        }
    }
    return false;
}

/**
 * \brief Explore the vertex at the base of \c initial_halfedge and fill the struct with a potential corner. /!\ The vertex field should be consistent with the halfedge
 * \param[in] initial_halfedge An half-edge going outwards the vertex to explore
 * \param[in] mesh_halfedges The half-edges interface of the mesh
 * \param[in] labeling The facet-to-label association, to detect boundary edges (= edges between 2 labels)
 */
void Corner::find_boundary_edges_in_vertex_ring(const MeshHalfedges::Halfedge& initial_halfedge, const MeshHalfedges& mesh_halfedges, const Attribute<int>& labeling) {

    // test if Corner::vertex is consistent with the initial_halfedge
    geo_assert(vertex == mesh_halfedges.mesh().facet_corners.vertex(initial_halfedge.corner));
    
    // prepare the vertex exploration
    std::vector<MeshHalfedges::Halfedge> boundary_edges_in_current_vertex_ring;
    int previous_label = labeling[initial_halfedge.facet];
    int current_label = -1;
    MeshHalfedges::Halfedge current_halfedge = initial_halfedge; // create another Halfedge that we can modify (we need to keep the initial one)
    
    // go around the vertex
    do {
        mesh_halfedges.move_to_next_around_vertex(current_halfedge); // moving to the next halfedge allows us to move to the next facet
        current_label = labeling[current_halfedge.facet]; // get the label of the facet associated to this current_halfedge
        if (previous_label != current_label) { // check if current_halfedge is a boundary, ie between facets of different label
            // We found a boundary edge
            // But because we are re-exploring vertices to take into account all vertex rings of each vertex
            // maybe we are re-exploring the same vertex ring.
            if(halfedge_is_in_boundary_edges(current_halfedge)) {
                return; // current vertex ring already explored
            }
            boundary_edges_in_current_vertex_ring.push_back(current_halfedge); // register the current halfedge as outgoing boundary edge
        }
        previous_label = current_label; // prepare next iteration
    } while (current_halfedge != initial_halfedge);

    if(!boundary_edges_in_current_vertex_ring.empty()) { // if we found boundary edges in this vertex ring
        boundary_edges.push_back(boundary_edges_in_current_vertex_ring); // store them in the Corner struct
    }
}

StaticLabelingGraph::StaticLabelingGraph(Mesh& mesh, const char* facet_attribute) : mesh_half_edges_(mesh) {

    // based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L795
    // but here we use a disjoint-set for step 1

    facet2chart_.resize(mesh.facets.nb()); // important: memory allocation allowing to call ds.getSetsMap() on the underlying array
    edge2boundary_.resize(mesh.facets.nb()*3,LabelingGraph::UNDEFINED); // contrary to facet2chart_ where all facets are associated to a chart, not all edges are associated to a boundary -> use UNDEFINED for edges that are not on a boundary
    vertex2corner_.resize(mesh.vertices.nb(),LabelingGraph::UNDEFINED); // contrary to facet2chart_ where all facets are associated to a chart, not all vertices are associated to a corner -> use UNDEFINED for vertices that are not on a corner

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
    // What I can do, with half-edges, is to iterate over all facet corners, get the vertex at this corner,
    // and go around the vertex ring with MeshHalfedges::move_to_next_around_vertex(). See Corner::find_boundary_edges_in_vertex_ring()
    // This works for shapes having several solids connected by a vertex only

    mesh_half_edges_.set_use_facet_region(facet_attribute); // could be useful when exploring boundaries (= borders for Geogram)

    vector<Corner> potential_corners(mesh.vertices.nb());
    for(index_t f: mesh.facets) { for(index_t c: mesh.facets.corners(f)) { // for each facet corner (f,c)
        
        MeshHalfedges::Halfedge initial_half_edge(f,c); // halfedge on facet f having corner c as base
        index_t current_vertex = mesh.facet_corners.vertex(initial_half_edge.corner); // get the vertex at the base of initial_half_edge
        potential_corners[current_vertex].vertex = current_vertex; // (re)write the corresponding vertex index
        
        // (re)explore this vertex, because maybe it was explored from another solid of the same shape (multiple vertex rings)
        potential_corners[current_vertex].find_boundary_edges_in_vertex_ring(initial_half_edge,mesh_half_edges_,labeling);
        
    }}

    // from potential corners to corners : keep only valence >= 3
    for(index_t v: mesh.vertices) {
        if(potential_corners[v].valence() >= 3) {
            corners_.push_back(potential_corners[v]); // register this corner
            geo_assert(v == potential_corners[v].vertex);
            vertex2corner_[v] = (index_t) corners_.size()-1; // link the underlying vertex to this corner
        }
    }

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

std::size_t StaticLabelingGraph::nb_vertices() const {
    return vertex2corner_.size();
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

index_t StaticLabelingGraph::vertex2corner(std::size_t index) const {
    return vertex2corner_.at(index);
}

const index_t* StaticLabelingGraph::get_facet2chart_ptr() const {
    return facet2chart_.data();
}