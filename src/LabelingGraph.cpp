#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/assert.h>

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/os.h>

#include <DisjointSet.hpp>

#include <utility> // for std::pair

#include "LabelingGraph.h"
#include "geometry.h"               // for other_axis(), HalfedgeCompare
#include "containers.h"             // for VECTOR_CONTAINS(), MAP_CONTAINS(), index_of_last(), += on std::vector
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges and CustomMeshHalfedges::Halfedge

std::ostream& operator<< (std::ostream &out, const Chart& data) {
    fmt::println(out,"\tlabel : {}",data.label);
    fmt::println(out,"\tfacets : {}",data.facets);
    fmt::println(out,"\tboundaries : {}",data.boundaries);
    return out;
}

std::size_t VertexRingWithBoundaries::valence() const {
    return boundary_edges.size();
}

bool VertexRingWithBoundaries::halfedge_is_in_boundary_edges(const CustomMeshHalfedges::Halfedge& halfedge) const {
    if(VECTOR_CONTAINS(boundary_edges,halfedge)) {
        return true;
    }
    return false;
}

/**
 * \brief Explore the vertex at the base of \c initial_halfedge and fill the struct with boundary edges encountered
 * \param[in] initial_halfedge An half-edge going outwards the vertex to explore
 * \param[in] mesh_halfedges The half-edges interface of the mesh
 * \param[in] labeling The facet-to-label association, to detect boundary edges (= edges between 2 labels)
 */
void VertexRingWithBoundaries::explore(const CustomMeshHalfedges::Halfedge& initial_halfedge,
                                       const CustomMeshHalfedges& mesh_halfedges) {
    
    // prepare the vertex exploration
    CustomMeshHalfedges::Halfedge current_halfedge = initial_halfedge; // create another Halfedge that we can modify (we need to keep the initial one)

    // go around the vertex, from boundary edge to boundary edge
    do {
        mesh_halfedges.move_to_next_around_vertex(current_halfedge,true);
        if(mesh_halfedges.halfedge_is_border(current_halfedge)) {
            boundary_edges.push_back(current_halfedge); // register the current halfedge as outgoing boundary edge for this vertex ring
        }
    } while ((current_halfedge != initial_halfedge));
}

void VertexRingWithBoundaries::check_boundary_edges(const CustomMeshHalfedges& mesh_halfedges) const {
    for(const auto& be: boundary_edges) {
        geo_assert(mesh_halfedges.halfedge_is_border(be));
    }
}

std::ostream& operator<< (std::ostream &out, const VertexRingWithBoundaries& data) {
    fmt::print(out,"{}",data.boundary_edges);
    return out;
}

std::size_t Corner::valence() const {
    // sum the valence of all vertex rings
    std::size_t valence = 0;
    for (const VertexRingWithBoundaries& vertex_ring_with_boundaries: vertex_rings_with_boundaries) {
        valence += vertex_ring_with_boundaries.valence();
    }
    return valence;
}

bool Corner::halfedge_is_in_boundary_edges(const CustomMeshHalfedges::Halfedge& halfedge) const {
    for(const VertexRingWithBoundaries& vertex_ring_with_boundaries: vertex_rings_with_boundaries) { // parse boundary edges grouped by vertex ring
        if(vertex_ring_with_boundaries.halfedge_is_in_boundary_edges(halfedge)) {
            return true;
        }
    }
    return false;
}

std::ostream& operator<< (std::ostream &out, const Corner& data) {
    fmt::println(out,"\tvertex : {}",data.vertex);
    fmt::println(out,"\t{} vertex ring(s) with boundaries : {}",data.vertex_rings_with_boundaries.size(),data.vertex_rings_with_boundaries);
    return out;
}

bool Boundary::empty() const {
    return (
        axis == -1 &&
        halfedges.empty() &&
        left_chart == LabelingGraph::UNDEFINED &&
        right_chart == LabelingGraph::UNDEFINED &&
        start_corner == LabelingGraph::UNDEFINED &&
        end_corner == LabelingGraph::UNDEFINED
    );
}

void Boundary::explore(const CustomMeshHalfedges::Halfedge& initial_halfedge,
                       const CustomMeshHalfedges& mesh_halfedges,
                       index_t index_of_self,
                       const std::vector<index_t>& facet2chart,
                       std::vector<index_t>& vertex2corner,
                       std::vector<Chart>& charts,
                       std::vector<Corner>& corners,
                       std::map<CustomMeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary,
                       std::vector<CustomMeshHalfedges::Halfedge>& boundary_edges_to_explore) {

    geo_assert(empty());

    // initialization
    index_t current_vertex = LabelingGraph::UNDEFINED;
    const Mesh& mesh = mesh_halfedges.mesh();
    CustomMeshHalfedges::Halfedge current_halfedge = initial_halfedge; // create a modifiable halfedge

    // Step 1 : store information we already have with the first boundary edge:
    // the start corner, the charts at the left and right, and the axis

    // fill start_corner field
    start_corner = vertex2corner[mesh.facet_corners.vertex(initial_halfedge.corner)];
    geo_assert(start_corner != LabelingGraph::UNDEFINED);
    
    // fill left_chart field
    left_chart = facet2chart[initial_halfedge.facet]; // link the current boundary to the chart at its left
    charts[left_chart].boundaries.emplace(index_of_self); // link the left chart to the current boundary
    halfedge2boundary[initial_halfedge] = {index_of_self,true}; // mark this halfedge, link it to the current boundary
    halfedges.push_back(initial_halfedge); // start the list of halfedges composing the current boundary

    // fill right_chart field
    mesh_halfedges.move_to_opposite(current_halfedge); // switch orientation
    right_chart = facet2chart[current_halfedge.facet]; // link the current boundary to the chart at its right
    charts[right_chart].boundaries.emplace(index_of_self); // link the left chart to the current boundary
    halfedge2boundary[current_halfedge] = {index_of_self,false}; // mark this halfedge, link it to the current boundary (opposite orientation)
    
    // fill axis field
    axis = other_axis(
        label2axis(charts[left_chart].label),
        label2axis(charts[right_chart].label)
    );

    // Step 2 : go from boundary edge to boundary edge until we found a corner

    do {
        current_vertex = mesh.facet_corners.vertex(current_halfedge.corner); // get the vertex at the base of halfedge
        VertexRingWithBoundaries current_vertex_ring;
        current_vertex_ring.explore(current_halfedge,mesh_halfedges); // explore the halfedges around current_vertex
        if(current_vertex_ring.valence() < 3) {
            // Not a corner. At least, from the point of view of this vertex ring

            mesh_halfedges.custom_move_to_next_around_border(current_halfedge); // cross current_vertex

            if(VECTOR_CONTAINS(halfedges,current_halfedge)) { // prevent infinite loops
                fmt::println(Logger::err("LabelingGraph"),"backtracked while exploring a boundary, no valence>3 vertex ring found :(");
                geo_abort();
            }

            halfedge2boundary[current_halfedge] = {index_of_self,true}; // mark this halfedge, link it to the current boundary
            halfedges.push_back(current_halfedge); // append to the list of halfedges composing the current boundary

            mesh_halfedges.move_to_opposite(current_halfedge); // switch orientation
            halfedge2boundary[current_halfedge] = {index_of_self,false}; // mark this halfedge, link it to the current boundary (opposite orientation)

            // the vertex at the base of halfedge will be explored in the next interation
        }
        else {
            // else : we found a corner !
            end_corner = vertex2corner[current_vertex];
            if(end_corner != LabelingGraph::UNDEFINED) {
                // this vertex has already been explored, maybe not by this vertex ring
                if(!corners[end_corner].halfedge_is_in_boundary_edges(current_vertex_ring.boundary_edges.front())) {
                    // this vertex has already been explored, but by another vertex ring
                    boundary_edges_to_explore += current_vertex_ring.boundary_edges;
                    corners[end_corner].vertex_rings_with_boundaries.push_back(current_vertex_ring);
                }
                // else : this vertex has already been explored by the same vertex ring
            }
            else {
                // we found an unexplored corner
                corners.push_back(Corner()); // create a new Corner
                corners.back().vertex = current_vertex; // link it to the current vertex
                corners.back().vertex_rings_with_boundaries.push_back(current_vertex_ring); // add it the just-explored vertex ring
                end_corner = (index_t) index_of_last(corners); // link the boundary to this corner
                vertex2corner[current_vertex] = end_corner; // link the current vertex to this corner
                boundary_edges_to_explore += current_vertex_ring.boundary_edges; // the boundary edges around this corner must be explored
            }
        }
    } while (end_corner == LabelingGraph::UNDEFINED); // continue the exploration until a corner is found
}

std::ostream& operator<< (std::ostream &out, const Boundary& data) {
    fmt::println(out,"\taxis : {}",data.axis);
    fmt::println(out,"\thalfedges : {}",data.halfedges);
    fmt::println(out,"\tleft_chart : {}",OPTIONAL_TO_STRING(data.left_chart));
    fmt::println(out,"\tright_chart : {}",OPTIONAL_TO_STRING(data.right_chart));
    fmt::println(out,"\tstart_corner : {}",OPTIONAL_TO_STRING(data.start_corner));
    fmt::println(out,"\tend_corner : {}",OPTIONAL_TO_STRING(data.end_corner));
    return out;
}

void StaticLabelingGraph::fill_from(Mesh& mesh, std::string facet_attribute) {

    // based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L795
    // but here we use a disjoint-set for step 1

    // expecting a surface triangle mesh
    geo_assert(mesh.facets.are_simplices());
    geo_assert(mesh.cells.nb()==0);

    clear();

    facet2chart.resize(mesh.facets.nb()); // important: memory allocation allowing to call ds.getSetsMap() on the underlying array
    vertex2corner.resize(mesh.vertices.nb(),LabelingGraph::UNDEFINED); // contrary to facet2chart where all facets are associated to a chart, not all vertices are associated to a corner -> use UNDEFINED for vertices that are not on a corner
    CustomMeshHalfedges mesh_half_edges_(mesh); // Half edges API

    // STEP 1 : Aggregete adjacent triangles of same label as chart

    Attribute<index_t> labeling(mesh.facets.attributes(),facet_attribute); // read and store the facets attribute corresponding to the labeling

    DisjointSet<index_t> ds(labeling.size()); // create a disjoint-set data structure to group facets by labels
    for(index_t f: mesh.facets) { // for each facet of the mesh
        for(index_t le = 0; le < 3; ++le) { // for each local edge of facet f in {0,1,2}
            index_t adjacent_facet = mesh.facets.adjacent(f,le); // get the facet id of the facet beyond le
            if(labeling[f] == labeling[adjacent_facet]) { // if facets f and adjacent_facet have the same label
                ds.mergeSets(f,adjacent_facet); // merge facets f and adjacent_facet
            }
        }
    }
    std::size_t nb_charts = ds.getSetsMap(facet2chart.data()); // get the map (facet id -> chart id) and the number of charts
    charts.resize(nb_charts);

    // fill the Chart objects
    for(index_t f: mesh.facets) { // for each facet of the mesh
        Chart& current_chart = charts[facet2chart[f]]; // get the chart associated to this facet
        current_chart.label = labeling[f]; // (re)define its label
        current_chart.facets.emplace(f); // register the facet
    }

    // STEP 2 : Explore corners and boundaries
    // Ideally, we have to iterate over all vertices, and,
    // for each of them, explore adjacent facets to compute valence (number of boundary edges).
    // But I don't know how to get the adjacent facets of a given vertex.
    // What I can do, with half-edges, is to iterate over all facet corners, get the vertex at this corner,
    // and go around the vertex ring with CustomMeshHalfedges::move_to_next_around_vertex(). See VertexRingWithBoundaries::explore()
    // This works for shapes having several solids connected by a vertex only

    mesh_half_edges_.set_use_facet_region(facet_attribute); // indicate to Geogram the attribute with the charts (= regions for Geogram) from which we want to find the boundaries (= borders for Geogram)
    std::vector<CustomMeshHalfedges::Halfedge> boundary_edges_to_explore; // vector of boundary edges we encountered, and that must be explored later

    index_t current_vertex = LabelingGraph::UNDEFINED;
    for(index_t f: mesh.facets) { for(index_t c: mesh.facets.corners(f)) { // for each facet corner (f,c)
        
        CustomMeshHalfedges::Halfedge halfedge(f,c); // halfedge on facet f having corner c as base
        current_vertex = mesh.facet_corners.vertex(halfedge.corner); // get the vertex at the base of halfedge
        
        // (re)explore this vertex, because maybe it was explored from another solid of the same shape (multiple vertex rings)
        VertexRingWithBoundaries current_vertex_ring;
        current_vertex_ring.explore(halfedge,mesh_half_edges_);
        if(current_vertex_ring.valence() < 3) {
            // Not a corner. At least, from the point of view of this vertex ring
            continue;            
        }
        // else : we found a corner !

        if(vertex2corner[current_vertex] != LabelingGraph::UNDEFINED) { // if a corner already exists on this vertex
            if(corners[vertex2corner[current_vertex]].halfedge_is_in_boundary_edges(halfedge)) { // if the current halfedge is already referenced in this corner
                continue; // nothing to do, skip to next vertex
            }
            else {
                // else : the corner on this vertex has never seen the current halfedge : add the boundary edges encountered to the boundary edges to explore
                boundary_edges_to_explore += current_vertex_ring.boundary_edges;
            }
        }
        else {
            // we found an unexplored corner
            corners.push_back(Corner()); // create a new Corner
            corners.back().vertex = current_vertex; // link it to the current vertex
            corners.back().vertex_rings_with_boundaries.push_back(current_vertex_ring); // add it the just-explored vertex ring
            vertex2corner[current_vertex] = (index_t) index_of_last(corners); // link the boundary to this corner
            boundary_edges_to_explore += current_vertex_ring.boundary_edges; // the boundary edges around this corner must be explored
        }        
        
        while(!boundary_edges_to_explore.empty()) { // while there still are boundaries to explore

            halfedge = boundary_edges_to_explore.back(); // pick an halfedge from the vector
            boundary_edges_to_explore.pop_back(); // remove this halfedge from the vector

            // Check if this halfedge has already been explored since its insertion in the vector
            if(MAP_CONTAINS(halfedge2boundary,halfedge)) {
                // this halfedge is already linked to a boundary, no need to re-explore the boundary
                continue;
            }

            boundaries.push_back(Boundary()); // create a new boundary
            // explore it, edge by edge
            boundaries.back().explore(halfedge,
                                       mesh_half_edges_,
                                       (index_t) index_of_last(boundaries), // get the index of this boundary
                                       facet2chart,
                                       vertex2corner,
                                       charts,
                                       corners,
                                       halfedge2boundary,
                                       boundary_edges_to_explore);
        }
    }}

    // STEP 3 : Find boundaries with no corners
    // TODO
}

void StaticLabelingGraph::clear() {
    charts.clear();
    boundaries.clear();
    corners.clear();
    facet2chart.clear();
    halfedge2boundary.clear();
    vertex2corner.clear();
}

std::size_t StaticLabelingGraph::nb_charts() const {
    return charts.size();
}

std::size_t StaticLabelingGraph::nb_boundaries() const {
    return boundaries.size();
}

std::size_t StaticLabelingGraph::nb_corners() const {
    return corners.size();
}

std::size_t StaticLabelingGraph::nb_facets() const {
    return facet2chart.size();
}

std::size_t StaticLabelingGraph::nb_vertices() const {
    return vertex2corner.size();
}

void StaticLabelingGraph::dump_to_file(const char* filename) const {
    auto out = fmt::output_file(filename);
    out.print("{}",(*this));
}

std::ostream& operator<< (std::ostream &out, const StaticLabelingGraph& data) {

    // write charts

    for(std::size_t chart_index = 0; chart_index < data.nb_charts(); ++chart_index) {
        fmt::println(out,"chart[{}]",chart_index);
        fmt::println(out,"{}",data.charts[chart_index]);
    }
    if(data.nb_charts()==0) fmt::println(out,"no charts");

    // write boundaries

    for(std::size_t boundary_index = 0; boundary_index < data.nb_boundaries(); ++boundary_index) {
        fmt::println(out,"boundaries[{}]",boundary_index);
        fmt::println(out,"{}",data.boundaries[boundary_index]);
    }
    if(data.nb_boundaries()==0) fmt::println(out,"no boundaries");

    // write corners

    for(std::size_t corner_index = 0; corner_index < data.nb_corners(); ++corner_index) {
        fmt::println(out,"corners[{}]",corner_index);
        fmt::println(out,"{}",data.corners[corner_index]);
    }
    if(data.nb_corners()==0) fmt::println(out,"no corners");

    // write facet2chart
    
    fmt::println(out,"facet2chart");
    for(std::size_t f = 0; f < data.nb_facets(); ++f) {
        fmt::println(out,"\t[{}] {}",f,data.facet2chart[f]);
    }

    // write vertex2corner
    
    fmt::println(out,"vertex2corner");
    for(std::size_t v = 0; v < data.nb_vertices(); ++v) {
        fmt::println(out,"\t[{}] {}",v,OPTIONAL_TO_STRING(data.vertex2corner[v]));
    }
    return out;
}