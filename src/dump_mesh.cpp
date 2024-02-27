#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/attributes.h>

#include <string>
#include <vector>

#include "dump_mesh.h"
#include "LabelingGraph.h"          // for Boundary
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

bool dump_vertex(std::string filename, vec3 coordinates) {
    Mesh out;
    out.vertices.create_vertices(1);
    out.vertices.point(0) = coordinates;
    return mesh_save(out,filename + ".geogram");
}

bool dump_vector(std::string filename, const vec3& origin, const vec3& vector) {
    Mesh out;
    index_t v0 = out.vertices.create_vertex(origin.data());
    geo_assert(v0 == 0);
    vec3 tip_of_vector = origin+vector;
    index_t v1 = out.vertices.create_vertex(tip_of_vector.data());
    geo_assert(v1 == 1);
    out.edges.create_edge(v0,v1);
    Attribute<bool> attr_index(out.vertices.attributes(),"index_in_edge");
    attr_index[v0] = 0; // 0 = origin vertex
    attr_index[v1] = 1; // 1 = extremity vertex
    return mesh_save(out,filename + ".geogram");
}

bool dump_edge(std::string filename, const Mesh& mesh, MeshHalfedges::Halfedge& halfedge) {
    Mesh out;
    out.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t edge_index = out.edges.create_edge();
    out.edges.set_vertex(edge_index,0,Geom::halfedge_vertex_index_from(mesh,halfedge));
    out.edges.set_vertex(edge_index,1,Geom::halfedge_vertex_index_to(mesh,halfedge));
    out.vertices.remove_isolated();
    return mesh_save(out,filename + ".geogram");
}

bool dump_edges(std::string filename, const Mesh& mesh, const std::set<std::pair<index_t,index_t>>& edges) {
    Mesh out;
    out.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t edge_index = out.edges.create_edges( (index_t) edges.size());
    for(const auto& e : edges) {
        out.edges.set_vertex(edge_index,0,e.first);
        out.edges.set_vertex(edge_index,1,e.second);
        edge_index++;
    }
    out.vertices.remove_isolated();
    return mesh_save(out,filename + ".geogram");
}

bool dump_boundary(std::string filename, const Mesh& mesh, const Boundary& boundary) {
    Mesh boundary_mesh;
    boundary_mesh.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t first_edge_index = boundary_mesh.edges.create_edges((index_t) boundary.halfedges.size()); // create as many edges as they are halfedges in `boundary`
    geo_assert(first_edge_index==0);
    FOR(halfedge_index,boundary.halfedges.size()) {
        boundary_mesh.edges.set_vertex(halfedge_index,0,Geom::halfedge_vertex_index_from(mesh,boundary.halfedges[halfedge_index]));
        boundary_mesh.edges.set_vertex(halfedge_index,1,Geom::halfedge_vertex_index_to(mesh,boundary.halfedges[halfedge_index]));
    }
    boundary_mesh.vertices.remove_isolated();
    return mesh_save(boundary_mesh,filename + ".geogram");
}

bool dump_boundary_with_halfedges_indices(std::string filename, const Mesh& mesh, const Boundary& boundary) {
    Mesh boundary_mesh;
    boundary_mesh.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t first_edge_index = boundary_mesh.edges.create_edges((index_t) boundary.halfedges.size()); // create as many edges as they are halfedges in `boundary`
    geo_assert(first_edge_index==0);
    Attribute<index_t> attr_index(boundary_mesh.edges.attributes(),"halfedge_index");
    FOR(halfedge_index,boundary.halfedges.size()) {
        boundary_mesh.edges.set_vertex(halfedge_index,0,Geom::halfedge_vertex_index_from(mesh,boundary.halfedges[halfedge_index]));
        boundary_mesh.edges.set_vertex(halfedge_index,1,Geom::halfedge_vertex_index_to(mesh,boundary.halfedges[halfedge_index]));
        attr_index[halfedge_index] = halfedge_index;
    }
    boundary_mesh.vertices.remove_isolated();
    return mesh_save(boundary_mesh,filename + ".geogram");
}

bool dump_chart_contour(std::string filename, const CustomMeshHalfedges& mesh_he, const StaticLabelingGraph& slg, index_t chart_index) {
    geo_assert(chart_index < slg.nb_charts());
    geo_assert(mesh_he.is_using_facet_region());
    const Mesh& mesh = mesh_he.mesh();
    std::vector<std::pair<index_t,bool>> counterclockwise_order;
    slg.charts[chart_index].counterclockwise_boundaries_order(mesh_he,slg.halfedge2boundary,slg.boundaries,counterclockwise_order);
    geo_assert(!counterclockwise_order.empty());
    std::map<std::pair<index_t,index_t>,index_t> counterclockwise_halfedges_around_chart;
    index_t halfedge_count = 0;
    for(const auto& [b,same_order] : counterclockwise_order) {
        const Boundary& boundary = slg.boundaries[b];
        if(same_order) {
            for(auto it = boundary.halfedges.begin(); it != boundary.halfedges.end(); ++it) {
                counterclockwise_halfedges_around_chart[
                    std::make_pair(
                        halfedge_vertex_index_from(mesh,*it),
                        halfedge_vertex_index_to(mesh,*it)
                    )
                ] = halfedge_count;
                halfedge_count++;
            }
        }
        else {
            // same but use reverse iterator
            for(auto it = boundary.halfedges.rbegin(); it != boundary.halfedges.rend(); ++it) {
                counterclockwise_halfedges_around_chart[
                    std::make_pair(
                        halfedge_vertex_index_from(mesh,*it),
                        halfedge_vertex_index_to(mesh,*it)
                    )
                ] = halfedge_count;
                halfedge_count++;
            }
        }
    }
    return dump_edges(filename,"counterclockwise_order",mesh,counterclockwise_halfedges_around_chart);
}

bool dump_all_boundaries(std::string filename, const Mesh& mesh, const std::vector<Boundary>& boundaries) {
    Mesh boundaries_mesh;
    boundaries_mesh.copy(mesh,false,MESH_VERTICES); // keep only vertices
    for(const Boundary& b : boundaries) { // for each boundary of the labeling graph
        index_t edge_index = boundaries_mesh.edges.create_edges((index_t) b.halfedges.size()); // create as many edges as they are halfedges in the current boundary
        for(const auto& he : b.halfedges) { // for each halfedge of the current boundary
            boundaries_mesh.edges.set_vertex(edge_index,0,Geom::halfedge_vertex_index_from(mesh,he)); // get the vertex at the base of 'halfedge', set as 1st vertex
            boundaries_mesh.edges.set_vertex(edge_index,1,Geom::halfedge_vertex_index_to(mesh,he)); // get the vertex at the base of 'halfedge', set as 2nd vertex
            edge_index++;
        }
    }
    boundaries_mesh.vertices.remove_isolated();
    return mesh_save(boundaries_mesh,filename + ".geogram");
}

bool dump_all_boundaries_with_indices_and_axes(std::string filename, const Mesh& mesh, const StaticLabelingGraph& slg) {
    // index_t edges_counter = 0;
    Mesh boundaries_mesh;
    boundaries_mesh.copy(mesh,false,MESH_VERTICES); // keep only vertices
    Attribute<index_t> attr_index(boundaries_mesh.edges.attributes(),"index");
    Attribute<int> attr_axis(boundaries_mesh.edges.attributes(),"axis");
    FOR(boundary_index,slg.boundaries.size()) { // for each boundary of the labeling graph
        const Boundary& b = slg.boundaries[boundary_index];
        index_t edge_index = boundaries_mesh.edges.create_edges((index_t) b.halfedges.size()); // create as many edges as they are halfedges in the current boundary
        for(const auto& he : b.halfedges) { // for each halfedge of the current boundary
            boundaries_mesh.edges.set_vertex(edge_index,0,Geom::halfedge_vertex_index_from(mesh,he)); // get the vertex at the base of 'halfedge', set as 1st vertex
            boundaries_mesh.edges.set_vertex(edge_index,1,Geom::halfedge_vertex_index_to(mesh,he)); // get the vertex at the base of 'halfedge', set as 2nd vertex
            attr_index[edge_index] = boundary_index;
            attr_axis[edge_index] = b.axis;
            edge_index++;
        }
    }
    boundaries_mesh.vertices.remove_isolated();
    return mesh_save(boundaries_mesh,filename + ".geogram");
}

bool dump_facets(std::string filename, const Mesh& mesh, std::set<index_t> facet_indices) {
    Mesh out;
    out.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t facet_index = out.facets.create_triangles( (index_t) facet_indices.size());
    for(index_t f : facet_indices) {
        out.facets.set_vertex(facet_index,0,mesh.facet_corners.vertex(mesh.facets.corner(f,0)));
        out.facets.set_vertex(facet_index,1,mesh.facet_corners.vertex(mesh.facets.corner(f,1)));
        out.facets.set_vertex(facet_index,2,mesh.facet_corners.vertex(mesh.facets.corner(f,2)));
        facet_index++;
    }
    out.vertices.remove_isolated();
    return mesh_save(out,filename + ".geogram");
}