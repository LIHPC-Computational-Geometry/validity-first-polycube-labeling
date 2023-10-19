#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/attributes.h>

#include <string>
#include <vector>

#include "dump_mesh.h"
#include "LabelingGraph.h"          // for Boundary
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

bool dump_vertex(std::string filename, const Mesh& mesh, index_t vertex_index) {
    Mesh out;
    out.vertices.create_vertices(1);
    out.vertices.point(0) = mesh.vertices.point(vertex_index);
    return mesh_save(out,filename + ".geogram");
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
    mesh_save(out,filename + ".geogram");
}

bool dump_all_boundaries(std::string filename, const Mesh& mesh, const CustomMeshHalfedges& mesh_he, const std::vector<Boundary>& boundaries) {
    Mesh boundaries_mesh;
    boundaries_mesh.copy(mesh,false,MESH_VERTICES); // keep only vertices
    for(const Boundary& b : boundaries) { // for each boundary of the labeling graph
        index_t edge_index = boundaries_mesh.edges.create_edges(b.halfedges.size()); // create as many edges as they are halfedges in the current boundary
        for(const auto& he : b.halfedges) { // for each halfedge of the current boundary
            CustomMeshHalfedges::Halfedge halfedge = he; //create a mutable copy
            boundaries_mesh.edges.set_vertex(edge_index,0,mesh.facet_corners.vertex(halfedge.corner)); // get the vertex at the base of 'halfedge', set as 1st vertex
            mesh_he.move_to_opposite(halfedge); // switch to the opposite halfedge
            boundaries_mesh.edges.set_vertex(edge_index,1,mesh.facet_corners.vertex(halfedge.corner)); // get the vertex at the base of 'halfedge', set as 2nd vertex
            edge_index++;
        }
    }
    boundaries_mesh.vertices.remove_isolated();
    return mesh_save(boundaries_mesh,filename + ".geogram");
}