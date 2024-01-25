// inspired by https://github.com/fprotais/robustPolycube/blob/main/lib/utils/trace.h

#pragma once

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_halfedges.h>    // for MeshHalfedges::Halfedge
#include <geogram/mesh/mesh_io.h>           // for mesh_save()

#include <string>
#include <set>
#include <map>
#include <utility>  // for std::pair

#include "LabelingGraph.h"          // for Boundary
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

using namespace GEO;

// 0D

bool dump_vertex(std::string filename, const Mesh& mesh, index_t vertex_index);

// 1D

bool dump_vector(std::string filename, const vec3& origin, const vec3& vector);

bool dump_vector(std::string filename, const Mesh& mesh, index_t origin_vertex, const vec3& vector);

bool dump_edge(std::string filename, const Mesh& mesh, MeshHalfedges::Halfedge& halfedge);

bool dump_edges(std::string filename, const Mesh& mesh, const std::set<std::pair<index_t,index_t>>& edges);

//idem but also write an edge attribute
template <typename T>
bool dump_edges(std::string filename, std::string attribute_name, const Mesh& mesh, const std::map<std::pair<index_t,index_t>,T>& edges_and_attributes) {
    Mesh out;
    out.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t edge_index = out.edges.create_edges( (index_t) edges_and_attributes.size());
    Attribute<T> attr(out.edges.attributes(),attribute_name);
    for(const auto& kv : edges_and_attributes) { // for each key-value pair (facet index -> attribute value)
        out.edges.set_vertex(edge_index,0,kv.first.first);
        out.edges.set_vertex(edge_index,1,kv.first.second);
        attr[edge_index] = kv.second;
        edge_index++;
    }
    out.vertices.remove_isolated();
    return mesh_save(out,filename + ".geogram");
}

bool dump_boundary(std::string filename, const Mesh& mesh, const Boundary& boundary);

bool dump_boundary_with_halfedges_indices(std::string filename, const Mesh& mesh, const Boundary& boundary);

bool dump_all_boundaries(std::string filename, const Mesh& mesh, const std::vector<Boundary>& boundaries);

bool dump_all_boundaries_with_indices_and_axes(std::string filename, const Mesh& mesh, const StaticLabelingGraph& slg);

// 2D

bool dump_facets(std::string filename, const Mesh& mesh, std::set<index_t> facet_indices);

//idem but also write a facet attribute
template <typename T>
bool dump_facets(std::string filename, std::string attribute_name, const Mesh& mesh, std::map<index_t,T> facet_indices_and_attributes) {
    Mesh out;
    out.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t facet_index = out.facets.create_triangles( (index_t) facet_indices_and_attributes.size());
    Attribute<T> attr(out.facets.attributes(),attribute_name);
    for(const auto& kv : facet_indices_and_attributes) { // for each key-value pair (facet index -> attribute value)
        out.facets.set_vertex(facet_index,0,mesh.facet_corners.vertex(mesh.facets.corner(kv.first,0)));
        out.facets.set_vertex(facet_index,1,mesh.facet_corners.vertex(mesh.facets.corner(kv.first,1)));
        out.facets.set_vertex(facet_index,2,mesh.facet_corners.vertex(mesh.facets.corner(kv.first,2)));
        attr[facet_index] = kv.second;
        facet_index++;
    }
    out.vertices.remove_isolated();
    return mesh_save(out,filename + ".geogram");
}