// inspired by https://github.com/fprotais/robustPolycube/blob/main/lib/utils/trace.h

#pragma once

#include <geogram/mesh/mesh.h>

#include <string>
#include <set>
#include <map>

#include "LabelingGraph.h"          // for Boundary
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

using namespace GEO;

bool dump_vertex(std::string filename, const Mesh& mesh, index_t vertex_index);

bool dump_facets(std::string filename, const Mesh& mesh, std::set<index_t> facet_indices);

//idem but also write a facet attribute
template <typename T>
bool dump_facets(std::string filename, std::string attribute_name, const Mesh& mesh, std::map<index_t,T> facet_indices_and_attributes) {
    Mesh out;
    out.copy(mesh,false,MESH_VERTICES); // keep only vertices
    index_t facet_index = out.facets.create_triangles( (index_t) facet_indices_and_attributes.size());
    Attribute<T> attr(out.facets.attributes(),attribute_name);
    for(auto kv : facet_indices_and_attributes) { // for each key-value pair (facet index -> attribute value)
        out.facets.set_vertex(facet_index,0,mesh.facet_corners.vertex(mesh.facets.corner(kv.first,0)));
        out.facets.set_vertex(facet_index,1,mesh.facet_corners.vertex(mesh.facets.corner(kv.first,1)));
        out.facets.set_vertex(facet_index,2,mesh.facet_corners.vertex(mesh.facets.corner(kv.first,2)));
        attr[facet_index] = kv.second;
        facet_index++;
    }
    out.vertices.remove_isolated();
    mesh_save(out,filename + ".geogram");
}

bool dump_all_boundaries(std::string filename, const Mesh& mesh, const CustomMeshHalfedges& mesh_he, const std::vector<Boundary>& boundaries);