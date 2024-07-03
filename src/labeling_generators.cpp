#include <geogram/basic/geometry.h> // for GEO::mat3
#include <geogram/mesh/mesh.h>      // for GEO::Mesh

#include <random>           // for std::random_device, std::uniform_int_distribution<>

#include "labeling_generators.h"
#include "labeling.h"
#include "containers_std.h" // for key_at_max_value(), fill_set_with_map_keys()
#include "GraphCutLabeling.h"

void random_labeling(GEO::Mesh& mesh, const char* attribute_name) {
    std::random_device rd; // https://en.cppreference.com/w/cpp/numeric/random/random_device
    std::uniform_int_distribution<index_t> dist(0, 5); // 0 and 5 included
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    for(index_t f: mesh.facets) { // for each facet
        label[f] = dist(rd);
    }
}

void naive_labeling(Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name) {

    // use GEO::Geom::triangle_normal_axis() instead ?

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    for(index_t f: mesh.facets) { // for each facet
        label[f] = nearest_label(normals[f]);
    }
}

void tweaked_naive_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name) {

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    std::set<index_t> facets_to_tilt;
    get_facets_to_tilt(mesh,normals,facets_to_tilt,(double) NAIVE_LABELING_TWEAK_SENSITIVITY);
    FOR(f,mesh.facets.nb()) { // for each facet
        label[f] = facets_to_tilt.contains(f) ? nearest_label(mult(rotation,normals[f])) : nearest_label(normals[f]);
    }
}

void smart_init_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    std::set<index_t> facets_to_tilt;
    get_facets_to_tilt(mesh,normals,facets_to_tilt,(double) NAIVE_LABELING_TWEAK_SENSITIVITY);
    std::vector<index_t> per_facet_group_index(mesh.facets.nb());
    index_t nb_groups = group_facets<std::vector>(mesh,facets_to_tilt,per_facet_group_index);
    // The group n°0 is facets we don't have to tilt
    // For groups we have to tilt, compute total surface area
    std::map<index_t,double> per_group_area;
    index_t current_facet_group_index = index_t(-1);
    double current_facet_area = 0.0;
    FOR(f,mesh.facets.nb()) {
        current_facet_group_index = per_facet_group_index[f];
        if(current_facet_group_index == 0) {
            continue; // do not compute the area of group n°0
        }
        current_facet_area = mesh_facet_area(mesh,f);
        per_group_area[current_facet_group_index] = (per_group_area.contains(current_facet_group_index) ? per_group_area[current_facet_group_index] : 0.0) + current_facet_area;
    }
    geo_assert(per_group_area.size() == nb_groups-1);

    double total_surface_area = mesh_area(mesh);
    index_t group_having_largest_area = key_at_max_value(per_group_area);

    std::set<index_t> groups_surrounded_by_feature_edges;
    fill_set_with_map_keys(per_group_area,groups_surrounded_by_feature_edges);
    index_t adjacent_facet = index_t(-1);
    FOR(f,mesh.facets.nb()) {
        FOR(le,3) {
            adjacent_facet = mesh.facets.adjacent(f,le);
            // local edge k is the one between local vertices k and (k+1)%3
            // see https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
            if(
                (!local_edge_is_on_feature_edge(mesh,f,le,feature_edges)) &&
                (per_facet_group_index[f] != per_facet_group_index[adjacent_facet])
            ) {
                groups_surrounded_by_feature_edges.erase(per_facet_group_index[f]);
                groups_surrounded_by_feature_edges.erase(per_facet_group_index[adjacent_facet]);
                if(groups_surrounded_by_feature_edges.empty()) {
                    goto break_both_loops;
                }
            }
        }
    }

break_both_loops:

    if(
        (per_group_area[group_having_largest_area] > total_surface_area * 0.01) ||
        !groups_surrounded_by_feature_edges.empty()
    ) {
        // the group having the largest area is bigger than 1% of the total surface
        fmt::println(Logger::out("labeling"),"Init labeling: *tweaked* naive labeling"); Logger::out("labeling").flush();
        tweaked_naive_labeling(mesh,normals,attribute_name); // TODO transmit `per_facet_group_index` to not recompute it
    }
    else {
        fmt::println(Logger::out("labeling"),"Init labeling: naive labeling"); Logger::out("labeling").flush();
        naive_labeling(mesh,normals,attribute_name);
    }
}

// https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L102
void graphcut_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, int compactness_coeff, int fidelity_coeff) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);
    auto gcl = GraphCutLabeling(mesh,normals);
    gcl.data_cost__set__fidelity_based(fidelity_coeff);
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(compactness_coeff);
    gcl.compute_solution(label);
}

void compute_per_facet_fidelity(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* labeling_attribute_name, const char* fidelity_attribute_name, BasicStats& stats) {
    // goal of the fidelity metric : measure how close the assigned direction (label) is from the triangle normal
    //   fidelity=1 (high fidelity) -> label equal to the normal      ex: triangle is oriented towards +X and we assign +X
    //   fidelity=0 (low fidelity)  -> label opposite to the normal   ex: triangle is oriented towards +X and we assign -X
    // 1. compute and normalize the normal
    // 2. compute the vector (ex: 0.0,1.0,0.0) of the label (ex: +Y) with label2vector[] (already normalized)
    // 3. compute the dot product, which is in [-1:1] : 1=equal, -1=opposite
    // 4. add 1 and divide by 2 so that the range is [0:1] : 1=equal, 0=opposite
    Attribute<index_t> label(mesh.facets.attributes(), labeling_attribute_name);
    Attribute<double> per_facet_fidelity(mesh.facets.attributes(), fidelity_attribute_name);
    stats.reset();
    FOR(f,mesh.facets.nb()) {
        per_facet_fidelity[f] = (GEO::dot(normals[f],label2vector[label[f]]) + 1.0)/2.0;
        stats.insert(per_facet_fidelity[f]);
    }
}