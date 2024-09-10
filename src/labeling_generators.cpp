#include <geogram/basic/geometry.h> // for GEO::mat3
#include <geogram/mesh/mesh.h>      // for GEO::Mesh

#include <random>           // for std::random_device, std::uniform_int_distribution<>
#include <cmath>            // for std::fpclassify()

#include "labeling_generators.h"
#include "labeling.h"
#include "containers_std.h" // for key_at_max_value(), fill_set_with_map_keys()
#include "labeling_graphcuts.h"

void random_labeling(const Mesh& mesh, Attribute<index_t>& labeling) {
    geo_debug_assert(labeling.is_bound());
    std::random_device rd; // https://en.cppreference.com/w/cpp/numeric/random/random_device
    std::uniform_int_distribution<index_t> dist(0, 5); // 0 and 5 included
    for(index_t f: mesh.facets) { // for each facet
        labeling[f] = dist(rd);
    }
}

void naive_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling) {
    geo_debug_assert(labeling.is_bound());
    // use GEO::Geom::triangle_normal_axis() instead ?
    for(index_t f: mesh_ext.facets) { // for each facet
        labeling[f] = nearest_label(mesh_ext.facet_normals[f]);
    }
}

void tweaked_naive_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling) {
    geo_debug_assert(labeling.is_bound());
    std::set<index_t> facets_to_tilt;
    get_facets_to_tilt(mesh_ext,facets_to_tilt,(double) NAIVE_LABELING_TWEAK_SENSITIVITY);
    FOR(f,mesh_ext.facets.nb()) { // for each facet
        labeling[f] = facets_to_tilt.contains(f) ? nearest_label(mult(ROTATION_MATRIX,mesh_ext.facet_normals[f])) : nearest_label(mesh_ext.facet_normals[f]);
    }
}

void smart_init_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling) {
    geo_debug_assert(labeling.is_bound());
    std::set<index_t> facets_to_tilt;
    get_facets_to_tilt(mesh_ext,facets_to_tilt,(double) NAIVE_LABELING_TWEAK_SENSITIVITY);
    std::vector<index_t> per_facet_group_index(mesh_ext.facets.nb());
    index_t nb_groups = group_facets<std::vector>(mesh_ext,facets_to_tilt,per_facet_group_index);
    // The group n°0 is facets we don't have to tilt
    // For groups we have to tilt, compute total surface area
    std::map<index_t,double> per_group_area;
    index_t current_facet_group_index = index_t(-1);
    double current_facet_area = 0.0;
    FOR(f,mesh_ext.facets.nb()) {
        current_facet_group_index = per_facet_group_index[f];
        if(current_facet_group_index == 0) {
            continue; // do not compute the area of group n°0
        }
        current_facet_area = mesh_facet_area(mesh_ext,f);
        per_group_area[current_facet_group_index] = (per_group_area.contains(current_facet_group_index) ? per_group_area[current_facet_group_index] : 0.0) + current_facet_area;
    }
    geo_assert(per_group_area.size() == nb_groups-1);

    double total_surface_area = mesh_area(mesh_ext);
    index_t group_having_largest_area = key_at_max_value(per_group_area);

    std::set<index_t> groups_surrounded_by_feature_edges;
    fill_set_with_map_keys(per_group_area,groups_surrounded_by_feature_edges);
    index_t adjacent_facet = index_t(-1);
    FOR(f,mesh_ext.facets.nb()) {
        FOR(le,3) {
            adjacent_facet = mesh_ext.facets.adjacent(f,le);
            // local edge k is the one between local vertices k and (k+1)%3
            // see https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
            if(
                (!mesh_ext.feature_edges.contain_facet_edge(f,le)) &&
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
        tweaked_naive_labeling(mesh_ext,labeling); // TODO transmit `per_facet_group_index` to not recompute it
    }
    else {
        fmt::println(Logger::out("labeling"),"Init labeling: naive labeling"); Logger::out("labeling").flush();
        naive_labeling(mesh_ext,labeling);
    }
}

// https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L102
void graphcut_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling, int compactness_coeff, int fidelity_coeff, double sensitivity, double angle_of_rotation) {
    geo_debug_assert(labeling.is_bound());
    auto gcl = GraphCutLabeling(mesh_ext);
    // TODO remove redundancy with app/graphcut_labeling.cpp, under ImGui::Button("Compute solution")
    if(std::fpclassify(sensitivity) == FP_ZERO) {
        // sensitivity is null, use the fidelity based data cost
        gcl.data_cost__set__fidelity_based(fidelity_coeff);
    }
    else {
        std::vector<int> custom_data_cost(mesh_ext.facets.nb()*6);
        vec3 normal;
        mat3 rotation_to_apply = rotation_matrix(angle_of_rotation);
        FOR(f,mesh_ext.facets.nb()) {
            normal = mesh_ext.facet_normals[f];
            if(is_a_facet_to_tilt(normal,sensitivity)) {
                normal = mult(rotation_to_apply,normal); // rotation of the normal
            }
            FOR(label,6) {
                double dot = (GEO::dot(normal,label2vector[label]) - 1.0)/0.2;
                double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
                custom_data_cost[f*6+label] = (int) (fidelity_coeff*100*cost);
            }
        }
        gcl.data_cost__set__all_at_once(custom_data_cost);
    }
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(compactness_coeff);
    gcl.compute_solution(labeling);
}