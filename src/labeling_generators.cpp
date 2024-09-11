#include <geogram/basic/geometry.h> // for GEO::mat3
#include <geogram/mesh/mesh.h>      // for GEO::Mesh

#include <random>           // for std::random_device, std::uniform_int_distribution<>
#include <cmath>            // for std::fpclassify()

#include "labeling_generators.h"
#include "labeling.h"
#include "containers_std.h" // for key_at_max_value(), fill_set_with_map_keys()
#include "labeling_graphcuts.h"
#include "geometry.h" // for rotation_matrix()

void random_labeling(const Mesh& mesh, Attribute<index_t>& labeling) {
    geo_debug_assert(labeling.is_bound());
    std::random_device rd; // https://en.cppreference.com/w/cpp/numeric/random/random_device
    std::uniform_int_distribution<index_t> dist(0, 5); // 0 and 5 included
    for(index_t f: mesh.facets) { // for each facet
        labeling[f] = dist(rd);
    }
}

void naive_labeling(const MeshExt& mesh_ext, Attribute<index_t>& labeling, double sensitivity, double angle_of_rotation) {
    geo_debug_assert(labeling.is_bound());
    // use GEO::Geom::triangle_normal_axis() instead ?

    std::set<index_t> facets_to_tilt;
    mat3 rotation_to_apply;
    if(std::fpclassify(sensitivity) != FP_ZERO) {
        get_facets_to_tilt(mesh_ext,facets_to_tilt,sensitivity);
        rotation_to_apply = rotation_matrix(angle_of_rotation);
    }
    // else: leave `facets_to_tilt` empty
    
    for(index_t f: mesh_ext.facets) { // for each facet
        labeling[f] = facets_to_tilt.contains(f) ? nearest_label(mult(rotation_to_apply,mesh_ext.facet_normals[f])) : nearest_label(mesh_ext.facet_normals[f]);
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