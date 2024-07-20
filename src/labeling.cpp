#include <geogram/basic/attributes.h>   // for GEO::Attribute

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <ranges>           // for std::ranges::view::keys()

#include "labeling.h"
#include "labeling_graph.h"
#include "containers_macros.h" // for VECTOR_CONTAINS(), MAP_CONTAINS()
#include "containers_std.h"    // for key_at_max_value()
#include "containers_Geogram.h"
#include "geometry_halfedges.h"


index_t nearest_label(const vec3& normal) {
    // from 3 signed to 6 unsigned components:
    std::array<double,6> weights = {
        normal.x < 0.0 ? 0.0 : normal.x,  // +X
        normal.x < 0.0 ? -normal.x : 0.0, // -X
        normal.y < 0.0 ? 0.0 : normal.y,  // +Y
        normal.y < 0.0 ? -normal.y : 0.0, // -Y
        normal.z < 0.0 ? 0.0 : normal.z,  // +Z
        normal.z < 0.0 ? -normal.z : 0.0  // -Z
    };
    // return index of max
    return (index_t) VECTOR_MAX_INDEX(weights);
}

index_t nearest_axis(const vec3& vector) {
    // find dimension with biggest coeff
    std::array<double,3> array({
        std::abs(vector.x),
        std::abs(vector.y),
        std::abs(vector.z)});
    return (index_t) VECTOR_MAX_INDEX(array); // 0=X, 1=Y, 2=Z
}

bool is_better_label(const vec3& facet_normal, index_t current_label, index_t new_label) {
    return dot(facet_normal,label2vector[new_label]) > dot(facet_normal,label2vector[current_label]);
}

bool are_orthogonal_labels(index_t label1, index_t label2) {
    geo_assert(label1 < 6);
    geo_assert(label2 < 6);
    return (label1/2) != (label2/2); // orthogonal if they don't have the same axis
}

index_t opposite_label(index_t label) {
    geo_assert(label < 6);
    return (label % 2 == 0) ? label+1 : label-1;
}

void flip_labeling(Mesh& mesh, Attribute<index_t>& labeling) {
    geo_debug_assert(labeling.is_bound());
    FOR(f,mesh.facets.nb()) { // for each facet
        labeling[f] = opposite_label(labeling[f]);
    }
}

index_t find_optimal_label(std::initializer_list<index_t> forbidden_axes, std::initializer_list<index_t> forbidden_labels, std::initializer_list<index_t> orthogonal_labels, vec3 close_vector) {
    std::map<index_t,double> candidates;
    FOR(label,6) {
        if(INIT_LIST_CONTAINS(forbidden_axes,label/2)) {
            goto ignore_this_label; // current `label` is on a forbidden axis
        }
        if(INIT_LIST_CONTAINS(forbidden_labels,label)) {
            goto ignore_this_label; // current `label` is forbidden
        }
        for(index_t orthogonal_label : orthogonal_labels) {
            if(!are_orthogonal_labels(label,orthogonal_label)) {
                goto ignore_this_label; // current `label` doesn't comply with one of the orthogonality constraints
            }
            // else: they are orthogonal, check other orthogonality constraints
        }
        // so `label` passed all the criteria
        if(close_vector == vec3(0.0,0.0,0.0)) { // if no `close_vector` provided
            candidates[label] = 0.0; // no dot product to compute
        }
        else {
            candidates[label] = dot(label2vector[label],normalize(close_vector));
        }
    ignore_this_label:
        continue;
    }
    if(candidates.empty()) {
        fmt::println(Logger::err("labeling"),"In find_optimal_label(), no label satisfy all constraints:");
        fmt::println(Logger::err("labeling")," - not on axes {}",forbidden_axes);
        fmt::println(Logger::err("labeling")," - not {}",forbidden_labels);
        fmt::println(Logger::err("labeling")," - orthogonal to {}",orthogonal_labels);
        Logger::err("labeling").flush();
        geo_assert_not_reached;
    }
    if(close_vector == vec3(0.0,0.0,0.0)) { // if no `close_vector` provided
        if(candidates.size() > 1) {
            fmt::println(Logger::err("labeling"),"In find_optimal_label(), several labels ({}) match the constraints,",std::ranges::views::keys(candidates));
            fmt::println(Logger::err("labeling"),"but no close vector was provided to prefer one");
            Logger::err("labeling").flush();
            geo_assert_not_reached;
        }
        return candidates.begin()->first;
    }
    // so a `close_vector` was provided
    index_t optimal_label = key_at_max_value(candidates); // return label with max dot product with `close_vector`
    if(candidates[optimal_label] <= 0.0) { // the dot product should be positive, else the label is not close to the `close_vector` the user asked for
        fmt::println(Logger::warn("labeling"),"In find_optimal_label(), max dot product is negative (label = {}, dot product = {})",optimal_label,candidates[optimal_label]);
        Logger::warn("labeling").flush();
    }
    return optimal_label;
}

void propagate_label(const Mesh& mesh, Attribute<index_t>& labeling, index_t new_label, const std::set<index_t>& facets_in, const std::set<index_t> facets_out, const std::vector<index_t>& facet2chart, index_t chart_index) {
    // change label to `new_label` for all facets in `facets_in` and their neighbors, step by step.
    // 2 limits for the propagation:
    //  - stay on `chart_index`
    //  - do not cross facets in `facets_out`
    geo_debug_assert(labeling.is_bound());
    std::vector<index_t> to_process(facets_in.begin(),facets_in.end()); // set -> vector data structure
    geo_assert(!to_process.empty());

    index_t current_facet = index_t(-1),
            adjacent_facet = index_t(-1);
    while (!to_process.empty()) {
        current_facet = to_process.back();
        to_process.pop_back();
        // don't check if `current_facet` is on `chart_index`, too restrictive is some cases
        if(labeling[current_facet] == new_label) {
            // facet already processed since insertion
            continue;
        }
        labeling[current_facet] = new_label;
        FOR(le,3) { // process adjacent triangles
            adjacent_facet = mesh.facets.adjacent(current_facet,le);
            if(facets_out.contains(adjacent_facet)) { // do not cross `facets_out` (walls)
                continue;
            }
            if(facet2chart[adjacent_facet] != chart_index) { // do not go away from `chart_index`
                continue;
            }
            to_process.push_back(adjacent_facet);
        }
    }
}

void compute_per_facet_fidelity(const MeshExt& mesh_ext, Attribute<index_t>& labeling, const char* fidelity_attribute_name, IncrementalStats& stats) {
    // goal of the fidelity metric : measure how close the assigned direction (label) is from the triangle normal
    //   fidelity=1 (high fidelity) -> label equal to the normal      ex: triangle is oriented towards +X and we assign +X
    //   fidelity=0 (low fidelity)  -> label opposite to the normal   ex: triangle is oriented towards +X and we assign -X
    // 1. compute and normalize the normal
    // 2. compute the vector (ex: 0.0,1.0,0.0) of the label (ex: +Y) with label2vector[] (already normalized)
    // 3. compute the dot product, which is in [-1:1] : 1=equal, -1=opposite
    // 4. add 1 and divide by 2 so that the range is [0:1] : 1=equal, 0=opposite
    geo_debug_assert(labeling.is_bound());
    Attribute<double> per_facet_fidelity(mesh_ext.facets.attributes(), fidelity_attribute_name);
    stats.reset();
    FOR(f,mesh_ext.facets.nb()) {
        per_facet_fidelity[f] = (GEO::dot(mesh_ext.facet_normals[f],label2vector[labeling[f]]) + 1.0)/2.0;
        stats.insert(per_facet_fidelity[f]);
    }
}

unsigned int count_lost_feature_edges(const MeshExt& mesh) {
    // parse all facet
    // parse all local vertex for each facet (= facet corners)
    // get halfedge
    // ignore if v1 < v0 (they are 2 halfeges for a given edge, keep the one where v0 < v1)
    // ignore if halfedge not on feature edge
    // fetch label at left and right
    // if same labeling, increment counter

    unsigned int nb_lost_feature_edges = 0;
    geo_assert(mesh.halfedges.is_using_facet_region()); // expecting the labeling to be bounded in `mesh.halfedges`
    MeshHalfedges::Halfedge current_halfedge;
    FOR(f,mesh.facets.nb()) {
        FOR(lv,3) { // for each local vertex of the current facet
            current_halfedge.facet = f;
            current_halfedge.corner = mesh.facets.corner(f,lv);
            geo_assert(mesh.halfedges.halfedge_is_valid(current_halfedge));
            if(halfedge_vertex_index_to(mesh,current_halfedge) < halfedge_vertex_index_from(mesh,current_halfedge)) {
                continue;
            }
            if(!mesh.feature_edges.contain_halfedge(current_halfedge)) {
                continue; // not on feature edge
            }
            if(!mesh.halfedges.halfedge_is_border(current_halfedge)) {
                // so our `current_halfedge`
                // - has v0 < v1 (to count each undirected edge once)
                // - is on a feature edge
                // - has the same label at left and right
                nb_lost_feature_edges++;
            }
        }
    }
    geo_assert(nb_lost_feature_edges <= mesh.feature_edges.nb());
    return nb_lost_feature_edges;
}

index_t adjacent_chart_in_common(const Boundary& b0, const Boundary& b1) {
    if( (b0.left_chart == b1.left_chart) || (b0.left_chart == b1.right_chart) ) {
        return b0.left_chart;
    }
    else if( (b0.right_chart == b1.left_chart) || (b0.right_chart == b1.right_chart) ) {
        return b0.right_chart;
    }
    else {
        // no chart in common
        geo_assert_not_reached;
    }
}