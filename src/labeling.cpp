#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/vecg.h>     // for vec3
#include <geogram/basic/matrix.h>   // for mat3
#include <geogram/basic/logger.h>   // for Logger::*

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/os.h>

#include <array>            // for std::array
#include <initializer_list> // for std::initializer_list
#include <algorithm>        // for std::max_element(), std::min(), std::max(), std::sort()
#include <iterator>         // for std::distance()
#include <tuple>            // for std::tuple
#include <queue>            // for std::queue
#include <cmath>            // for std::pow(), std::abs()
#include <set>              // for std::set
#include <map>              // for std::map()
#include <ranges>           // for std::ranges::view::keys(), std::ranges::sort()

#include "labeling.h"
#include "LabelingGraph.h"
#include "containers.h"
#include "GraphCutLabeling.h"
#include "CustomMeshHalfedges.h"
#include "basic_stats.h"
#include "dump_mesh.h"

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name) {

    //open the file
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        fmt::println(Logger::out("I/O"),"Could not open labeling file '{}'",filename); 
        return false;
    }

    fmt::println(Logger::out("I/O"),"Loading file {}...",filename); Logger::out("I/O").flush();

    std::string current_line;
    unsigned long current_label; // necessary type for string to unsigned int conversion (stoul)
    GEO::index_t current_line_number = 0;
    GEO::index_t facets_number = mesh.facets.nb(); // expected line number (one line per facet)

    // Add an attribute on facets. see https://github.com/BrunoLevy/geogram/wiki/Mesh#attributes
    // The type must be Attribute<index_t> to be used with geogram/mesh/mesh_halfedges.h . See MeshHalfedges::facet_region_
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);

    //fill the attribute
    while (ifs >> current_line) {
        try {
            current_label = std::stoul(current_line);
            if(current_label > 5) {
                fmt::println(Logger::err("I/O"),"In load_labeling(), each line must be a label in [0,5]\nBut found {}",current_label); Logger::err("I/O").flush();
                return false;
            }
        }
        catch (const std::exception& e) { // if stoul() failed
            fmt::println(Logger::err("I/O"),"In load_labeling(), each line must be an unsigned integer\nBut found '{}'\nException message : {}",current_line,e.what()); Logger::err("I/O").flush();
            return false;
        }
        if(current_line_number >= facets_number) {
            fmt::println(Logger::err("I/O"),"In load_labeling(), the number of labels is greater than the number of facets");
            fmt::println(Logger::err("I/O"),"Number of labels so far = {}",current_line_number+1);
            fmt::println(Logger::err("I/O"),"Number of facets = {}",facets_number); Logger::err("I/O").flush();
            return false;
        }
        label[current_line_number] = (index_t) current_label;
        current_line_number++;
    }

    ifs.close();

    //compare with expected size
    if (current_line_number != facets_number){
        fmt::println(Logger::err("I/O"),"In load_labeling(), the number of labels is lesser than the number of facets");
        fmt::println(Logger::err("I/O"),"Number of labels = {}",current_line_number+1);
        fmt::println(Logger::err("I/O"),"Number of facets = {}",facets_number); Logger::err("I/O").flush();
        return false;
    }

    return true;
}

bool save_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // fetch the labeling
    auto out = fmt::output_file(filename);
    FOR(f,mesh.facets.nb()) {
        out.print("{}",label[f]);
        if(f != mesh.facets.nb()-1)
            out.print("\n");
    }
    return true;
}

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

void flip_labeling(Mesh& mesh, const char* attribute_name) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    FOR(f,mesh.facets.nb()) { // for each facet
        label[f] = opposite_label(label[f]);
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

void propagate_label(const Mesh& mesh, const char* attribute_name, index_t new_label, const std::set<index_t>& facets_in, const std::set<index_t> facets_out, const std::vector<index_t>& facet2chart, index_t chart_index) {
    // change label to `new_label` for all facets in `facets_in` and their neighbors, step by step.
    // 2 limits for the propagation:
    //  - stay on `chart_index`
    //  - do not cross facets in `facets_out`
    std::vector<index_t> to_process(facets_in.begin(),facets_in.end()); // set -> vector data structure
    geo_assert(!to_process.empty());
    Attribute<index_t> label(mesh.facets.attributes(),attribute_name);

    index_t current_facet = index_t(-1),
            adjacent_facet = index_t(-1);
    while (!to_process.empty()) {
        current_facet = to_process.back();
        to_process.pop_back();
        // don't check if `current_facet` is on `chart_index`, too restrictive is some cases
        if(label[current_facet] == new_label) {
            // facet already processed since insertion
            continue;
        }
        label[current_facet] = new_label;
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

size_t remove_surrounded_charts(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);

    // Get charts surronded by only 1 label
    // Broader fix than just looking at charts having 1 boundary

    size_t nb_invalid_charts_processed = 0;
    for(index_t chart_index : slg.invalid_charts) {
        auto boundary_iterator = slg.charts[chart_index].boundaries.begin();
        index_t chart_at_other_side = slg.boundaries[*boundary_iterator].chart_at_other_side(chart_index);
        index_t surrounding_label = slg.charts[chart_at_other_side].label;
        if(slg.boundaries[*boundary_iterator].on_feature_edge) {
            goto skip_modification; // by propagating the `surrounding_label`, we will lost a feature edge
        }

        boundary_iterator++; // go to next boundary
        for(;boundary_iterator != slg.charts[chart_index].boundaries.end(); ++boundary_iterator) {
            if(slg.boundaries[*boundary_iterator].on_feature_edge) {
                goto skip_modification; // by propagating the `surrounding_label`, we will lost a feature edge
            }
            chart_at_other_side = slg.boundaries[*boundary_iterator].chart_at_other_side(chart_index);
            if(surrounding_label != slg.charts[chart_at_other_side].label) {
                goto skip_modification; // this chart has several labels around (goto because break in nested loops is a worse idea)
            }
        }

        // if we are here, the goto was not used, so the only label around is surrounding_label
        for(index_t facet_index : slg.charts[chart_index].facets) {
            label[facet_index] = surrounding_label;
        }
        nb_invalid_charts_processed++;


        skip_modification:
            ; // just end this loop interation
    }

    return nb_invalid_charts_processed;
}

bool fix_an_invalid_boundary(
    GEO::Mesh& mesh,
    const char* attribute_name,
    StaticLabelingGraph& slg,
    const std::vector<vec3>& facet_normals,
    const std::set<std::pair<index_t,index_t>>& feature_edges,
    const std::vector<std::vector<index_t>>& adj_facets
) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_he(mesh); // create an halfedges interface for this mesh
    mesh_he.set_use_facet_region(attribute_name);

    // For each invalid boundary,
    // Get facets at left and right
    // Find best side to put the new chart
    // Change label along this side
    //
    // TODO distribute the new label on the adjacent chart whole surface
    // so that the new chart is not that narrow
    // See MAMBO B29 for example
    // https://gitlab.com/franck.ledoux/mambo

    index_t new_label = index_t(-1);
    vec3 new_label_as_vector;
    MeshHalfedges::Halfedge current_halfedge;
    for(index_t boundary_index : slg.invalid_boundaries) { // for each invalid boundary

        // get ref to boundary object and retreive geometry and neighborhood info

        const Boundary& current_boundary = slg.boundaries[boundary_index];
        MeshHalfedges::Halfedge an_halfedge_of_the_invalid_boundary = current_boundary.halfedges[0];
        const Chart& left_chart = slg.charts[current_boundary.left_chart];
        const Chart& right_chart = slg.charts[current_boundary.right_chart];
        new_label = find_optimal_label(
            {},
            {},
            {left_chart.label,right_chart.label}, // new label must be orthogonal to the labels at left/right of the boundary
            current_boundary.average_normal // new label must be close to the boundary normal
        );
        new_label_as_vector = label2vector[new_label];
        std::set<index_t> left_facets_along_boundary;
        std::set<index_t> right_facets_along_boundary;
        index_t a_facet_on_new_chart = index_t(-1);
        current_boundary.get_adjacent_facets(mesh,left_facets_along_boundary,OnlyLeft,slg.facet2chart,1); // get triangles at left and at distance 0 or 1 from the boundary
        current_boundary.get_adjacent_facets(mesh,right_facets_along_boundary,OnlyRight,slg.facet2chart,1); // get triangles at right and at distance 0 or 1 from the boundary

        // Find out if the chart to insert on the boundary is better suited on the left side or the right side
        // Don't manage case where both sides are equally suited,
        // because the boundary could be on a feature edge, and we shouldn't de-capture it
        
        double average_dot_product_left = average_angle(facet_normals,left_facets_along_boundary,new_label_as_vector);
        double average_dot_product_right = average_angle(facet_normals,right_facets_along_boundary,new_label_as_vector);
        // comparison
        if (average_dot_product_left < average_dot_product_right) { // on average, the left side is a better place to put the new chart
            for(index_t f : left_facets_along_boundary) {
                label[f] = new_label;
            }
            a_facet_on_new_chart = *left_facets_along_boundary.begin();
        }
        else { // on average, the right side is a better place to put the new chart
            for(index_t f : right_facets_along_boundary) {
                label[f] = new_label;
            }
            a_facet_on_new_chart = *right_facets_along_boundary.begin();
        }

        // If the new chart is invalid
        // by having only one adjacent chart (the one at the other side of what was the invalid boundary)
        // extend it from the two corners until we found boundaries
        // and change the label of all facets below
        // See MAMBO B13 model for example
        // https://gitlab.com/franck.ledoux/mambo

        geo_assert(!adj_facets.empty()); // we need adjacency between vertices and facets for this part

        slg.fill_from(mesh,attribute_name,feature_edges);
        mesh_he.set_use_facet_region(attribute_name); // update facet regions

        geo_assert(slg.boundaries[slg.halfedge2boundary[an_halfedge_of_the_invalid_boundary].first].axis != -1);
        index_t axis_of_just_fixed_boundary = (index_t) slg.boundaries[slg.halfedge2boundary[an_halfedge_of_the_invalid_boundary].first].axis;

        const Chart& created_chart = slg.charts[slg.facet2chart[a_facet_on_new_chart]];

        if(created_chart.boundaries.size() > 2) {
            // nothing more to do
            return true;
        }

        // get the non-monotone boundary among the 2 boundaries around the created chart
        // it should have 2 turning-points
        // find the one toward positive `axis_of_just_fixed_boundary` and trace a path in the same direction
        // change label on one side
        // find the one toward negative `axis_of_just_fixed_boundary` and trace a path in the same direction
        // change label on one side
        index_t non_monotone_boundary = index_t(-1);
        TurningPoint tp0;
        TurningPoint tp1;
        if(slg.boundaries[*created_chart.boundaries.begin()].turning_points.size() == 2) {
            non_monotone_boundary = *created_chart.boundaries.begin();
            tp0 = slg.boundaries[non_monotone_boundary].turning_points[0];
            tp1 = slg.boundaries[non_monotone_boundary].turning_points[1];
        }
        else if(slg.boundaries[*created_chart.boundaries.rbegin()].turning_points.size() == 2) {
            non_monotone_boundary = *created_chart.boundaries.rbegin();
            tp0 = slg.boundaries[non_monotone_boundary].turning_points[0];
            tp1 = slg.boundaries[non_monotone_boundary].turning_points[1];
        }
        else {
            fmt::println(Logger::err("fix validity"),"In fix_as_much_invalid_boundaries_as_possible(), did not find the boundary with 2 turning-points in the contour of the created chart");
            continue;
        }
        const TurningPoint& turning_point_at_max_coordinate_on_axis = 
            mesh_vertex(mesh,tp0.vertex(slg.boundaries[non_monotone_boundary],mesh))[axis_of_just_fixed_boundary] >
            mesh_vertex(mesh,tp1.vertex(slg.boundaries[non_monotone_boundary],mesh))[axis_of_just_fixed_boundary] ?
            tp0 : tp1;
        const TurningPoint& turning_point_at_min_coordinate_on_axis = turning_point_at_max_coordinate_on_axis == tp0 ?
            tp1 : tp0;

        std::vector<MeshHalfedges::Halfedge> path;
        std::set<index_t> facets_at_left;
        std::set<index_t> facets_at_right;
        trace_path_on_chart(
            mesh_he,
            adj_facets,
            slg.facet2chart,
            slg.turning_point_vertices,
            turning_point_at_max_coordinate_on_axis.vertex(slg.boundaries[non_monotone_boundary],mesh),
            label2vector[axis_of_just_fixed_boundary*2], // positive direction
            facets_at_left,
            facets_at_right,
            path
        );

        if(dot(normalize(halfedge_midpoint_to_left_facet_tip_vector(mesh,path[0])),label2vector[created_chart.label]) >
            dot(normalize(halfedge_midpoint_to_right_facet_tip_vector(mesh,path[0])),label2vector[created_chart.label])
        ) {
            // facets_at_right -> wall
            // facets_at_left -> new chart
            propagate_label(
                mesh,
                attribute_name,
                created_chart.label,
                facets_at_left,
                facets_at_right,
                slg.facet2chart,
                slg.facet2chart[*facets_at_left.begin()]
            );
        }
        else {
            // facets_at_right -> new chart
            // facets_at_left -> wall
            propagate_label(
                mesh,
                attribute_name,
                created_chart.label,
                facets_at_right,
                facets_at_left,
                slg.facet2chart,
                slg.facet2chart[*facets_at_right.begin()]
            );
        }

        path.clear();
        facets_at_left.clear();
        facets_at_right.clear();
        trace_path_on_chart(
            mesh_he,
            adj_facets,
            slg.facet2chart,
            slg.turning_point_vertices,
            turning_point_at_min_coordinate_on_axis.vertex(slg.boundaries[non_monotone_boundary],mesh),
            label2vector[axis_of_just_fixed_boundary*2+1], // negative direction
            facets_at_left,
            facets_at_right,
            path
        );

        if(dot(normalize(halfedge_midpoint_to_left_facet_tip_vector(mesh,path[0])),label2vector[created_chart.label]) >
            dot(normalize(halfedge_midpoint_to_right_facet_tip_vector(mesh,path[0])),label2vector[created_chart.label])
        ) {
            // facets_at_right -> wall
            // facets_at_left -> new chart
            propagate_label(
                mesh,
                attribute_name,
                created_chart.label,
                facets_at_left,
                facets_at_right,
                slg.facet2chart,
                slg.facet2chart[*facets_at_left.begin()]
            );
        }
        else {
            // facets_at_right -> new chart
            // facets_at_left -> wall
            propagate_label(
                mesh,
                attribute_name,
                created_chart.label,
                facets_at_right,
                facets_at_left,
                slg.facet2chart,
                slg.facet2chart[*facets_at_right.begin()]
            );
        }
        return true;
    }

    return false;
}

size_t fix_as_much_invalid_corners_as_possible(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<index_t>& facet2chart, const std::vector<std::vector<index_t>>& adj_facets) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_he(mesh); // create an halfedges interface for this mesh
    mesh_he.set_use_facet_region(attribute_name);

    // Replace the labels around invalid corners
    // by the nearest one from the vertex normal

    unsigned int nb_invalid_corners_processed = 0;
    index_t new_label;
    vec3 vertex_normal = {0,0,0};
    std::vector<MeshHalfedges::Halfedge> outgoing_halfedges_on_feature_edge;
    for(index_t corner_index : slg.invalid_corners) { // for each invalid corner

        // check if there is a feature edge along the corner
        slg.corners[corner_index].get_outgoing_halfedges_on_feature_edge(mesh,feature_edges,outgoing_halfedges_on_feature_edge);
        // also compute the evenness of angles between adj boundaries
        double sd_boundary_angles = slg.corners[corner_index].sd_boundary_angles(mesh);
        if(outgoing_halfedges_on_feature_edge.empty() || sd_boundary_angles < 0.02) {

            // cone-like, or pyramid-like (MAMBO B20) invalid corner
            // compute the normal by adding normals of adjacent facets

            geo_assert(!adj_facets.empty());
            vertex_normal = {0,0,0};
            for(index_t f : adj_facets[slg.corners[corner_index].vertex]) {
                vertex_normal += normals[f];
            }
            vertex_normal /= (double) adj_facets[slg.corners[corner_index].vertex].size();

            // compute new label
            new_label = nearest_label(vertex_normal);

            // change the label of adjacent facets
            for(index_t f : adj_facets[slg.corners[corner_index].vertex]) {
                label[f] = new_label;
            }

            nb_invalid_corners_processed++;
        }
        else {
            // there is some kind of feature-edges pinching

            MeshHalfedges::Halfedge halfedge;
            if(vertex_has_lost_feature_edge_in_neighborhood(mesh_he,adj_facets,feature_edges,slg.corners[corner_index].vertex,halfedge)) {

                // There is an outgoing feature edge not captured (ie same label on both sides),
                // Follow the feature edge.
                // If we arrive at another invalid corner, trace a chart between the two invalid corners
                // Else do nothing
                // See MAMBO S9 for example
                // https://gitlab.com/franck.ledoux/mambo/

                index_t current_chart = slg.facet2chart[halfedge_facet_left(mesh,halfedge)];
                index_t label_of_current_chart = slg.charts[current_chart].label;
                index_t new_label = index_t(-1);
                geo_assert(label_of_current_chart == label[halfedge_facet_right(mesh,halfedge)]);
                std::set<index_t> facets_at_left;
                std::set<index_t> facets_at_right;
                index_t adjacent_facet = index_t(-1);
                halfedge = follow_feature_edge_on_chart(mesh_he,halfedge,feature_edges,slg.facet2chart,facets_at_left,facets_at_right);
                index_t corner_found = slg.vertex2corner[halfedge_vertex_index_to(mesh,halfedge)]; // can be -1
                if( (corner_found != index_t(-1)) && VECTOR_CONTAINS(slg.invalid_corners,corner_found) ) {
                    // we found a lost feature edge between 2 invalid corners
                    vec3 avg_normal_at_left = average_facets_normal(normals,facets_at_left);
                    vec3 avg_normal_at_right = average_facets_normal(normals,facets_at_right);
                    if(dot(avg_normal_at_left,label2vector[label_of_current_chart]) < dot(avg_normal_at_right,label2vector[label_of_current_chart])) {
                        // better to change the label of `facets_at_left`
                        new_label = find_optimal_label(
                            {},
                            {},
                            {label_of_current_chart}, // must be orthogonal to the current label
                            avg_normal_at_left // should be close to the `facets_at_left` normals
                        );
                        geo_assert(new_label != label_of_current_chart);
                        // change the label of the facets at the left & their neighbors
                        for(auto f : facets_at_left) {
                            label[f] = new_label;
                            FOR(le,3) { // for each local edge of facet f
                                adjacent_facet = mesh.facets.adjacent(f,le);
                                if(facets_at_right.contains(adjacent_facet)) {
                                    continue;
                                }
                                if(slg.facet2chart[adjacent_facet] == current_chart) {
                                    label[adjacent_facet] = new_label; // also change the label of the adjacent facet
                                }
                            }
                        }
                    }
                    else {
                        // better to change the label of `facets_at_right`
                        new_label = find_optimal_label(
                            {},
                            {},
                            {label_of_current_chart}, // must be orthogonal to the current label
                            avg_normal_at_right // should be close to the `facets_at_right` normals
                        );
                        geo_assert(new_label != label_of_current_chart);
                        // change the label of the facets at the left & their neighbors
                        for(auto f : facets_at_right) {
                            label[f] = new_label;
                            FOR(le,3) { // for each local edge of facet f
                                adjacent_facet = mesh.facets.adjacent(f,le);
                                if(facets_at_left.contains(adjacent_facet)) {
                                    continue;
                                }
                                if(slg.facet2chart[adjacent_facet] == current_chart) {
                                    label[adjacent_facet] = new_label; // also change the label of the adjacent facet
                                }
                            }
                        }
                    }

                    slg.fill_from(mesh,attribute_name,feature_edges);

                    // we need to stop fix_as_much_invalid_corners_as_possible() now
                    // because the loop iterating over `slg.invalid_corners` is no longer up to date
                    return nb_invalid_corners_processed;
                }
                // else : don't know how to fix, ignore
            }
            else {

                if(outgoing_halfedges_on_feature_edge.size() != 4) {
                    continue;
                }

                // A problematic corner on feature edge like on MAMBO S24 model https://gitlab.com/franck.ledoux/mambo/
                // Get the 4 outgoing boundaries, all of which being on feature edges
                // Keep only the 2 shortest
                // Store the chart they have in common
                // Change the label on their side where it's not their chart in common, with the label of the chart in common
                //
                // Only changes facets at distance of 0 or 1 from the 2 shortest boundaries.
                // Works for MAMBO S24 but it's not generalizable

                std::vector<std::pair<index_t,double>> adjacent_boundaries_and_their_length;
                for(const auto& halfedge : outgoing_halfedges_on_feature_edge) {
                    index_t b = slg.halfedge2boundary.at(halfedge).first;
                    adjacent_boundaries_and_their_length.push_back(std::make_pair(
                        b, // boundary index
                        slg.boundaries[b].length(mesh) // boundary length
                    ));
                }
                std::ranges::sort(
                    adjacent_boundaries_and_their_length,
                    [](const std::pair<index_t,double>& a, const std::pair<index_t,double>& b) { return a.second < b.second; } // sort by second item in the std::pair (the length)
                );
                // get the 2 shortest boundaries
                const Boundary& b0 = slg.boundaries[adjacent_boundaries_and_their_length[0].first];
                const Boundary& b1 = slg.boundaries[adjacent_boundaries_and_their_length[1].first];
                index_t chart_in_common = adjacent_chart_in_common(b0,b1);
                index_t label_of_chart_in_common = slg.charts[chart_in_common].label;

                std::set<index_t> facets_to_edit;
                b0.get_adjacent_facets(
                    mesh,
                    facets_to_edit,
                    chart_in_common == b0.left_chart ? OnlyRight : OnlyLeft,
                    facet2chart,
                    1 // include facet at a distance of 1 (touch the boundary from a vertex)
                );
                // `facets_to_edit` is not cleared inside get_adjacent_facets(), we can call it a second time on top of the existing set values
                b1.get_adjacent_facets(
                    mesh,
                    facets_to_edit,
                    chart_in_common == b1.left_chart ? OnlyRight : OnlyLeft,
                    facet2chart,
                    1 // include facet at a distance of 1 (touch the boundary from a vertex)
                );
                for(index_t f : facets_to_edit) {
                    label[f] = label_of_chart_in_common;
                }
                nb_invalid_corners_processed++;
            }
        }
    }

    return nb_invalid_corners_processed;
}

// from https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/labeling_ops.cpp removeChartMutation()
bool remove_invalid_charts(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg) {

    if(slg.invalid_charts.empty()) {
        fmt::println(Logger::warn("fix_labeling"),"remove_invalid_charts canceled because there are no invalid charts"); Logger::warn("fix_labeling").flush();
        return false;
    }

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute

    // Fill invalid charts using a Graph-Cut optimization,
    // preventing existing label from being re-applied

    // compactness = 1
    GraphCutLabeling gcl(mesh,normals);
    gcl.data_cost__set__locked_labels(label); // start by locking all the labels, so valid charts will not be modified
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(1);

    unsigned int nb_charts_to_remove = 0;
    for(index_t chart_index : slg.invalid_charts) { // for each invalid chart

        // prevent removal of chart surrounded by feature edges
        // TODO : increase their valence, by boundary duplication or boundary shift
        if(slg.charts[chart_index].is_surrounded_by_feature_edges(slg.boundaries)) {
            fmt::println(Logger::warn("fix_labeling"),"Cannot remove chart n°{} because it is surrounded by feature edges",chart_index); Logger::warn("fix_labeling").flush();
            continue;
        }

        for(index_t facet_index : slg.charts[chart_index].facets) { // for each facet inside this chart
            gcl.data_cost__change_to__fidelity_based(facet_index,1);
            // if facet next to a boundary, lower the cost of assigning the neighboring label
            FOR(le,3) { // for each local edge
                index_t adjacent_facet = mesh.facets.adjacent(facet_index,le);
                if(label[adjacent_facet] != label[facet_index]) {
                    gcl.data_cost__change_to__scaled(facet_index,label[adjacent_facet],0.5f); // halve the cost
                }
            }
            gcl.data_cost__change_to__forbidden_polycube_label(facet_index,label[facet_index]); // prevent the label from staying the same
        }
        nb_charts_to_remove++;
    }

    if(nb_charts_to_remove == 0) {
        // all invalid charts are surrounded by feature edges
        return true;
    }

    gcl.compute_solution(label);
    return false;
}

void remove_charts_around_invalid_boundaries(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg) {
    
    if(slg.invalid_boundaries.empty()) {
        fmt::println(Logger::out("fix_labeling"),"Warning : operation canceled because there are no invalid boundaries"); Logger::out("fix_labeling").flush();
        return;
    }

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute

    // In involved charts, preventing existing label from being re-applied
    // Lock labels in other charts

    std::set<index_t> charts_to_remove;
    for(index_t b : slg.invalid_boundaries) {
        charts_to_remove.insert(slg.boundaries[b].left_chart);
        charts_to_remove.insert(slg.boundaries[b].right_chart);
    }

    GraphCutLabeling gcl(mesh,normals);
    gcl.data_cost__set__locked_labels(label); // start by locking all the labels, so other charts will not be modified
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(1);

    for(index_t chart_index : charts_to_remove) { // for each chart to remove

        for(index_t facet_index : slg.charts[chart_index].facets) { // for each facet inside this chart
            gcl.data_cost__change_to__fidelity_based(facet_index,1);
            // if facet next to a boundary, lower the cost of assigning the neighboring label
            FOR(le,3) { // for each local edge
                index_t adjacent_facet = mesh.facets.adjacent(facet_index,le);
                if(label[adjacent_facet] != label[facet_index]) {
                    gcl.data_cost__change_to__scaled(facet_index,label[adjacent_facet],0.5f); // halve the cost
                }
            }
            gcl.data_cost__change_to__forbidden_polycube_label(facet_index,label[facet_index]); // prevent the label from staying the same
        }
    }

    gcl.compute_solution(label);

}

bool increase_chart_valence(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, StaticLabelingGraph& slg, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges) {

    FOR(invalid_chart_index,slg.invalid_charts.size()) {
        index_t chart_index = slg.invalid_charts[invalid_chart_index];
        const Chart& chart = slg.charts[chart_index];
        if(!chart.is_surrounded_by_feature_edges(slg.boundaries)) {
            continue; // check next invalid chart
        }
        Attribute<index_t> label(mesh.facets.attributes(), attribute_name);
        CustomMeshHalfedges mesh_he(mesh);
        mesh_he.set_use_facet_region(attribute_name);

        // transform the set of boundaries around `chart` into a vector

        // among all boundaries around `chart`, find
        //  - a corner between 2 boundaries that are associated to the same axis (kind of 2D non-orthogonal separation = 2D invalidity)
        // or
        //  - a boundary having a turning-point on a feature edge
        //    assume if there is a turning point, there is only one turning point in the chart contour

        index_t problematic_corner = index_t(-1);
        index_t problematic_non_monotone_boundary = index_t(-1);
        index_t problematic_vertex = index_t(-1);
        Boundary downward_boundary;
        Boundary upward_boundary;
        index_t current_axis = index_t(-1);
        index_t axis_to_insert = index_t(-1);

        std::vector<std::pair<index_t,bool>> counterclockwise_order;
        chart.counterclockwise_boundaries_order(mesh_he,slg.halfedge2boundary,slg.boundaries,counterclockwise_order);

        FOR(lb,counterclockwise_order.size()) { // for each local boundary index (relative to the chart contour)
            auto [b,b_is_same_direction] = counterclockwise_order[lb]; // b is a boundary index
            if( !slg.boundaries[b].turning_points.empty() && (slg.boundaries[b].axis != -1) ) {
                // naive labeling of MAMBO S36 has 2 turning-points around one of its invalid chart surrounded by feature edges
                // only one of them is relevant for this operator -> the closest to the boundary midpoint
                index_t ltp = 0; // local turning-point index
                if(slg.boundaries[b].turning_points.size() > 1) {
                    size_t max_dist_to_corners = 0; // in number of edges
                    FOR(current_ltp,slg.boundaries[b].turning_points.size()) { // for each turning-point of the boundary `b`
                        size_t current_dist = std::min(
                            (size_t) slg.boundaries[b].turning_points[current_ltp].outgoing_local_halfedge_index_,
                            slg.boundaries[b].halfedges.size() - ((size_t) slg.boundaries[b].turning_points[current_ltp].outgoing_local_halfedge_index_) + 1
                        );
                        if(current_dist > max_dist_to_corners) {
                            max_dist_to_corners = current_dist;
                            ltp = current_ltp;
                        }
                    }
                }
                current_axis = (index_t) slg.boundaries[b].axis;
                problematic_non_monotone_boundary = b;
                Boundary current_boundary_as_counterclockwise;
                if(b_is_same_direction) {
                    current_boundary_as_counterclockwise = slg.boundaries[b];
                }
                else {
                    slg.boundaries[b].get_flipped(mesh_he,current_boundary_as_counterclockwise);
                    ltp = ((index_t) slg.boundaries[b].turning_points.size()) - 1 - ltp;
                }
                current_boundary_as_counterclockwise.split_at_turning_point(mesh_he,downward_boundary,upward_boundary,ltp);
                axis_to_insert = nearest_axis_of_edges(mesh,{
                    slg.boundaries[b].halfedges[slg.boundaries[b].turning_points[ltp].outgoing_local_halfedge_index_],
                    slg.boundaries[b].halfedges[slg.boundaries[b].turning_points[ltp].outgoing_local_halfedge_index_-1]
                },{current_axis});
                problematic_vertex = current_boundary_as_counterclockwise.turning_points[ltp].vertex(current_boundary_as_counterclockwise,mesh);
                break;
            }
            geo_assert(slg.boundaries[b].start_corner != index_t(-1));
            geo_assert(slg.boundaries[b].end_corner != index_t(-1));
            auto [next_b,next_b_is_same_direction] = counterclockwise_order[(lb+1) % counterclockwise_order.size()];
            if(slg.boundaries[b].axis == slg.boundaries[next_b].axis) {
                geo_assert(slg.boundaries[b].axis != -1);
                current_axis = (index_t) slg.boundaries[b].axis;
                if(b_is_same_direction) {
                    slg.boundaries[b].get_flipped(mesh_he,downward_boundary);
                    problematic_corner = slg.boundaries[b].end_corner;
                }
                else {
                    downward_boundary = slg.boundaries[b];
                    problematic_corner = slg.boundaries[b].start_corner;
                }
                // if this is a pyramid-like invalid corner (see MAMBO B20),
                // don't apply increase_chart_valence() but fix_as_much_invalid_corners_as_possible()
                double sd_boundary_angles = slg.corners[problematic_corner].sd_boundary_angles(mesh);
                if(
                    (slg.corners[problematic_corner].valence() == 4) && 
                    slg.corners[problematic_corner].all_adjacent_boundary_edges_are_on_feature_edges(mesh,feature_edges) &&
                    sd_boundary_angles < 0.02 // required to distinguish B20 from S36, S33
                ) {
                    problematic_corner = index_t(-1);
                    problematic_vertex = index_t(-1);
                    continue;
                }
                if(next_b_is_same_direction) {
                    upward_boundary = slg.boundaries[next_b];
                }
                else {
                    slg.boundaries[next_b].get_flipped(mesh_he,upward_boundary);
                }
                problematic_vertex = slg.corners[problematic_corner].vertex;
                axis_to_insert = nearest_axis_of_edges(mesh,{
                    downward_boundary.halfedges[0],
                    upward_boundary.halfedges[0]
                },{current_axis});
                break;
            }
        }

        if(problematic_vertex == index_t(-1)) {
            continue; // propagate 'continue'
        }
        geo_assert( (problematic_corner != index_t(-1)) || (problematic_non_monotone_boundary != index_t(-1)) );
        
        #ifndef NDEBUG
            dump_vertex("problematic_vertex",mesh,problematic_vertex);
            dump_boundary_with_halfedges_indices("downward_boundary",mesh,downward_boundary);
            dump_boundary_with_halfedges_indices("upward_boundary",mesh,upward_boundary);
        #endif

        std::vector<double> downward_boundary_cumulative_cost_for_current_axis;
        std::vector<double> downward_boundary_cumulative_cost_for_axis_to_insert;
        std::vector<double> upward_boundary_cumulative_cost_for_current_axis;
        std::vector<double> upward_boundary_cumulative_cost_for_axis_to_insert;
        downward_boundary.per_edges_cumulative_axis_assignement_cost(mesh,current_axis,downward_boundary_cumulative_cost_for_current_axis,false);
        downward_boundary.per_edges_cumulative_axis_assignement_cost(mesh,axis_to_insert,downward_boundary_cumulative_cost_for_axis_to_insert,true);
        upward_boundary.per_edges_cumulative_axis_assignement_cost(mesh,current_axis,upward_boundary_cumulative_cost_for_current_axis,false);
        upward_boundary.per_edges_cumulative_axis_assignement_cost(mesh,axis_to_insert,upward_boundary_cumulative_cost_for_axis_to_insert,true);

        // find equilibrium along `downward_boundary`

        index_t downward_boundary_equilibrium_vertex_index = 0;
        while(downward_boundary_cumulative_cost_for_axis_to_insert[downward_boundary_equilibrium_vertex_index] <= downward_boundary_cumulative_cost_for_current_axis[downward_boundary_equilibrium_vertex_index]) {
            downward_boundary_equilibrium_vertex_index++;
            if(downward_boundary_equilibrium_vertex_index == downward_boundary.halfedges.size()) {
                break;
            }
        }
        index_t downward_boundary_equilibrium_vertex = (downward_boundary_equilibrium_vertex_index == downward_boundary.halfedges.size()) ?
            halfedge_vertex_index_to(mesh,downward_boundary.halfedges[downward_boundary_equilibrium_vertex_index-1]) :
            halfedge_vertex_index_from(mesh,downward_boundary.halfedges[downward_boundary_equilibrium_vertex_index]);
        #ifndef NDEBUG
            dump_vertex("downward_boundary_equilibrium",mesh,downward_boundary_equilibrium_vertex);
        #endif

        // find equilibrium along `upward_boundary`

        index_t upward_boundary_equilibrium_vertex_index = 0;
        while(upward_boundary_cumulative_cost_for_axis_to_insert[upward_boundary_equilibrium_vertex_index] <= upward_boundary_cumulative_cost_for_current_axis[upward_boundary_equilibrium_vertex_index]) {
            upward_boundary_equilibrium_vertex_index++;
            if(upward_boundary_equilibrium_vertex_index == upward_boundary.halfedges.size()) {
                break;
            }
        }
        index_t upward_boundary_equilibrium_vertex = (upward_boundary_equilibrium_vertex_index == upward_boundary.halfedges.size()) ?
            halfedge_vertex_index_to(mesh,upward_boundary.halfedges[upward_boundary_equilibrium_vertex_index-1]) :
            halfedge_vertex_index_from(mesh,upward_boundary.halfedges[upward_boundary_equilibrium_vertex_index]);
        #ifndef NDEBUG
            dump_vertex("upward_boundary_equilibrium",mesh,upward_boundary_equilibrium_vertex);
        #endif

        geo_assert(downward_boundary_equilibrium_vertex != upward_boundary_equilibrium_vertex);

        // If neither the downward equilibrium point nor the upward one are on the `problematic_vertex`,
        // move the closest onto the `problematic_vertex`,
        // so that we always have a boundary starting from an equilibrium point and one starting from the `problematic_vertex`
        // Avoids strange behavior on MAMBO B41
        if(
            (downward_boundary_equilibrium_vertex != problematic_vertex) && 
            (upward_boundary_equilibrium_vertex != problematic_vertex)
        ) {
            if(
                upward_boundary_equilibrium_vertex_index < downward_boundary_equilibrium_vertex_index
            ) {
                // `upward_boundary_equilibrium_vertex` is the closest
                upward_boundary_equilibrium_vertex = problematic_vertex;
            }
            else {
                downward_boundary_equilibrium_vertex = problematic_vertex;
            }
        }

        // Compute the direction to follow
        // Between
        //   label2vector[chart.label] <=> same direction as the label of the current `chart` (like with MAMBO B49)
        // and
        //   label2vector[opposite_label(chart.label)] <=> opposite direction (like with MAMBO B21)
        // choose the one for which the most aligned outgoing halfedge is the closer
        // Check for both
        //   `downward_boundary_equilibrium_vertex`
        // and
        //   `upward_boundary_equilibrium_vertex`

        MeshHalfedges::Halfedge outgoing_halfedge;
        vec3 direction;
        
        // neighborhood of `downward_boundary_equilibrium_vertex`
        outgoing_halfedge = get_an_outgoing_halfedge_of_vertex(mesh,adj_facets,downward_boundary_equilibrium_vertex);
        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(outgoing_halfedge,mesh_he,label2vector[chart.label]);
        double max_angle_same_direction = angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[chart.label]);

        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(outgoing_halfedge,mesh_he,label2vector[opposite_label(chart.label)]);
        double max_angle_opposite_direction = angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[opposite_label(chart.label)]);

        // neighborhood of `upward_boundary_equilibrium_vertex`
        outgoing_halfedge = get_an_outgoing_halfedge_of_vertex(mesh,adj_facets,upward_boundary_equilibrium_vertex);
        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(outgoing_halfedge,mesh_he,label2vector[chart.label]);
        max_angle_same_direction = std::max(max_angle_same_direction,angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[chart.label]));

        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(outgoing_halfedge,mesh_he,label2vector[opposite_label(chart.label)]);
        max_angle_opposite_direction = std::max(max_angle_opposite_direction,angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[opposite_label(chart.label)]));

        if(max_angle_same_direction < max_angle_opposite_direction) {
            // better to go in the same direction as the label of the current chart
            direction = label2vector[chart.label];
        }
        else {
            // better to go in the opposite direction
            direction = label2vector[opposite_label(chart.label)];
        }

        // Trace a boundary from each equilibrium point
        // If one of the equilibrium point is on a corner, re-use the boundary going out of the chart

        std::vector<MeshHalfedges::Halfedge> path;
        std::set<index_t> facets_of_new_chart;
        std::set<index_t> walls;

        if( (downward_boundary_equilibrium_vertex == problematic_vertex) && (problematic_corner != index_t(-1)) ) {
            // This equilibrium point is on the problematic corner (between 2 boundaries assigned to the same axis)
            // So there is already an outgoing boundary, we just need to fetch facets at its left and right
            MeshHalfedges::Halfedge halfedge = slg.corners[problematic_corner].get_most_aligned_boundary_halfedge(mesh,direction);
            geo_assert(halfedge_vertex_index_from(mesh,halfedge) == downward_boundary_equilibrium_vertex);
            auto [boundary_to_reuse, same_direction] = slg.halfedge2boundary[halfedge];
            geo_assert(boundary_to_reuse != index_t(-1));
            if(same_direction) {
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyLeft,slg.facet2chart);
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyRight,slg.facet2chart);
            }
            else {
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyRight,slg.facet2chart);
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyLeft,slg.facet2chart);
            }
        }
        else {
            // trace path
            trace_path_on_chart(mesh_he,adj_facets,slg.facet2chart,slg.turning_point_vertices,downward_boundary_equilibrium_vertex,direction*10,facets_of_new_chart,walls,path);
        }

        if( (upward_boundary_equilibrium_vertex == problematic_vertex) && (problematic_corner != index_t(-1)) ) {
            // This equilibrium point is on the problematic corner (between 2 boundaries assigned to the same axis)
            // So there is already an outgoing boundary, we just need to fetch facets at its left and right
            MeshHalfedges::Halfedge halfedge = slg.corners[problematic_corner].get_most_aligned_boundary_halfedge(mesh,direction);
            geo_assert(halfedge_vertex_index_from(mesh,halfedge) == upward_boundary_equilibrium_vertex);
            auto [boundary_to_reuse, same_direction] = slg.halfedge2boundary[halfedge];
            geo_assert(boundary_to_reuse != index_t(-1));
            if(same_direction) {
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyLeft,slg.facet2chart);
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyRight,slg.facet2chart);
            }
            else {
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyRight,slg.facet2chart);
                slg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyLeft,slg.facet2chart);
            }
        }
        else {
            // trace path
            trace_path_on_chart(mesh_he,adj_facets,slg.facet2chart,slg.turning_point_vertices,upward_boundary_equilibrium_vertex,direction*10,walls,facets_of_new_chart,path);
        }

        index_t chart_on_wich_the_new_chart_will_be = slg.facet2chart[*facets_of_new_chart.begin()];

        #ifndef NDEBUG
            dump_facets("facets_of_new_chart",mesh,facets_of_new_chart);
            dump_facets("walls",mesh,walls);
        #endif

        index_t label_to_insert = find_optimal_label(
            { // 2 forbidden axes:
                chart.label/2, // the axis of the chart we're going to increase the valence
                slg.charts[chart_on_wich_the_new_chart_will_be].label/2 // the axis of the chart on which we traced boundaries
            },
            {},
            {}, // no need to specify orthogonality constraints, becase we already forbid 2 axes over 3
            average_facets_normal(normals,facets_of_new_chart) // among the 2 remaining labels, choose the one the closest to the avg normal of `facets_of_new_chart`
        );

        propagate_label(mesh,attribute_name,label_to_insert,facets_of_new_chart,walls,slg.facet2chart,chart_on_wich_the_new_chart_will_be);

        return true; // an invalid chart has been processed
    }
    return false;
}

bool auto_fix_validity(Mesh& mesh, std::vector<vec3>& normals, const char* attribute_name, StaticLabelingGraph& slg, unsigned int max_nb_loop, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<vec3>& facet_normals, const std::vector<std::vector<index_t>>& adj_facets) {
    unsigned int nb_loops = 0;
    size_t nb_processed = 0;
    std::set<std::array<std::size_t,7>> set_of_labeling_features_combinations_encountered;
    while(!slg.is_valid() && nb_loops <= max_nb_loop) { // until valid labeling OR too much steps
        nb_loops++;

        // as much as possible, remove isolated (surrounded) charts
        do {
            nb_processed = remove_surrounded_charts(mesh,attribute_name,slg);
            // update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
            slg.fill_from(mesh,attribute_name,feature_edges);
        } while(nb_processed != 0);

        if(slg.is_valid())
            return true;

        do {
            nb_processed = (size_t) increase_chart_valence(mesh,normals,attribute_name,slg,adj_facets,feature_edges);
            slg.fill_from(mesh,attribute_name,feature_edges);
        } while(nb_processed != 0);

        if(slg.is_valid())
            return true;

        do {
            nb_processed = (size_t) fix_an_invalid_boundary(mesh,attribute_name,slg,facet_normals,feature_edges,adj_facets);
            slg.fill_from(mesh,attribute_name,feature_edges);
        }
        while (nb_processed != 0);

        if(slg.is_valid())
            return true;
        
        do {
            nb_processed = fix_as_much_invalid_corners_as_possible(mesh,normals,attribute_name,slg,feature_edges,slg.facet2chart,adj_facets);
            slg.fill_from(mesh,attribute_name,feature_edges);
        }
        while (nb_processed != 0);

        if(slg.is_valid())
            return true;

        set_of_labeling_features_combinations_encountered.clear();
        set_of_labeling_features_combinations_encountered.insert({
            slg.nb_charts(),
            slg.nb_boundaries(),
            slg.nb_corners(),
            slg.nb_invalid_charts(),
            slg.nb_invalid_boundaries(),
            slg.nb_invalid_corners(),
            slg.nb_turning_points()
        });

        while(1) {
            remove_invalid_charts(mesh,normals,attribute_name,slg);
            slg.fill_from(mesh,attribute_name,feature_edges);

            if(slg.is_valid())
                return true;

            std::array<std::size_t,7> features_combination = {
                slg.nb_charts(),
                slg.nb_boundaries(),
                slg.nb_corners(),
                slg.nb_invalid_charts(),
                slg.nb_invalid_boundaries(),
                slg.nb_invalid_corners(),
                slg.nb_turning_points()
            };

            if(VECTOR_CONTAINS(set_of_labeling_features_combinations_encountered,features_combination)) { // we can use VECTOR_CONTAINS() on sets because they also have find(), cbegin() and cend()
                // we backtracked
                // There is probably small charts that we can remove to help the fixing routine
                remove_charts_around_invalid_boundaries(mesh,normals,attribute_name,slg);
                slg.fill_from(mesh,attribute_name,feature_edges);
                break; // go back to the beginning of the loop, with other fix operators
            }
            else {
                set_of_labeling_features_combinations_encountered.insert(features_combination); // store the current combination of number of features
            }
        }
        
    }

    if(!slg.is_valid()) {
        fmt::println(Logger::out("fix_labeling"),"auto fix validity stopped (max nb loops reached), no valid labeling found"); Logger::out("fix_labeling").flush();
        return false;
    }
    
    return true;
}

size_t move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges) {

    size_t nb_labels_changed = 0;
    size_t nb_turning_points_moved = 0;

    if(slg.non_monotone_boundaries.empty()) {
        fmt::println(Logger::out("fix_labeling"),"Warning : operation canceled because all boundaries are monotone"); Logger::out("fix_labeling").flush();
        return nb_labels_changed;
    }

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_halfedges(mesh); // create an halfedges interface for this mesh
    mesh_halfedges.set_use_facet_region(attribute_name);
    MeshHalfedges::Halfedge initial_halfedge, current_halfedge;
    index_t new_label = index_t(-1);

    // Get turning points, change the label if this doesnt break charts connectivity

    for(auto b : slg.non_monotone_boundaries) {

        if(slg.boundaries[b].on_feature_edge) {
            continue; // don't move this boundary away from the feature edge
        }

        for(auto tp : slg.boundaries[b].turning_points) {

            initial_halfedge = slg.boundaries[b].halfedges[tp.outgoing_local_halfedge_index_];
            geo_assert(mesh_halfedges.halfedge_is_border(initial_halfedge));
            geo_assert(MAP_CONTAINS(slg.halfedge2boundary,initial_halfedge));

            // if one of the halfedges, the one before or the one after,
            // is on a feature edge, don't move this boundary

            if(halfedge_is_on_feature_edge(mesh,initial_halfedge,feature_edges)) {
                break; // don't move this boundary away from the feature edge
            }
            if(halfedge_is_on_feature_edge(mesh,
                slg.boundaries[b].halfedges[(tp.outgoing_local_halfedge_index_+1) % slg.boundaries[b].halfedges.size()],
                feature_edges
            )) {
                break; // don't move this boundary away from the feature edge
            }

            // get the vertex index
            index_t current_vertex = mesh.facet_corners.vertex(initial_halfedge.corner);
            geo_assert(slg.turning_point_vertices.contains(current_vertex));
            geo_assert(slg.vertex2corner[current_vertex] == index_t(-1)); // a turning point should not be a corner
            // test if the valence of current_vertex is 2
            VertexRingWithBoundaries vr;
            vr.explore(initial_halfedge,mesh_halfedges);
            geo_assert(vr.valence() == 2); // should have only 2 charts

            // new label according to towards which chart the turning point is
            new_label = slg.charts[tp.is_towards_left() ? slg.boundaries[b].left_chart : slg.boundaries[b].right_chart].label;

            current_halfedge = initial_halfedge;
            // go around the vertex and assign new_label to adjacent facets
            do {
                mesh_halfedges.move_counterclockwise_around_vertex(current_halfedge,true);
                if(label[current_halfedge.facet] != new_label) {
                    label[current_halfedge.facet] = new_label;
                    nb_labels_changed++;
                }
            } while (current_halfedge != initial_halfedge);

            nb_turning_points_moved++;
        }
    }
    return nb_turning_points_moved;
}

void straighten_boundary_with_GCO(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg, index_t boundary_index) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    const Boundary& current_boundary = slg.boundaries[boundary_index];
    index_t left_chart_index = current_boundary.left_chart;
    index_t right_chart_index = current_boundary.right_chart;
    const Chart& left_chart = slg.charts[left_chart_index];
    const Chart& right_chart = slg.charts[right_chart_index];
    CustomMeshHalfedges mesh_he(mesh);

    #ifndef NDEBUG
        fmt::println("Working on boundary {} -> charts {} and {}",boundary_index,left_chart_index,right_chart_index);
    #endif

    if(current_boundary.halfedges.size()<=4) {
        // small boundary -> skip
        #ifndef NDEBUG
            fmt::println("Skipped (4 or less edges)");
        #endif
        return;
    }

    #ifndef NDEBUG
        dump_all_boundaries("boundaries",mesh,slg.boundaries);

        std::map<index_t,bool> facets_of_the_2_charts;
        for(auto f : left_chart.facets) {
            facets_of_the_2_charts[f] = 0; // insert this facet, associated to 0=left
        }
        for(auto f : right_chart.facets) {
            facets_of_the_2_charts[f] = 1; // insert this facet, associated to 1=right
        }
        dump_facets("2_charts","on_wich_side",mesh,facets_of_the_2_charts);
    #endif

    std::map<index_t,unsigned int> distance_to_boundary; // map a facet index to a distance. only defined for facet inside the 2 charts
    for(auto f : left_chart.facets) {
        distance_to_boundary[f] = (unsigned int) -1;
    }
    for(auto f : right_chart.facets) {
        distance_to_boundary[f] = (unsigned int) -1;
    }
    // Go through the boundary and set distance of adjacent triangle to 0
    auto boundary_halfedge = current_boundary.halfedges.cbegin(); // get an iterator pointing at the first halfedge
    boundary_halfedge++; // go to the second halfedge (the vertex at the base of the first one has facets of other charts next to it)
    CustomMeshHalfedges::Halfedge current_halfedge;
    for(; boundary_halfedge != current_boundary.halfedges.cend(); ++boundary_halfedge) { // walk along the boundary, from (boundary) halfedge to (boundary) halfedge
        current_halfedge = (*boundary_halfedge); // copy the boundary halfedge into a mutable variable
        do { // for each facet around the vertex at the base of the current boundary halfedge
            distance_to_boundary[current_halfedge.facet] = 0; // set distance to 0 (the current facet touch the boundary by an edge or a vertex)
            mesh_he.move_counterclockwise_around_vertex(current_halfedge,true);
        } while (current_halfedge != *boundary_halfedge); // go around the vertex, until we are back on the initial boundary edge
    }
    per_facet_distance(mesh,distance_to_boundary);

    geo_assert(distance_to_boundary.size() == left_chart.facets.size()+right_chart.facets.size());

    #ifndef NDEBUG
        dump_facets("per_facet_dist","dist",mesh,distance_to_boundary);
        // Export the contour (facets in the perimeter of the union of the 2 charts)
        std::set<index_t> contour;
        for(const auto& kv : distance_to_boundary) { // for each facet in the 2 charts
            if ( !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,0)) || 
                 !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,1)) || 
                 !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,2)) ) { // if one of the neighbors is not inside the 2 charts
                contour.insert(kv.first);
            }
        }
        dump_facets("contour",mesh,contour);
    #endif

    auto gcl = GraphCutLabeling(mesh,normals, (index_t) distance_to_boundary.size(),{0,1,2,3,4,5}); // graph-cut only on the two charts
    for(const auto& kv : distance_to_boundary) {
        gcl.add_facet(kv.first);
    }
    gcl.data_cost__set__fidelity_based(1);
    // tweak data cost based on fidelity
    // the further the triangle is from the boundary, the lower is the cost of assigning it to its current label
    // also : if the facet is on the contour of the 2 charts, lock the label (prevent modification). Will lock the 2 corners
    for(const auto& kv : distance_to_boundary) {
        if ( !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,0)) || 
             !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,1)) || 
             !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,2)) ) {
            // this facet is on the contour of the 2 charts -> lock its label
            gcl.data_cost__change_to__locked_polycube_label(kv.first,label[kv.first]);
        }
        else {
            // penalize modification proportionally to the distance to the boundary (only close facets should be modifiable)
            // kv.second is the distance
            // -> multiply the cost of re-assigning by (0.8)^distance
            gcl.data_cost__change_to__scaled(kv.first,label[kv.first],(float) std::pow(0.8,kv.second));
        }
    }
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(1);
    if(current_boundary.axis != -1) {
        // penalize boundary tracing through halfedges poorly aligned with the boundary axis (X, Y or Z)
        // parse all undirected edge inside the 2 charts surface and tweak the cost
        #ifndef NDEBUG
            std::map<std::pair<index_t,index_t>,double> per_edge_dot_product;
        #endif
        index_t vertex1 = index_t(-1);
        index_t vertex2 = index_t(-1);
        index_t adjacent_facet = index_t(-1);
        vec3 edge_vector;
        double dot_product = 0.0; // dot product of the edge vector & boundary axis
        std::set<std::set<index_t>> already_processed; // use a set and not a pair so that the facets are sorted -> unordered pairs
        for(const auto& kv : distance_to_boundary) {
            FOR(le,3) { // for each local edge of the current facet
                adjacent_facet = mesh.facets.adjacent(kv.first,le);
                if(already_processed.contains({kv.first,adjacent_facet})) {
                    continue; // we already processed this edge (indices flipped)
                }
                if(distance_to_boundary.contains(adjacent_facet)) { // the edge at le is between 2 facets of the 2 charts union
                    // the vertices at the end point of the local edge le are le and (le+1)%3
                    // https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
                    vertex1 = mesh.facet_corners.vertex(mesh.facets.corner(kv.first,le));
                    vertex2 = mesh.facet_corners.vertex(mesh.facets.corner(kv.first,(le+1)%3));
                    edge_vector = normalize(mesh.vertices.point(vertex2)-mesh.vertices.point(vertex1));
                    dot_product = std::abs(dot(edge_vector,label2vector[current_boundary.axis*2])); // =1 -> //, =0 -> ⟂
                    gcl.neighbors__change_to__shifted(kv.first,adjacent_facet,(float) (1-dot_product)*500); // add a penalty the more ⟂ the edge is (// -> 0, ⟂ -> 300)
                    already_processed.insert({kv.first, adjacent_facet});
                    #ifndef NDEBUG
                        per_edge_dot_product[std::make_pair(vertex1,vertex2)] = dot_product;
                    #endif
                }
            }
        }
        #ifndef NDEBUG
            dump_edges("dot_products","dot_product",mesh,per_edge_dot_product);
        #endif
    }
    gcl.compute_solution(label);
}

bool straighten_boundary(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, index_t boundary_index, const std::vector<std::vector<index_t>>& adj_facets) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    const Boundary& current_boundary = slg.boundaries[boundary_index];

    #ifndef NDEBUG
        dump_boundary_with_halfedges_indices("boundary_to_straighten",mesh,current_boundary);
    #endif

    if(current_boundary.halfedges.size() <= 3) {
        return true; // too small, no need to straighten
    }

    // trace a new path for `current_boundary` between its two corners
    // by keeping its first and last halfedges, and, starting from the extremity of the first halfedge (= 2nd vertex),
    // move edge by edge toward the last halfedge, choosing the best-aligned edges

    index_t left_chart_index = current_boundary.left_chart;
    index_t right_chart_index = current_boundary.right_chart;
    const Chart& left_chart = slg.charts[left_chart_index];
    const Chart& right_chart = slg.charts[right_chart_index];
    CustomMeshHalfedges mesh_he(mesh);

    index_t target_vertex = halfedge_vertex_index_from(mesh,*current_boundary.halfedges.rbegin()); // origin of laft halfedge
    vec3 target_point = mesh_vertex(mesh,target_vertex);

    MeshHalfedges::Halfedge previous_halfedge;
    MeshHalfedges::Halfedge current_halfedge = current_boundary.halfedges[0]; // first halfedge
    
    // add edges one by one, aiming at `target_point`
    std::set<index_t> facets_at_left, facets_at_right;
    facets_at_left.insert(halfedge_facet_left(mesh,current_boundary.halfedges[0]));
    facets_at_left.insert(halfedge_facet_left(mesh,*current_boundary.halfedges.rbegin()));
    facets_at_right.insert(halfedge_facet_right(mesh,current_boundary.halfedges[0]));
    facets_at_right.insert(halfedge_facet_right(mesh,*current_boundary.halfedges.rbegin()));
    while (halfedge_vertex_index_to(mesh,current_halfedge) != target_vertex) {
        previous_halfedge = current_halfedge;
        mesh_he.move_to_opposite(previous_halfedge); // flip `previous_halfedge` so that its origin vertex is the origin vertex of the next halfedge
        current_halfedge = get_most_aligned_halfedge_around_vertex(previous_halfedge,mesh_he,target_point - halfedge_vertex_from(mesh,previous_halfedge));
        if(current_halfedge == previous_halfedge) {
            // we are backtracking
            fmt::println(Logger::warn("monotonicity"),"Cannot straighten boundary {} (backtracking)",boundary_index); Logger::warn("monotonicity").flush();
            return false;
        }
        if(!slg.vertex_is_only_surrounded_by(halfedge_vertex_index_to(mesh,current_halfedge),{left_chart_index,right_chart_index},adj_facets)) {
            // the `current_halfedge` is leading us away from the two charts on which the boundary must stay, we are going to cross another boundary
            // -> we have to straighten the other boundary first

            // for now, do not straighten the `current_boundary` 
            fmt::println(Logger::warn("monotonicity"),"Cannot straighten boundary {} (new path encountered another boundary)",boundary_index); Logger::warn("monotonicity").flush();
            return false;
        }
        facets_at_left.insert(halfedge_facet_left(mesh,current_halfedge));
        facets_at_right.insert(halfedge_facet_right(mesh,current_halfedge));
    }

    std::vector<index_t> left_facets_to_process(facets_at_left.begin(),facets_at_left.end());
    std::vector<index_t> right_facets_to_process(facets_at_right.begin(),facets_at_right.end());

    index_t current_facet = index_t(-1),
            adjacent_facet = index_t(-1);

    // process facets at left of new boundary path

    while (!left_facets_to_process.empty()) {
        current_facet = left_facets_to_process.back();
        left_facets_to_process.pop_back();
        if(label[current_facet] == left_chart.label) {
            // facet already processed since insertion
            continue;
        }
        label[current_facet] = left_chart.label;
        FOR(le,3) { // process adjacent triangles
            adjacent_facet = mesh.facets.adjacent(current_facet,le);
            if(label[adjacent_facet] == left_chart.label) {
                continue;
            }
            if(facets_at_right.contains(adjacent_facet)) {
                continue;
            }
            if(!left_chart.facets.contains(adjacent_facet) && !right_chart.facets.contains(adjacent_facet)) {
                continue;
            }
            left_facets_to_process.push_back(adjacent_facet);
        }
    }

    // process facets at right of new boundary path

    while (!right_facets_to_process.empty()) {
        current_facet = right_facets_to_process.back();
        right_facets_to_process.pop_back();
        if(label[current_facet] == right_chart.label) {
            // facet already processed since insertion
            continue;
        }
        label[current_facet] = right_chart.label;
        FOR(le,3) { // process adjacent triangles
            adjacent_facet = mesh.facets.adjacent(current_facet,le);
            if(label[adjacent_facet] == right_chart.label) {
                continue;
            }
            if(facets_at_left.contains(adjacent_facet)) {
                continue;
            }
            if(!left_chart.facets.contains(adjacent_facet) && !right_chart.facets.contains(adjacent_facet)) {
                continue;
            }
            right_facets_to_process.push_back(adjacent_facet);
        }
    }

    return true;
}

void straighten_boundaries(GEO::Mesh& mesh, const char* attribute_name, StaticLabelingGraph& slg, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    if(slg.boundaries.empty()) {
        fmt::println(Logger::out("monotonicity"),"No boundaries, operation canceled"); Logger::out("monotonicity").flush();
        return;
    }

    geo_assert(!adj_facets.empty())

    // Because the boundaries may not have the same index before and after slg.fill_from()
    // and because we need to call slg.fill_from() after that one boundary is processed to update charts,
    // we first gather the first boundary edge of all boundaries, so at each step we can get the current index of the
    // associated boundary with slg.halfedge2boundary
    // The fist boundary edges should not move after a call of straighten_boundary()...
    std::deque<MeshHalfedges::Halfedge> boundary_edges_to_process;
    boundary_edges_to_process.resize(slg.nb_boundaries()); // preallocation
    FOR(b,slg.boundaries.size()) {
        if (slg.boundaries[b].on_feature_edge) {
            fmt::println(Logger::out("monotonicity"),"Boundary {} skipped for straighten_boundary() because it is on a feature edge",b); Logger::out("monotonicity").flush();
            continue; // do not straighten boundaries surrounded by feature edges
        }
        geo_assert(!slg.boundaries[b].halfedges.empty());
        boundary_edges_to_process.push_back(slg.boundaries[b].halfedges[0]);
    }
    MeshHalfedges::Halfedge current_boundary_edge;
    index_t boundary_index = index_t(-1);
    bool boundary_in_same_direction = false;
    unsigned int count_iterations = 0;
    while (!boundary_edges_to_process.empty()) {
        current_boundary_edge = boundary_edges_to_process.back();
        boundary_edges_to_process.pop_back();
        std::tie(boundary_index,boundary_in_same_direction) = slg.halfedge2boundary[current_boundary_edge];
        geo_assert(boundary_index != index_t(-1));
        if(straighten_boundary(mesh,attribute_name,slg,boundary_index,adj_facets)) 
        {
            slg.fill_from(mesh,attribute_name,feature_edges);
        }
        else {
            boundary_edges_to_process.push_front(current_boundary_edge); // re-process this boundary edge later
        }
        count_iterations++;
        if(count_iterations > 100) {
            // if straighten_boundary() failed because backtracking and not because we encountered another boundary,
            // we could end up in an infinite loop where we process again and again a boundary for which we cannot reach the end corner...
            fmt::println(Logger::err("monotonicity"),"straighten_boundaries() stopped, reached max nb iter"); Logger::err("monotonicity").flush();
            break;
        }
    }
}

void move_corners(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_he(mesh);
    mesh_he.set_use_facet_region(attribute_name);
    std::set<index_t> adjacent_charts_of_corner; // indices of adjacent charts
    std::set<index_t> adjacent_charts_of_new_vertex;
    vec3 average_coordinates;
    index_t nearest_vertex_of_average_coordinates = index_t(-1);
    MeshHalfedges::Halfedge displacement_halfedge;
    unsigned int nb_corner_moved = 0;

    // we can process all corners at once, because moving one will not interfere with others
    FOR(c,slg.corners.size()) { // for each corner 
        // 1. Check if all boundary edges adjacent to this corner are on feature edges. If so, skip (continue with next corner)
        // 2. Store set of adjacent charts
        // 3. Compute average coordinates of the next vertex along each boundary. Maybe also take into account 2nd vertices.
        // 4. Find the nearest vertex of these average coordinates. If still the same vertex, skip (continue with next corner)
        // 5. Store set of adjacent charts of the nearest vertex. If different charts, skip (continue with next corner), because the new corner position will break validity
        // 6. Adjust labeling around the vertex found, so that the next call of StaticLabelingGraph::fill_from() will move the corner and neighboring boundary edges

        const Corner& current_corner = slg.corners[c];

        // 1.

        if(current_corner.all_adjacent_boundary_edges_are_on_feature_edges(mesh,feature_edges)) {
            fmt::println(Logger::out("monotonicity"),"Corner {} ignored because surrounded by feature edges",c); Logger::out("monotonicity").flush();
            continue; // do not move this corner
        }
        // TODO if some adjacent boundary edges are on feature edges but not all,
        // constraint the new position of the corner to be on a feature edge

        // 2.

        slg.get_adjacent_charts_of_vertex(current_corner.vertex,adj_facets,adjacent_charts_of_corner);

        // 3.

        average_coordinates = current_corner.average_coordinates_of_neighborhood(mesh,slg,false,1);

        // 4.

        nearest_vertex_of_average_coordinates = get_nearest_vertex_of_coordinates(mesh_he,adj_facets,average_coordinates,current_corner.vertex,1);
        if(nearest_vertex_of_average_coordinates == current_corner.vertex) {
            fmt::println(Logger::out("monotonicity"),"Corner {} ignored because new position is on the same vertex",c); Logger::out("monotonicity").flush();
            continue;
        }

        // 5.
        
        slg.get_adjacent_charts_of_vertex(nearest_vertex_of_average_coordinates,adj_facets,adjacent_charts_of_new_vertex);
        if(adjacent_charts_of_new_vertex != adjacent_charts_of_corner) {
            fmt::println(Logger::out("monotonicity"),"Corner {} ignored because new position will break validity",c); Logger::out("monotonicity").flush();
            continue;
        }

        // 6.

        // Because we impose a max distance of 1 for `get_nearest_vertex_of_coordinates()`,
        // and beacause `nearest_vertex_of_average_coordinates` != `current_corner.vertex`,
        // we know there is only one edge between `nearest_vertex_of_average_coordinates` and `current_corner.vertex`
        displacement_halfedge = get_halfedge_between_two_vertices(mesh_he,adj_facets,current_corner.vertex,nearest_vertex_of_average_coordinates);
        geo_assert(halfedge_vertex_index_from(mesh,displacement_halfedge) == current_corner.vertex);
        geo_assert(halfedge_vertex_index_to(mesh,displacement_halfedge) == nearest_vertex_of_average_coordinates);

        mesh_he.move_clockwise_around_vertex(displacement_halfedge,true);
        label[halfedge_facet_left(mesh,displacement_halfedge)] = label[halfedge_facet_right(mesh,displacement_halfedge)]; // copy label
        mesh_he.move_counterclockwise_around_vertex(displacement_halfedge,true);
        mesh_he.move_counterclockwise_around_vertex(displacement_halfedge,true);
        label[halfedge_facet_right(mesh,displacement_halfedge)] = label[halfedge_facet_left(mesh,displacement_halfedge)]; // copy label

        nb_corner_moved++;
    }

    fmt::println(Logger::out("monotonicity"),"{} corners moved",nb_corner_moved); Logger::out("monotonicity").flush();
}

bool merge_a_turning_point_and_its_closest_corner(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets) {
    
    if(slg.non_monotone_boundaries.empty()) {
        return false;
    }
    if(feature_edges.empty()) {
        return false;
    }

    // 1. Parse all non-monotone boundaries
    // 2. For each of them, parse turning-points and count how many are on feature edges
    //     Didn't find one -> skip this boundary
    //     Found at least one -> store total number and continue with the first one
    //     Should we expect only one turning-points per non-monotone boundary?
    // 3. Find which of the two corners of the boundary is the closest to the turning point
    // 4. From the turning point "orientation" (toward left or right),
    //    find which boundary to move among the ones around the closest corner
    // 5. Find which label to apply between the new boundary and the boundary to move
    // 6. Move edge by edge from the turning point following the same direction as the boundary to move,
    //    and store facets at left/right of this path
    // 7. Facet by facet, change the labels on one side of the path, effectively moving the boundary

    // 1. 

    FOR(non_monotone_boundary_index,slg.non_monotone_boundaries.size()) { // for each non-monotone boundary
        const Boundary& non_monotone_boundary = slg.boundaries[slg.non_monotone_boundaries[non_monotone_boundary_index]];

        // 2. count turning-points which are on feature edges

        unsigned int nb_turning_points_on_feature_edges = 0;
        TurningPoint first_turning_point_on_feature_edge = non_monotone_boundary.turning_points[0];
        index_t ingoing_local_halfedge_index = index_t(-1);
        index_t outgoing_local_halfedge_index = index_t(-1);
        MeshHalfedges::Halfedge ingoing_local_halfedge;
        MeshHalfedges::Halfedge outgoing_local_halfedge;
        for(const TurningPoint& turning_point : non_monotone_boundary.turning_points) {
            outgoing_local_halfedge_index = turning_point.outgoing_local_halfedge_index_;
            ingoing_local_halfedge_index = (outgoing_local_halfedge_index-1) % (index_t) non_monotone_boundary.halfedges.size();
            ingoing_local_halfedge = non_monotone_boundary.halfedges[ingoing_local_halfedge_index];
            outgoing_local_halfedge = non_monotone_boundary.halfedges[outgoing_local_halfedge_index];
            if(halfedge_is_on_feature_edge(mesh,ingoing_local_halfedge,feature_edges) || halfedge_is_on_feature_edge(mesh,outgoing_local_halfedge,feature_edges)) {
                nb_turning_points_on_feature_edges++;
                if(nb_turning_points_on_feature_edges == 1) {
                    first_turning_point_on_feature_edge = turning_point;
                }
            }
        }

        if(nb_turning_points_on_feature_edges == 0) {
            // all turning-points are on smooth edges
            continue; // check next non-monotone boundary
        }

        // get labeling attribute and instanciate the halfedges API
        Attribute<index_t> label(mesh.facets.attributes(), attribute_name);
        CustomMeshHalfedges mesh_he(mesh);
        mesh_he.set_use_facet_region(attribute_name);
        index_t first_turning_point_on_feature_edge_vertex = first_turning_point_on_feature_edge.vertex(non_monotone_boundary,mesh);

        // 3. Only process `first_turning_point_on_feature_edge`. Find which corner of `non_monotone_boundary` is the closest.

        #ifndef NDEBUG
            dump_vertex("current_tp",mesh_vertex(mesh,first_turning_point_on_feature_edge_vertex));
        #endif
        
        index_t closest_corner = first_turning_point_on_feature_edge.get_closest_corner(non_monotone_boundary,mesh_he);
        #ifndef NDEBUG
            dump_vertex("closest_corner",mesh,slg.corners[closest_corner].vertex);
        #endif

        // 4. Find which boundary around `closest_corner` is closest to `first_turning_point_on_feature_edge`, excluding `non_monotone_boundary`
        //    Also compute the vector between its corners

        const Boundary& boundary_to_move = slg.boundaries[non_monotone_boundary.get_closest_boundary_of_turning_point(first_turning_point_on_feature_edge,closest_corner,mesh_he,slg.halfedge2boundary,slg.corners)];
        #ifndef NDEBUG
            dump_boundary_with_halfedges_indices("boundary_to_move",mesh,boundary_to_move);
        #endif
        vec3 boundary_to_move_vector = boundary_to_move.vector_between_corners(mesh,slg.corners);
        // flip the vector if `boundary_to_move.end_corner` is adjacent to `non_monotone_boundary` and not `boundary_to_move.start_corner`
        if(slg.corners[boundary_to_move.start_corner].vertex != slg.corners[closest_corner].vertex) {
            geo_assert(slg.corners[boundary_to_move.end_corner].vertex == slg.corners[closest_corner].vertex);
            boundary_to_move_vector *= -1.0;
        }
        #ifndef NDEBUG
            dump_vector("direction_for_new_boundary",mesh,first_turning_point_on_feature_edge.vertex(non_monotone_boundary,mesh_he.mesh()),boundary_to_move_vector);
        #endif

        // 5. find which label to assign between `boundary_to_move` and its "parallel" passing by the turning point

        // find the adjacent chart in common with `non_monotone_boundary` and 
        index_t chart_on_which_the_new_boundary_will_be = (
            ((non_monotone_boundary.left_chart == boundary_to_move.left_chart) || (non_monotone_boundary.left_chart == boundary_to_move.right_chart) ) ? non_monotone_boundary.left_chart : (
            ((non_monotone_boundary.right_chart == boundary_to_move.left_chart) || (non_monotone_boundary.right_chart == boundary_to_move.right_chart) ) ? non_monotone_boundary.right_chart : (
            index_t(-1) // in case no match
        )));
        geo_assert(chart_on_which_the_new_boundary_will_be != index_t(-1));
        const std::set<index_t>& facets_of_the_chart_on_which_the_new_boundary_will_be = slg.charts[chart_on_which_the_new_boundary_will_be].facets;
        #ifndef NDEBUG
            dump_facets("facets_of_the_chart_on_which_the_new_boundary_will_be",mesh,facets_of_the_chart_on_which_the_new_boundary_will_be);
        #endif
        index_t new_label = index_t(-1);
        if(chart_on_which_the_new_boundary_will_be == boundary_to_move.right_chart) {
            new_label = slg.charts[boundary_to_move.left_chart].label;
        }
        else {
            new_label = slg.charts[boundary_to_move.right_chart].label;
        }

        // 6. Start from the turning-point, move edge by edge in the direction of `boundary_to_move_vector`

        std::vector<MeshHalfedges::Halfedge> path;
        std::set<index_t> facets_at_left, facets_at_right;

        // quick fix for MAMBO B39-like models
        // after auto_fix_validity(): some corners are misplaced & turning-points are very close to one of the corners,
        // separated by a lost feature edge
        MeshHalfedges::Halfedge halfedge_on_lost_feature_edge;
        if(vertex_has_lost_feature_edge_in_neighborhood(mesh_he,adj_facets,feature_edges,first_turning_point_on_feature_edge_vertex,halfedge_on_lost_feature_edge)) {
            if(dot(normalize(halfedge_vector(mesh,halfedge_on_lost_feature_edge)),normalize(boundary_to_move_vector)) < 0.5) { // if the direction of the lost feature edge is far from the direction of the boundary to move
                goto default_behavior; // S9 also has a lost feature edge but we must use the default behavior
            }
            follow_feature_edge_on_chart(mesh_he,halfedge_on_lost_feature_edge,feature_edges,slg.facet2chart,facets_at_left,facets_at_right);
            const std::set<index_t>& facets_to_re_label = 
                (closest_corner == non_monotone_boundary.end_corner) ?
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_left : facets_at_right) : 
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_right : facets_at_left);
            const std::set<index_t>& wall =
                (closest_corner == non_monotone_boundary.end_corner) ?
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_right : facets_at_left) :
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_left : facets_at_right);
            propagate_label(mesh,attribute_name,new_label,facets_to_re_label,wall,slg.facet2chart,chart_on_which_the_new_boundary_will_be);
            return true;
        }

    default_behavior:

        trace_path_on_chart(mesh_he,adj_facets,slg.facet2chart,slg.turning_point_vertices,first_turning_point_on_feature_edge_vertex,boundary_to_move_vector,facets_at_left,facets_at_right,path);
        #ifndef NDEBUG
            dump_facets("facets_at_left",mesh,facets_at_left);
            dump_facets("facets_at_right",mesh,facets_at_right);
        #endif

        // 7. change the labeling on one side of the just-traced path

        bool the_wall_is_left_facets = (first_turning_point_on_feature_edge.is_towards_left() && closest_corner == non_monotone_boundary.end_corner) || 
                                    (first_turning_point_on_feature_edge.is_towards_right() && closest_corner == non_monotone_boundary.start_corner);
        const std::set<index_t>& wall_facets = the_wall_is_left_facets ? facets_at_left : facets_at_right;
        const std::set<index_t>& facets_to_process = the_wall_is_left_facets ? facets_at_right : facets_at_left;

        // Issue : by propagating the `new_label` on `facets_to_process`, we could make the adjacent chart invalid if it is narrow
        // see MAMBO S9 for example
        // https://gitlab.com/franck.ledoux/mambo
        // Solution : Go through all vertices in `path` except the extremities, and apply the supporting chart label on all adjacent facets
        index_t vertex = index_t(-1);
        for(auto halfedge = path.begin()+1 ; halfedge != path.end(); ++halfedge) {
            vertex = halfedge_vertex_index_from(mesh,*halfedge);
            for(index_t facet : adj_facets[vertex]) {
                label[facet] = slg.charts[chart_on_which_the_new_boundary_will_be].label;
            }
        }

        propagate_label(mesh,attribute_name,new_label,facets_to_process,wall_facets,slg.facet2chart,chart_on_which_the_new_boundary_will_be);

        return true;

    }
    // so none of the non-monotone boundary can be processed
    return false;
}

// returns true if a chart has been created
// returns false if any of the turning-points can be processed with this operator
bool join_turning_points_pair_with_new_chart(GEO::Mesh& mesh, const char* attribute_name, StaticLabelingGraph& slg, const std::vector<vec3>& normals, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets) {
    if(slg.non_monotone_boundaries.empty()) {
        return false;
    }
    if(feature_edges.empty()) {
        return false; // nothing to do, there are no feature edges
    }
    geo_assert(!adj_facets.empty()); // `adj_facets` must have been filled in a parent function
    geo_assert(!normals.empty()); // `normals` must have been filled in a parent function

    // Parse all turning-points
    // Look at their neighboring halfedges
    // If there is an outgoing feature edge not captured (ie same label on both sides),
    // follow the feature edge.
    // If we arrive at another turning point, trace a chart between the two turning-points
    // Else do nothing

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);
    CustomMeshHalfedges mesh_he(mesh);
    mesh_he.set_use_facet_region(attribute_name);
    MeshHalfedges::Halfedge halfedge;
    index_t chart_on_which_the_lost_feature_edge_is = index_t(-1);
    index_t label_on_which_the_lost_feature_edge_is = index_t(-1);
    std::set<index_t> facets_at_left;
    std::set<index_t> facets_at_right;
    index_t vertex_at_tip_of_feature_edge = index_t(-1);
    index_t new_label = index_t(-1);
    index_t adjacent_facet = index_t(-1);
    index_t boundary_of_the_second_turning_point = index_t(-1);

    for(index_t b : slg.non_monotone_boundaries) { // for each non-monotone boundary (b is an index of slg.boundaries)
        for(const auto& tp : slg.boundaries[b].turning_points) { // for each turning-point of the current boundary
            if(vertex_has_lost_feature_edge_in_neighborhood(mesh_he,adj_facets,feature_edges,tp.vertex(slg.boundaries[b],mesh),halfedge)) {
                chart_on_which_the_lost_feature_edge_is = slg.facet2chart[halfedge_facet_left(mesh,halfedge)];
                geo_assert(chart_on_which_the_lost_feature_edge_is == slg.facet2chart[halfedge_facet_right(mesh,halfedge)]);
                label_on_which_the_lost_feature_edge_is = slg.charts[chart_on_which_the_lost_feature_edge_is].label;
                halfedge = follow_feature_edge_on_chart(mesh_he,halfedge,feature_edges,slg.facet2chart,facets_at_left,facets_at_right);
                // so move unsuccessful (no more halfedges on feature edge), or we left the chart (found a boundary / corner / turning-point)
                vertex_at_tip_of_feature_edge = halfedge_vertex_index_to(mesh,halfedge);
                if(slg.turning_point_vertices.contains(vertex_at_tip_of_feature_edge)) {
                    geo_assert(slg.turning_point_vertices[vertex_at_tip_of_feature_edge].size() == 1); // assert only one turning-point associated to this vertex
                    boundary_of_the_second_turning_point = slg.non_monotone_boundaries[slg.turning_point_vertices[vertex_at_tip_of_feature_edge][0].first];
                    new_label = find_optimal_label(
                        {},
                        { // forbidden labels
                            label_on_which_the_lost_feature_edge_is,
                            slg.boundaries[b].other_label(slg.charts,label_on_which_the_lost_feature_edge_is),
                            slg.boundaries[boundary_of_the_second_turning_point].other_label(slg.charts,label_on_which_the_lost_feature_edge_is)
                        }, 
                        {},
                        (average_facets_normal(normals,facets_at_left) + average_facets_normal(normals,facets_at_right)) / 2.0 // new label must be close to the facet normals
                    );
                    bool keep_left_facet = average_dot_product(normals,facets_at_left,label2vector[new_label]) > average_dot_product(normals,facets_at_right,label2vector[new_label]);
                    const std::set<index_t>& facets_to_edit = keep_left_facet ? facets_at_left : facets_at_right;
                    const std::set<index_t>& wall = keep_left_facet ? facets_at_right : facets_at_left;
                    // change label to `new_label` for `facets_to_edit` and their adjacent facets
                    for(auto f : facets_to_edit) {
                        label[f] = new_label;
                        FOR(le,3) { // for each local edge of facet f
                            adjacent_facet = mesh.facets.adjacent(f,le);
                            if(wall.contains(adjacent_facet)) {
                                continue;
                            }
                            if(slg.facet2chart[adjacent_facet] == chart_on_which_the_lost_feature_edge_is) {
                                label[adjacent_facet] = new_label; // also change the label of the adjacent facet
                            }
                        }
                    }
                    return true;
                }
                // else: clear variables & continue
                halfedge = MeshHalfedges::Halfedge(NO_FACET,NO_CORNER);
                chart_on_which_the_lost_feature_edge_is = index_t(-1);
                label_on_which_the_lost_feature_edge_is = index_t(-1);
                facets_at_left.clear();
                facets_at_right.clear();
                vertex_at_tip_of_feature_edge = index_t(-1);
                new_label = index_t(-1);
                adjacent_facet = index_t(-1);
                boundary_of_the_second_turning_point = index_t(-1);
            }
            // else : continue
        }
    }

    return false;
}

bool auto_fix_monotonicity(Mesh& mesh, const char* attribute_name, StaticLabelingGraph& slg, std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<vec3>& normals) {
    
    size_t nb_processed = 0;

    if (slg.non_monotone_boundaries.empty()) {
        return true;
    }

    // compute vertex-to-facet adjacency if not already done
    // required for
    //   join_turning_points_pair_with_new_chart(),
    //   merge_a_turning_point_and_its_closest_corner(),
    //   straighten_boundaries()
    
    if(adj_facets.empty()) {
        compute_adjacent_facets_of_vertices(mesh,adj_facets);
    }

    do {
        nb_processed = (size_t) join_turning_points_pair_with_new_chart(mesh,attribute_name,slg,normals,feature_edges,adj_facets);
        slg.fill_from(mesh,attribute_name,feature_edges);
    } while (nb_processed != 0);

    if (slg.non_monotone_boundaries.empty()) {
        return true;
    }

    size_t nb_iter = 0;
    do {
        nb_processed = (size_t) merge_a_turning_point_and_its_closest_corner(mesh,attribute_name,slg,feature_edges,adj_facets);
        slg.fill_from(mesh,attribute_name,feature_edges);
        nb_iter++;
        if(nb_iter > 20) {
            break;
        }
    } while (nb_processed != 0);

    if (slg.non_monotone_boundaries.empty()) {
        return true;
    }

    move_boundaries_near_turning_points(mesh,attribute_name,slg,feature_edges);
    slg.fill_from(mesh,attribute_name,feature_edges);

    if (slg.non_monotone_boundaries.empty()) {
        return true;
    }

    straighten_boundaries(mesh,attribute_name,slg,adj_facets,feature_edges);
    // `slg` already updated in straighten_boundaries()

    if (slg.non_monotone_boundaries.empty()) {
        return true;
    }

    fmt::println(Logger::out("monotonicity"),"auto fix monotonicity stopped didn't reach all-monotone boundaries"); Logger::out("monotonicity").flush();
    return false;
}

unsigned int count_lost_feature_edges(const CustomMeshHalfedges& mesh_he, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    // parse all facet
    // parse all local vertex for each facet (= facet corners)
    // get halfedge
    // ignore if v1 < v0 (they are 2 halfeges for a given edge, keep the one where v0 < v1)
    // ignore if halfedge not on feature edge
    // fetch label at left and right
    // if same labeling, increment counter

    unsigned int nb_lost_feature_edges = 0;
    geo_assert(mesh_he.is_using_facet_region()); // expecting the labeling to be bounded in `mesh_he`
    const Mesh& mesh = mesh_he.mesh();
    MeshHalfedges::Halfedge current_halfedge;
    FOR(f,mesh.facets.nb()) {
        FOR(lv,3) { // for each local vertex of the current facet
            current_halfedge.facet = f;
            current_halfedge.corner = mesh.facets.corner(f,lv);
            geo_assert(mesh_he.halfedge_is_valid(current_halfedge));
            if(halfedge_vertex_index_to(mesh,current_halfedge) < halfedge_vertex_index_from(mesh,current_halfedge)) {
                continue;
            }
            if(!halfedge_is_on_feature_edge(mesh,current_halfedge,feature_edges)) {
                continue; // not on feature edge
            }
            if(!mesh_he.halfedge_is_border(current_halfedge)) {
                // so our `current_halfedge`
                // - has v0 < v1 (to count each undirected edge once)
                // - is on a feature edge
                // - has the same label at left and right
                nb_lost_feature_edges++;
            }
        }
    }
    geo_assert(nb_lost_feature_edges <= feature_edges.size());
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