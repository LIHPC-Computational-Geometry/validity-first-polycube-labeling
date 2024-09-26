#include <ranges>           // for std::ranges::sort()

#include "labeling_operators_on_invalidity.h"
#include "labeling.h"
#include "containers_macros.h" // for VECTOR_CONTAINS(), MAP_CONTAINS()
#include "labeling_graphcuts.h"
#include "io_dump.h"

size_t remove_surrounded_charts(Attribute<index_t>& labeling, const LabelingGraph& lg) {
    geo_debug_assert(labeling.is_bound());

    // Get charts surrounded by only 1 label
    // Broader fix than just looking at charts having 1 boundary

    size_t nb_invalid_charts_processed = 0;
    for(index_t chart_index : lg.invalid_charts) {
        auto boundary_iterator = lg.charts[chart_index].boundaries.begin();
        index_t chart_at_other_side = lg.boundaries[*boundary_iterator].chart_at_other_side(chart_index);
        index_t surrounding_label = lg.charts[chart_at_other_side].label;
        if(lg.boundaries[*boundary_iterator].on_feature_edge) {
            goto skip_modification; // by propagating the `surrounding_label`, we will lost a feature edge
        }

        boundary_iterator++; // go to next boundary
        for(;boundary_iterator != lg.charts[chart_index].boundaries.end(); ++boundary_iterator) {
            if(lg.boundaries[*boundary_iterator].on_feature_edge) {
                goto skip_modification; // by propagating the `surrounding_label`, we will lost a feature edge
            }
            chart_at_other_side = lg.boundaries[*boundary_iterator].chart_at_other_side(chart_index);
            if(surrounding_label != lg.charts[chart_at_other_side].label) {
                goto skip_modification; // this chart has several labels around (goto because break in nested loops is a worse idea)
            }
        }

        // if we are here, the goto was not used, so the only label around is surrounding_label
        for(index_t facet_index : lg.charts[chart_index].facets) {
            labeling[facet_index] = surrounding_label;
        }
        nb_invalid_charts_processed++;


        skip_modification:
            ; // just end this loop interation
    }

    return nb_invalid_charts_processed;
}

bool fix_an_invalid_boundary(
    const MeshExt& mesh,
    Attribute<index_t>& labeling,
    LabelingGraph& lg
) {
    geo_debug_assert(labeling.is_bound());
    geo_assert(mesh.halfedges.is_using_facet_region());

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
    for(index_t boundary_index : lg.invalid_boundaries) { // for each invalid boundary

        // get ref to boundary object and retrieve geometry and neighborhood info

        const Boundary& current_boundary = lg.boundaries[boundary_index];
        MeshHalfedges::Halfedge an_halfedge_of_the_invalid_boundary = current_boundary.halfedges[0];
        const Chart& left_chart = lg.charts[current_boundary.left_chart];
        const Chart& right_chart = lg.charts[current_boundary.right_chart];
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
        current_boundary.get_adjacent_facets(mesh,left_facets_along_boundary,OnlyLeft,lg.facet2chart,1); // get triangles at left and at distance 0 or 1 from the boundary
        current_boundary.get_adjacent_facets(mesh,right_facets_along_boundary,OnlyRight,lg.facet2chart,1); // get triangles at right and at distance 0 or 1 from the boundary

        // Find out if the chart to insert on the boundary is better suited on the left side or the right side
        // Don't manage case where both sides are equally suited,
        // because the boundary could be on a feature edge, and we shouldn't de-capture it
        
        double average_dot_product_left = average_angle(mesh.facet_normals.as_vector(),left_facets_along_boundary,new_label_as_vector);
        double average_dot_product_right = average_angle(mesh.facet_normals.as_vector(),right_facets_along_boundary,new_label_as_vector);
        // comparison
        if (average_dot_product_left < average_dot_product_right) { // on average, the left side is a better place to put the new chart
            for(index_t f : left_facets_along_boundary) {
                labeling[f] = new_label;
            }
            a_facet_on_new_chart = *left_facets_along_boundary.begin();
        }
        else { // on average, the right side is a better place to put the new chart
            for(index_t f : right_facets_along_boundary) {
                labeling[f] = new_label;
            }
            a_facet_on_new_chart = *right_facets_along_boundary.begin();
        }

        // If the new chart is invalid
        // by having only one adjacent chart (the one at the other side of what was the invalid boundary)
        // extend it from the two corners until we found boundaries
        // and change the label of all facets below
        // See MAMBO B13 model for example
        // https://gitlab.com/franck.ledoux/mambo

        geo_assert(mesh.adj_facet_corners.size_matches_nb_vertices()); // we need adjacency between vertices and facets for this part

        lg.fill_from(mesh,labeling);
        // should mesh.halfedges regions be updated?

        geo_assert(lg.boundaries[lg.halfedge2boundary[an_halfedge_of_the_invalid_boundary].first].axis != -1);
        index_t axis_of_just_fixed_boundary = (index_t) lg.boundaries[lg.halfedge2boundary[an_halfedge_of_the_invalid_boundary].first].axis;

        const Chart& created_chart = lg.charts[lg.facet2chart[a_facet_on_new_chart]];

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
        if(lg.boundaries[*created_chart.boundaries.begin()].turning_points.size() == 2) {
            non_monotone_boundary = *created_chart.boundaries.begin();
            tp0 = lg.boundaries[non_monotone_boundary].turning_points[0];
            tp1 = lg.boundaries[non_monotone_boundary].turning_points[1];
        }
        else if(lg.boundaries[*created_chart.boundaries.rbegin()].turning_points.size() == 2) {
            non_monotone_boundary = *created_chart.boundaries.rbegin();
            tp0 = lg.boundaries[non_monotone_boundary].turning_points[0];
            tp1 = lg.boundaries[non_monotone_boundary].turning_points[1];
        }
        else {
            fmt::println(Logger::err("fix validity"),"In fix_as_much_invalid_boundaries_as_possible(), did not find the boundary with 2 turning-points in the contour of the created chart");
            continue;
        }
        const TurningPoint& turning_point_at_max_coordinate_on_axis = 
            mesh_vertex(mesh,tp0.vertex(lg.boundaries[non_monotone_boundary],mesh))[axis_of_just_fixed_boundary] >
            mesh_vertex(mesh,tp1.vertex(lg.boundaries[non_monotone_boundary],mesh))[axis_of_just_fixed_boundary] ?
            tp0 : tp1;
        const TurningPoint& turning_point_at_min_coordinate_on_axis = turning_point_at_max_coordinate_on_axis == tp0 ?
            tp1 : tp0;

        std::vector<MeshHalfedges::Halfedge> path;
        std::set<index_t> facets_at_left;
        std::set<index_t> facets_at_right;
        trace_path_on_chart(
            mesh,
            lg.facet2chart,
            lg.turning_point_vertices,
            turning_point_at_max_coordinate_on_axis.vertex(lg.boundaries[non_monotone_boundary],mesh),
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
                labeling,
                created_chart.label,
                facets_at_left,
                facets_at_right,
                lg.facet2chart,
                lg.facet2chart[*facets_at_left.begin()]
            );
        }
        else {
            // facets_at_right -> new chart
            // facets_at_left -> wall
            propagate_label(
                mesh,
                labeling,
                created_chart.label,
                facets_at_right,
                facets_at_left,
                lg.facet2chart,
                lg.facet2chart[*facets_at_right.begin()]
            );
        }

        path.clear();
        facets_at_left.clear();
        facets_at_right.clear();
        trace_path_on_chart(
            mesh,
            lg.facet2chart,
            lg.turning_point_vertices,
            turning_point_at_min_coordinate_on_axis.vertex(lg.boundaries[non_monotone_boundary],mesh),
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
                labeling,
                created_chart.label,
                facets_at_left,
                facets_at_right,
                lg.facet2chart,
                lg.facet2chart[*facets_at_left.begin()]
            );
        }
        else {
            // facets_at_right -> new chart
            // facets_at_left -> wall
            propagate_label(
                mesh,
                labeling,
                created_chart.label,
                facets_at_right,
                facets_at_left,
                lg.facet2chart,
                lg.facet2chart[*facets_at_right.begin()]
            );
        }
        return true;
    }

    return false;
}

size_t fix_as_much_invalid_corners_as_possible(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg) {
    geo_debug_assert(labeling.is_bound());
    geo_assert(mesh.halfedges.is_using_facet_region());

    // Replace the labels around invalid corners
    // by the nearest one from the vertex normal

    unsigned int nb_invalid_corners_processed = 0;
    index_t new_label;
    vec3 vertex_normal = {0,0,0};
    std::vector<MeshHalfedges::Halfedge> outgoing_halfedges_on_feature_edge;
    for(index_t corner_index : lg.invalid_corners) { // for each invalid corner

        // check if there is a feature edge along the corner
        lg.corners[corner_index].get_outgoing_halfedges_on_feature_edge(mesh,outgoing_halfedges_on_feature_edge);
        // also compute the evenness of angles between adj boundaries
        double sd_boundary_angles = lg.corners[corner_index].sd_boundary_angles(mesh);
        if(outgoing_halfedges_on_feature_edge.empty() || sd_boundary_angles < 0.02) {

            // cone-like, or pyramid-like (MAMBO B20) invalid corner
            // compute the normal by adding normals of adjacent facets

            geo_assert(mesh.adj_facet_corners.size_matches_nb_vertices());
            vertex_normal = {0,0,0};
            for(const auto& facet_corner : mesh.adj_facet_corners.of_vertex(lg.corners[corner_index].vertex)) {
                vertex_normal += mesh.facet_normals[facet_corner.facet_index];
            }
            vertex_normal /= (double) mesh.adj_facet_corners.of_vertex(lg.corners[corner_index].vertex).size();

            // compute new label
            new_label = nearest_label(vertex_normal);

            // change the label of adjacent facets
            for(const auto& facet_corner : mesh.adj_facet_corners.of_vertex(lg.corners[corner_index].vertex)) {
                labeling[facet_corner.facet_index] = new_label;
            }

            nb_invalid_corners_processed++;
        }
        else {
            // there is some kind of feature-edges pinching

            MeshHalfedges::Halfedge halfedge;
            if(vertex_has_lost_feature_edge_in_neighborhood(mesh,lg.corners[corner_index].vertex,halfedge)) {

                // There is an outgoing feature edge not captured (ie same label on both sides),
                // Follow the feature edge.
                // If we arrive at another invalid corner, trace a chart between the two invalid corners
                // Else do nothing
                // See MAMBO S9 for example
                // https://gitlab.com/franck.ledoux/mambo/

                index_t current_chart = lg.facet2chart[halfedge_facet_left(mesh,halfedge)];
                index_t label_of_current_chart = lg.charts[current_chart].label;
                index_t new_label = index_t(-1);
                if(label_of_current_chart != labeling[halfedge_facet_right(mesh,halfedge)]) {
                    return nb_invalid_corners_processed; // cannot process this configuration. Seems to arise on MAMBO M8
                }
                std::set<index_t> facets_at_left;
                std::set<index_t> facets_at_right;
                index_t adjacent_facet = index_t(-1);
                halfedge = follow_feature_edge_on_chart(mesh,halfedge,lg.facet2chart,facets_at_left,facets_at_right);
                index_t corner_found = lg.vertex2corner[halfedge_vertex_index_to(mesh,halfedge)]; // can be -1
                if( (corner_found != index_t(-1)) && VECTOR_CONTAINS(lg.invalid_corners,corner_found) ) {
                    // we found a lost feature edge between 2 invalid corners
                    vec3 avg_normal_at_left = average_facets_normal(mesh.facet_normals.as_vector(),facets_at_left);
                    vec3 avg_normal_at_right = average_facets_normal(mesh.facet_normals.as_vector(),facets_at_right);
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
                            labeling[f] = new_label;
                            FOR(le,3) { // for each local edge of facet f
                                adjacent_facet = mesh.facets.adjacent(f,le);
                                if(facets_at_right.contains(adjacent_facet)) {
                                    continue;
                                }
                                if(lg.facet2chart[adjacent_facet] == current_chart) {
                                    labeling[adjacent_facet] = new_label; // also change the label of the adjacent facet
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
                            labeling[f] = new_label;
                            FOR(le,3) { // for each local edge of facet f
                                adjacent_facet = mesh.facets.adjacent(f,le);
                                if(facets_at_left.contains(adjacent_facet)) {
                                    continue;
                                }
                                if(lg.facet2chart[adjacent_facet] == current_chart) {
                                    labeling[adjacent_facet] = new_label; // also change the label of the adjacent facet
                                }
                            }
                        }
                    }

                    lg.fill_from(mesh,labeling);

                    // we need to stop fix_as_much_invalid_corners_as_possible() now
                    // because the loop iterating over `lg.invalid_corners` is no longer up to date
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
                    index_t b = lg.halfedge2boundary.at(halfedge).first;
                    adjacent_boundaries_and_their_length.push_back(std::make_pair(
                        b, // boundary index
                        lg.boundaries[b].length(mesh) // boundary length
                    ));
                }
                std::ranges::sort(
                    adjacent_boundaries_and_their_length,
                    [](const std::pair<index_t,double>& a, const std::pair<index_t,double>& b) { return a.second < b.second; } // sort by second item in the std::pair (the length)
                );
                // get the 2 shortest boundaries
                const Boundary& b0 = lg.boundaries[adjacent_boundaries_and_their_length[0].first];
                const Boundary& b1 = lg.boundaries[adjacent_boundaries_and_their_length[1].first];
                index_t chart_in_common = adjacent_chart_in_common(b0,b1);
                index_t label_of_chart_in_common = lg.charts[chart_in_common].label;

                std::set<index_t> facets_to_edit;
                b0.get_adjacent_facets(
                    mesh,
                    facets_to_edit,
                    chart_in_common == b0.left_chart ? OnlyRight : OnlyLeft,
                    lg.facet2chart,
                    1 // include facet at a distance of 1 (touch the boundary from a vertex)
                );
                // `facets_to_edit` is not cleared inside get_adjacent_facets(), we can call it a second time on top of the existing set values
                b1.get_adjacent_facets(
                    mesh,
                    facets_to_edit,
                    chart_in_common == b1.left_chart ? OnlyRight : OnlyLeft,
                    lg.facet2chart,
                    1 // include facet at a distance of 1 (touch the boundary from a vertex)
                );
                for(index_t f : facets_to_edit) {
                    labeling[f] = label_of_chart_in_common;
                }
                nb_invalid_corners_processed++;
            }
        }
    }

    return nb_invalid_corners_processed;
}

// from https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/labeling_ops.cpp removeChartMutation()
bool remove_invalid_charts(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg) {
    geo_debug_assert(labeling.is_bound());

    if(lg.invalid_charts.empty()) {
        fmt::println(Logger::warn("fix_labeling"),"remove_invalid_charts canceled because there are no invalid charts"); Logger::warn("fix_labeling").flush();
        return false;
    }

    // Fill invalid charts using a Graph-Cut optimization,
    // preventing existing label from being re-applied

    // compactness = 1
    GraphCutLabeling gcl(mesh);
    gcl.data_cost__set__locked_labels(labeling); // start by locking all the labels, so valid charts will not be modified
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(1);

    unsigned int nb_charts_to_remove = 0;
    for(index_t chart_index : lg.invalid_charts) { // for each invalid chart

        // prevent removal of chart surrounded by feature edges
        // TODO : increase their valence, by boundary duplication or boundary shift
        if(lg.charts[chart_index].is_surrounded_by_feature_edges(lg.boundaries)) {
            fmt::println(Logger::warn("fix_labeling"),"Cannot remove chart n°{} because it is surrounded by feature edges",chart_index); Logger::warn("fix_labeling").flush();
            continue;
        }

        for(index_t facet_index : lg.charts[chart_index].facets) { // for each facet inside this chart
            gcl.data_cost__change_to__fidelity_based(facet_index,1);
            // if facet next to a boundary, lower the cost of assigning the neighboring label
            FOR(le,3) { // for each local edge
                index_t adjacent_facet = mesh.facets.adjacent(facet_index,le);
                if(labeling[adjacent_facet] != labeling[facet_index]) {
                    gcl.data_cost__change_to__scaled(facet_index,labeling[adjacent_facet],0.5f); // halve the cost
                }
            }
            gcl.data_cost__change_to__forbidden_polycube_label(facet_index,labeling[facet_index]); // prevent the label from staying the same
        }
        nb_charts_to_remove++;
    }

    if(nb_charts_to_remove == 0) {
        // all invalid charts are surrounded by feature edges
        return true;
    }

    gcl.compute_solution(labeling);
    return false;
}

void remove_charts_around_invalid_boundaries(const MeshExt& mesh, Attribute<index_t>& labeling, const LabelingGraph& lg) {
    geo_debug_assert(labeling.is_bound());
    
    if(lg.invalid_boundaries.empty()) {
        fmt::println(Logger::out("fix_labeling"),"Warning : operation canceled because there are no invalid boundaries"); Logger::out("fix_labeling").flush();
        return;
    }

    // In involved charts, preventing existing label from being re-applied
    // Lock labels in other charts

    std::set<index_t> charts_to_remove;
    for(index_t b : lg.invalid_boundaries) {
        charts_to_remove.insert(lg.boundaries[b].left_chart);
        charts_to_remove.insert(lg.boundaries[b].right_chart);
    }

    GraphCutLabeling gcl(mesh);
    gcl.data_cost__set__locked_labels(labeling); // start by locking all the labels, so other charts will not be modified
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(1);

    for(index_t chart_index : charts_to_remove) { // for each chart to remove

        for(index_t facet_index : lg.charts[chart_index].facets) { // for each facet inside this chart
            gcl.data_cost__change_to__fidelity_based(facet_index,1);
            // if facet next to a boundary, lower the cost of assigning the neighboring label
            FOR(le,3) { // for each local edge
                index_t adjacent_facet = mesh.facets.adjacent(facet_index,le);
                if(labeling[adjacent_facet] != labeling[facet_index]) {
                    gcl.data_cost__change_to__scaled(facet_index,labeling[adjacent_facet],0.5f); // halve the cost
                }
            }
            gcl.data_cost__change_to__forbidden_polycube_label(facet_index,labeling[facet_index]); // prevent the label from staying the same
        }
    }

    gcl.compute_solution(labeling);

}

bool increase_chart_valence(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg) {
    geo_debug_assert(labeling.is_bound());

    FOR(invalid_chart_index,lg.invalid_charts.size()) {
        index_t chart_index = lg.invalid_charts[invalid_chart_index];
        const Chart& chart = lg.charts[chart_index];
        if(!chart.is_surrounded_by_feature_edges(lg.boundaries)) {
            continue; // check next invalid chart
        }
        geo_assert(mesh.halfedges.is_using_facet_region());

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
        chart.counterclockwise_boundaries_order(mesh.halfedges,lg.halfedge2boundary,lg.boundaries,counterclockwise_order);

        FOR(lb,counterclockwise_order.size()) { // for each local boundary index (relative to the chart contour)
            auto [b,b_is_same_direction] = counterclockwise_order[lb]; // b is a boundary index
            if( !lg.boundaries[b].turning_points.empty() && (lg.boundaries[b].axis != -1) ) {
                // naive labeling of MAMBO S36 has 2 turning-points around one of its invalid chart surrounded by feature edges
                // only one of them is relevant for this operator -> the closest to the boundary midpoint
                index_t ltp = 0; // local turning-point index
                if(lg.boundaries[b].turning_points.size() > 1) {
                    size_t max_dist_to_corners = 0; // in number of edges
                    FOR(current_ltp,lg.boundaries[b].turning_points.size()) { // for each turning-point of the boundary `b`
                        size_t current_dist = std::min(
                            (size_t) lg.boundaries[b].turning_points[current_ltp].outgoing_local_halfedge_index_,
                            lg.boundaries[b].halfedges.size() - ((size_t) lg.boundaries[b].turning_points[current_ltp].outgoing_local_halfedge_index_) + 1
                        );
                        if(current_dist > max_dist_to_corners) {
                            max_dist_to_corners = current_dist;
                            ltp = current_ltp;
                        }
                    }
                }
                current_axis = (index_t) lg.boundaries[b].axis;
                problematic_non_monotone_boundary = b;
                Boundary current_boundary_as_counterclockwise;
                if(b_is_same_direction) {
                    current_boundary_as_counterclockwise = lg.boundaries[b];
                }
                else {
                    lg.boundaries[b].get_flipped(mesh.halfedges,current_boundary_as_counterclockwise);
                    ltp = ((index_t) lg.boundaries[b].turning_points.size()) - 1 - ltp;
                }
                current_boundary_as_counterclockwise.split_at_turning_point(mesh.halfedges,downward_boundary,upward_boundary,ltp);
                axis_to_insert = nearest_axis_of_edges(mesh,{
                    lg.boundaries[b].halfedges[lg.boundaries[b].turning_points[ltp].outgoing_local_halfedge_index_],
                    lg.boundaries[b].halfedges[lg.boundaries[b].turning_points[ltp].outgoing_local_halfedge_index_-1]
                },{current_axis});
                problematic_vertex = current_boundary_as_counterclockwise.turning_points[ltp].vertex(current_boundary_as_counterclockwise,mesh);
                break;
            }
            geo_assert(lg.boundaries[b].start_corner != index_t(-1));
            geo_assert(lg.boundaries[b].end_corner != index_t(-1));
            auto [next_b,next_b_is_same_direction] = counterclockwise_order[(lb+1) % counterclockwise_order.size()];
            if(lg.boundaries[b].axis == lg.boundaries[next_b].axis) {
                geo_assert(lg.boundaries[b].axis != -1);
                current_axis = (index_t) lg.boundaries[b].axis;
                if(b_is_same_direction) {
                    lg.boundaries[b].get_flipped(mesh.halfedges,downward_boundary);
                    problematic_corner = lg.boundaries[b].end_corner;
                }
                else {
                    downward_boundary = lg.boundaries[b];
                    problematic_corner = lg.boundaries[b].start_corner;
                }
                // if this is a pyramid-like invalid corner (see MAMBO B20),
                // don't apply increase_chart_valence() but fix_as_much_invalid_corners_as_possible()
                double sd_boundary_angles = lg.corners[problematic_corner].sd_boundary_angles(mesh);
                if(
                    (lg.corners[problematic_corner].valence() == 4) && 
                    lg.corners[problematic_corner].all_adjacent_boundary_edges_are_on_feature_edges(mesh) &&
                    sd_boundary_angles < 0.02 // required to distinguish B20 from S36, S33
                ) {
                    problematic_corner = index_t(-1);
                    problematic_vertex = index_t(-1);
                    continue;
                }
                if(next_b_is_same_direction) {
                    upward_boundary = lg.boundaries[next_b];
                }
                else {
                    lg.boundaries[next_b].get_flipped(mesh.halfedges,upward_boundary);
                }
                problematic_vertex = lg.corners[problematic_corner].vertex;
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
        downward_boundary.per_edges_cumulative_axis_assignment_cost(mesh,current_axis,downward_boundary_cumulative_cost_for_current_axis,false);
        downward_boundary.per_edges_cumulative_axis_assignment_cost(mesh,axis_to_insert,downward_boundary_cumulative_cost_for_axis_to_insert,true);
        upward_boundary.per_edges_cumulative_axis_assignment_cost(mesh,current_axis,upward_boundary_cumulative_cost_for_current_axis,false);
        upward_boundary.per_edges_cumulative_axis_assignment_cost(mesh,axis_to_insert,upward_boundary_cumulative_cost_for_axis_to_insert,true);

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
        outgoing_halfedge = get_an_outgoing_halfedge_of_vertex(mesh,downward_boundary_equilibrium_vertex);
        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(mesh,outgoing_halfedge,label2vector[chart.label]);
        double max_angle_same_direction = angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[chart.label]);

        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(mesh,outgoing_halfedge,label2vector[opposite_label(chart.label)]);
        double max_angle_opposite_direction = angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[opposite_label(chart.label)]);

        // neighborhood of `upward_boundary_equilibrium_vertex`
        outgoing_halfedge = get_an_outgoing_halfedge_of_vertex(mesh,upward_boundary_equilibrium_vertex);
        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(mesh,outgoing_halfedge,label2vector[chart.label]);
        max_angle_same_direction = std::max(max_angle_same_direction,angle(normalize(halfedge_vector(mesh,outgoing_halfedge)),label2vector[chart.label]));

        outgoing_halfedge = get_most_aligned_halfedge_around_vertex(mesh,outgoing_halfedge,label2vector[opposite_label(chart.label)]);
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
            MeshHalfedges::Halfedge halfedge = lg.corners[problematic_corner].get_most_aligned_boundary_halfedge(mesh,direction);
            geo_assert(halfedge_vertex_index_from(mesh,halfedge) == downward_boundary_equilibrium_vertex);
            auto [boundary_to_reuse, same_direction] = lg.halfedge2boundary[halfedge];
            geo_assert(boundary_to_reuse != index_t(-1));
            if(same_direction) {
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyLeft,lg.facet2chart);
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyRight,lg.facet2chart);
            }
            else {
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyRight,lg.facet2chart);
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyLeft,lg.facet2chart);
            }
        }
        else {
            // trace path
            trace_path_on_chart(mesh,lg.facet2chart,lg.turning_point_vertices,downward_boundary_equilibrium_vertex,direction*10,facets_of_new_chart,walls,path);
        }

        if( (upward_boundary_equilibrium_vertex == problematic_vertex) && (problematic_corner != index_t(-1)) ) {
            // This equilibrium point is on the problematic corner (between 2 boundaries assigned to the same axis)
            // So there is already an outgoing boundary, we just need to fetch facets at its left and right
            MeshHalfedges::Halfedge halfedge = lg.corners[problematic_corner].get_most_aligned_boundary_halfedge(mesh,direction);
            geo_assert(halfedge_vertex_index_from(mesh,halfedge) == upward_boundary_equilibrium_vertex);
            auto [boundary_to_reuse, same_direction] = lg.halfedge2boundary[halfedge];
            geo_assert(boundary_to_reuse != index_t(-1));
            if(same_direction) {
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyLeft,lg.facet2chart);
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyRight,lg.facet2chart);
            }
            else {
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,walls,OnlyRight,lg.facet2chart);
                lg.boundaries[boundary_to_reuse].get_adjacent_facets(mesh,facets_of_new_chart,OnlyLeft,lg.facet2chart);
            }
        }
        else {
            // trace path
            trace_path_on_chart(mesh,lg.facet2chart,lg.turning_point_vertices,upward_boundary_equilibrium_vertex,direction*10,walls,facets_of_new_chart,path);
        }

        index_t chart_on_which_the_new_chart_will_be = lg.facet2chart[*facets_of_new_chart.begin()];

        #ifndef NDEBUG
            dump_facets("facets_of_new_chart",mesh,facets_of_new_chart);
            dump_facets("walls",mesh,walls);
        #endif

        index_t label_to_insert = find_optimal_label(
            { // 2 forbidden axes:
                chart.label/2, // the axis of the chart we're going to increase the valence
                lg.charts[chart_on_which_the_new_chart_will_be].label/2 // the axis of the chart on which we traced boundaries
            },
            {},
            {}, // no need to specify orthogonality constraints, because we already forbid 2 axes over 3
            average_facets_normal(mesh.facet_normals.as_vector(),facets_of_new_chart) // among the 2 remaining labels, choose the one the closest to the avg normal of `facets_of_new_chart`
        );

        propagate_label(mesh,labeling,label_to_insert,facets_of_new_chart,walls,lg.facet2chart,chart_on_which_the_new_chart_will_be);

        return true; // an invalid chart has been processed
    }
    return false;
}

bool auto_fix_validity(const MeshExt& mesh, Attribute<index_t>& labeling, LabelingGraph& lg, unsigned int max_nb_loop) {
    geo_debug_assert(labeling.is_bound());
    unsigned int nb_loops = 0;
    size_t nb_processed = 0;
    std::set<std::array<std::size_t,7>> set_of_labeling_features_combinations_encountered;
    while(!lg.is_valid() && nb_loops <= max_nb_loop) { // until valid labeling OR too much steps
        nb_loops++;

        // as much as possible, remove isolated (surrounded) charts
        do {
            nb_processed = remove_surrounded_charts(labeling,lg);
            // update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
            lg.fill_from(mesh,labeling);
        } while(nb_processed != 0);

        if(lg.is_valid())
            return true;

        do {
            nb_processed = (size_t) increase_chart_valence(mesh,labeling,lg);
            lg.fill_from(mesh,labeling);
        } while(nb_processed != 0);

        if(lg.is_valid())
            return true;

        unsigned int nb_iter_fix_invalid_boundary = 0; // on some models like MAMBO M8, the loop of this operator is called indefinitely -> add a max number of iterations
        do {
            nb_processed = (size_t) fix_an_invalid_boundary(mesh,labeling,lg);
            lg.fill_from(mesh,labeling);
            nb_iter_fix_invalid_boundary++;
        }
        while ((nb_processed != 0) && nb_iter_fix_invalid_boundary < 20);

        if(lg.is_valid())
            return true;
        
        do {
            nb_processed = fix_as_much_invalid_corners_as_possible(mesh,labeling,lg);
            lg.fill_from(mesh,labeling);
        }
        while (nb_processed != 0);

        if(lg.is_valid())
            return true;

        set_of_labeling_features_combinations_encountered.clear();
        set_of_labeling_features_combinations_encountered.insert({
            lg.nb_charts(),
            lg.nb_boundaries(),
            lg.nb_corners(),
            lg.nb_invalid_charts(),
            lg.nb_invalid_boundaries(),
            lg.nb_invalid_corners(),
            lg.nb_turning_points()
        });

        while(1) {
            remove_invalid_charts(mesh,labeling,lg);
            lg.fill_from(mesh,labeling);

            if(lg.is_valid())
                return true;

            std::array<std::size_t,7> features_combination = {
                lg.nb_charts(),
                lg.nb_boundaries(),
                lg.nb_corners(),
                lg.nb_invalid_charts(),
                lg.nb_invalid_boundaries(),
                lg.nb_invalid_corners(),
                lg.nb_turning_points()
            };

            if(VECTOR_CONTAINS(set_of_labeling_features_combinations_encountered,features_combination)) { // we can use VECTOR_CONTAINS() on sets because they also have find(), cbegin() and cend()
                // we backtracked
                // There is probably small charts that we can remove to help the fixing routine
                remove_charts_around_invalid_boundaries(mesh,labeling,lg);
                lg.fill_from(mesh,labeling);
                break; // go back to the beginning of the loop, with other fix operators
            }
            else {
                set_of_labeling_features_combinations_encountered.insert(features_combination); // store the current combination of number of features
            }
        }
        
    }

    if(!lg.is_valid()) {
        fmt::println(Logger::out("fix_labeling"),"auto fix validity stopped (max nb loops reached), no valid labeling found"); Logger::out("fix_labeling").flush();
        return false;
    }
    
    return true;
}