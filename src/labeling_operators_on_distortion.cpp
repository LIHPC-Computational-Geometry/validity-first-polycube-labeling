#include "labeling_operators_on_distortion.h"
#include "labeling.h"
#include "containers_macros.h" // for MAP_CONTAINS()
#include "io_dump.h"
#include "labeling_graphcuts.h"

#include <deque>

size_t move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges) {

    size_t nb_labels_changed = 0;
    size_t nb_turning_points_moved = 0;

    if(lg.non_monotone_boundaries.empty()) {
        fmt::println(Logger::out("fix_labeling"),"Warning : operation canceled because all boundaries are monotone"); Logger::out("fix_labeling").flush();
        return nb_labels_changed;
    }

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    MeshHalfedgesExt mesh_halfedges(mesh); // create an halfedges interface for this mesh
    mesh_halfedges.set_use_facet_region(attribute_name);
    MeshHalfedges::Halfedge initial_halfedge, current_halfedge;
    index_t new_label = index_t(-1);

    // Get turning points, change the label if this doesnt break charts connectivity

    for(auto b : lg.non_monotone_boundaries) {

        if(lg.boundaries[b].on_feature_edge) {
            continue; // don't move this boundary away from the feature edge
        }

        for(auto tp : lg.boundaries[b].turning_points) {

            initial_halfedge = lg.boundaries[b].halfedges[tp.outgoing_local_halfedge_index_];
            geo_assert(mesh_halfedges.halfedge_is_border(initial_halfedge));
            geo_assert(MAP_CONTAINS(lg.halfedge2boundary,initial_halfedge));

            // if one of the halfedges, the one before or the one after,
            // is on a feature edge, don't move this boundary

            if(halfedge_is_on_feature_edge(mesh,initial_halfedge,feature_edges)) {
                break; // don't move this boundary away from the feature edge
            }
            if(halfedge_is_on_feature_edge(mesh,
                lg.boundaries[b].halfedges[(tp.outgoing_local_halfedge_index_+1) % lg.boundaries[b].halfedges.size()],
                feature_edges
            )) {
                break; // don't move this boundary away from the feature edge
            }

            // get the vertex index
            index_t current_vertex = mesh.facet_corners.vertex(initial_halfedge.corner);
            geo_assert(lg.turning_point_vertices.contains(current_vertex));
            geo_assert(lg.vertex2corner[current_vertex] == index_t(-1)); // a turning point should not be a corner
            // test if the valence of current_vertex is 2
            VertexRingWithBoundaries vr;
            vr.explore(initial_halfedge,mesh_halfedges);
            geo_assert(vr.valence() == 2); // should have only 2 charts

            // new label according to towards which chart the turning point is
            new_label = lg.charts[tp.is_towards_left() ? lg.boundaries[b].left_chart : lg.boundaries[b].right_chart].label;

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

void straighten_boundary_with_GCO(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const LabelingGraph& lg, index_t boundary_index) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    const Boundary& current_boundary = lg.boundaries[boundary_index];
    index_t left_chart_index = current_boundary.left_chart;
    index_t right_chart_index = current_boundary.right_chart;
    const Chart& left_chart = lg.charts[left_chart_index];
    const Chart& right_chart = lg.charts[right_chart_index];
    MeshHalfedgesExt mesh_he(mesh);

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
        dump_all_boundaries("boundaries",mesh,lg.boundaries);

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
    MeshHalfedgesExt::Halfedge current_halfedge;
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

bool straighten_boundary(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, index_t boundary_index, const std::vector<std::vector<index_t>>& adj_facets) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    const Boundary& current_boundary = lg.boundaries[boundary_index];

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
    const Chart& left_chart = lg.charts[left_chart_index];
    const Chart& right_chart = lg.charts[right_chart_index];
    MeshHalfedgesExt mesh_he(mesh);

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
        if(!lg.vertex_is_only_surrounded_by(halfedge_vertex_index_to(mesh,current_halfedge),{left_chart_index,right_chart_index},adj_facets)) {
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

void straighten_boundaries(GEO::Mesh& mesh, const char* attribute_name, LabelingGraph& lg, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    if(lg.boundaries.empty()) {
        fmt::println(Logger::out("monotonicity"),"No boundaries, operation canceled"); Logger::out("monotonicity").flush();
        return;
    }

    geo_assert(!adj_facets.empty())

    // Because the boundaries may not have the same index before and after lg.fill_from()
    // and because we need to call lg.fill_from() after that one boundary is processed to update charts,
    // we first gather the first boundary edge of all boundaries, so at each step we can get the current index of the
    // associated boundary with lg.halfedge2boundary
    // The fist boundary edges should not move after a call of straighten_boundary()...
    std::deque<MeshHalfedges::Halfedge> boundary_edges_to_process;
    boundary_edges_to_process.resize(lg.nb_boundaries()); // preallocation
    FOR(b,lg.boundaries.size()) {
        if (lg.boundaries[b].on_feature_edge) {
            fmt::println(Logger::out("monotonicity"),"Boundary {} skipped for straighten_boundary() because it is on a feature edge",b); Logger::out("monotonicity").flush();
            continue; // do not straighten boundaries surrounded by feature edges
        }
        geo_assert(!lg.boundaries[b].halfedges.empty());
        boundary_edges_to_process.push_back(lg.boundaries[b].halfedges[0]);
    }
    MeshHalfedges::Halfedge current_boundary_edge;
    index_t boundary_index = index_t(-1);
    bool boundary_in_same_direction = false;
    unsigned int count_iterations = 0;
    while (!boundary_edges_to_process.empty()) {
        current_boundary_edge = boundary_edges_to_process.back();
        boundary_edges_to_process.pop_back();
        std::tie(boundary_index,boundary_in_same_direction) = lg.halfedge2boundary[current_boundary_edge];
        geo_assert(boundary_index != index_t(-1));
        if(straighten_boundary(mesh,attribute_name,lg,boundary_index,adj_facets)) 
        {
            lg.fill_from(mesh,attribute_name,feature_edges);
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

void move_corners(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    MeshHalfedgesExt mesh_he(mesh);
    mesh_he.set_use_facet_region(attribute_name);
    std::set<index_t> adjacent_charts_of_corner; // indices of adjacent charts
    std::set<index_t> adjacent_charts_of_new_vertex;
    vec3 average_coordinates;
    index_t nearest_vertex_of_average_coordinates = index_t(-1);
    MeshHalfedges::Halfedge displacement_halfedge;
    unsigned int nb_corner_moved = 0;

    // we can process all corners at once, because moving one will not interfere with others
    FOR(c,lg.corners.size()) { // for each corner 
        // 1. Check if all boundary edges adjacent to this corner are on feature edges. If so, skip (continue with next corner)
        // 2. Store set of adjacent charts
        // 3. Compute average coordinates of the next vertex along each boundary. Maybe also take into account 2nd vertices.
        // 4. Find the nearest vertex of these average coordinates. If still the same vertex, skip (continue with next corner)
        // 5. Store set of adjacent charts of the nearest vertex. If different charts, skip (continue with next corner), because the new corner position will break validity
        // 6. Adjust labeling around the vertex found, so that the next call of LabelingGraph::fill_from() will move the corner and neighboring boundary edges

        const Corner& current_corner = lg.corners[c];

        // 1.

        if(current_corner.all_adjacent_boundary_edges_are_on_feature_edges(mesh,feature_edges)) {
            fmt::println(Logger::out("monotonicity"),"Corner {} ignored because surrounded by feature edges",c); Logger::out("monotonicity").flush();
            continue; // do not move this corner
        }
        // TODO if some adjacent boundary edges are on feature edges but not all,
        // constraint the new position of the corner to be on a feature edge

        // 2.

        lg.get_adjacent_charts_of_vertex(current_corner.vertex,adj_facets,adjacent_charts_of_corner);

        // 3.

        average_coordinates = current_corner.average_coordinates_of_neighborhood(mesh,lg,false,1);

        // 4.

        nearest_vertex_of_average_coordinates = get_nearest_vertex_of_coordinates(mesh_he,adj_facets,average_coordinates,current_corner.vertex,1);
        if(nearest_vertex_of_average_coordinates == current_corner.vertex) {
            fmt::println(Logger::out("monotonicity"),"Corner {} ignored because new position is on the same vertex",c); Logger::out("monotonicity").flush();
            continue;
        }

        // 5.
        
        lg.get_adjacent_charts_of_vertex(nearest_vertex_of_average_coordinates,adj_facets,adjacent_charts_of_new_vertex);
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

bool merge_a_turning_point_and_its_closest_corner(GEO::Mesh& mesh, const char* attribute_name, const LabelingGraph& lg, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets) {
    
    if(lg.non_monotone_boundaries.empty()) {
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

    FOR(non_monotone_boundary_index,lg.non_monotone_boundaries.size()) { // for each non-monotone boundary
        const Boundary& non_monotone_boundary = lg.boundaries[lg.non_monotone_boundaries[non_monotone_boundary_index]];

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
        MeshHalfedgesExt mesh_he(mesh);
        mesh_he.set_use_facet_region(attribute_name);
        index_t first_turning_point_on_feature_edge_vertex = first_turning_point_on_feature_edge.vertex(non_monotone_boundary,mesh);

        // 3. Only process `first_turning_point_on_feature_edge`. Find which corner of `non_monotone_boundary` is the closest.

        #ifndef NDEBUG
            dump_vertex("current_tp",mesh_vertex(mesh,first_turning_point_on_feature_edge_vertex));
        #endif
        
        index_t closest_corner = first_turning_point_on_feature_edge.get_closest_corner(non_monotone_boundary,mesh_he);
        #ifndef NDEBUG
            dump_vertex("closest_corner",mesh,lg.corners[closest_corner].vertex);
        #endif

        // 4. Find which boundary around `closest_corner` is closest to `first_turning_point_on_feature_edge`, excluding `non_monotone_boundary`
        //    Also compute the vector between its corners

        const Boundary& boundary_to_move = lg.boundaries[non_monotone_boundary.get_closest_boundary_of_turning_point(first_turning_point_on_feature_edge,closest_corner,mesh_he,lg.halfedge2boundary,lg.corners)];
        #ifndef NDEBUG
            dump_boundary_with_halfedges_indices("boundary_to_move",mesh,boundary_to_move);
        #endif
        vec3 boundary_to_move_vector = boundary_to_move.vector_between_corners(mesh,lg.corners);
        // flip the vector if `boundary_to_move.end_corner` is adjacent to `non_monotone_boundary` and not `boundary_to_move.start_corner`
        if(lg.corners[boundary_to_move.start_corner].vertex != lg.corners[closest_corner].vertex) {
            geo_assert(lg.corners[boundary_to_move.end_corner].vertex == lg.corners[closest_corner].vertex);
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
        #ifndef NDEBUG
            const std::set<index_t>& facets_of_the_chart_on_which_the_new_boundary_will_be = lg.charts[chart_on_which_the_new_boundary_will_be].facets;
            dump_facets("facets_of_the_chart_on_which_the_new_boundary_will_be",mesh,facets_of_the_chart_on_which_the_new_boundary_will_be);
        #endif
        index_t new_label = index_t(-1);
        if(chart_on_which_the_new_boundary_will_be == boundary_to_move.right_chart) {
            new_label = lg.charts[boundary_to_move.left_chart].label;
        }
        else {
            new_label = lg.charts[boundary_to_move.right_chart].label;
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
            follow_feature_edge_on_chart(mesh_he,halfedge_on_lost_feature_edge,feature_edges,lg.facet2chart,facets_at_left,facets_at_right);
            const std::set<index_t>& facets_to_re_label = 
                (closest_corner == non_monotone_boundary.end_corner) ?
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_left : facets_at_right) : 
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_right : facets_at_left);
            const std::set<index_t>& wall =
                (closest_corner == non_monotone_boundary.end_corner) ?
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_right : facets_at_left) :
                (first_turning_point_on_feature_edge.is_towards_right_ ? facets_at_left : facets_at_right);
            propagate_label(mesh,attribute_name,new_label,facets_to_re_label,wall,lg.facet2chart,chart_on_which_the_new_boundary_will_be);
            return true;
        }

    default_behavior:

        trace_path_on_chart(mesh_he,adj_facets,lg.facet2chart,lg.turning_point_vertices,first_turning_point_on_feature_edge_vertex,boundary_to_move_vector,facets_at_left,facets_at_right,path);
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
                label[facet] = lg.charts[chart_on_which_the_new_boundary_will_be].label;
            }
        }

        propagate_label(mesh,attribute_name,new_label,facets_to_process,wall_facets,lg.facet2chart,chart_on_which_the_new_boundary_will_be);

        return true;

    }
    // so none of the non-monotone boundary can be processed
    return false;
}

// returns true if a chart has been created
// returns false if any of the turning-points can be processed with this operator
bool join_turning_points_pair_with_new_chart(GEO::Mesh& mesh, const char* attribute_name, LabelingGraph& lg, const std::vector<vec3>& normals, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<std::vector<index_t>>& adj_facets) {
    if(lg.non_monotone_boundaries.empty()) {
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
    MeshHalfedgesExt mesh_he(mesh);
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

    for(index_t b : lg.non_monotone_boundaries) { // for each non-monotone boundary (b is an index of lg.boundaries)
        for(const auto& tp : lg.boundaries[b].turning_points) { // for each turning-point of the current boundary
            if(vertex_has_lost_feature_edge_in_neighborhood(mesh_he,adj_facets,feature_edges,tp.vertex(lg.boundaries[b],mesh),halfedge)) {
                chart_on_which_the_lost_feature_edge_is = lg.facet2chart[halfedge_facet_left(mesh,halfedge)];
                geo_assert(chart_on_which_the_lost_feature_edge_is == lg.facet2chart[halfedge_facet_right(mesh,halfedge)]);
                label_on_which_the_lost_feature_edge_is = lg.charts[chart_on_which_the_lost_feature_edge_is].label;
                halfedge = follow_feature_edge_on_chart(mesh_he,halfedge,feature_edges,lg.facet2chart,facets_at_left,facets_at_right);
                // so move unsuccessful (no more halfedges on feature edge), or we left the chart (found a boundary / corner / turning-point)
                vertex_at_tip_of_feature_edge = halfedge_vertex_index_to(mesh,halfedge);
                if(lg.turning_point_vertices.contains(vertex_at_tip_of_feature_edge)) {
                    geo_assert(lg.turning_point_vertices[vertex_at_tip_of_feature_edge].size() == 1); // assert only one turning-point associated to this vertex
                    boundary_of_the_second_turning_point = lg.non_monotone_boundaries[lg.turning_point_vertices[vertex_at_tip_of_feature_edge][0].first];
                    new_label = find_optimal_label(
                        {},
                        { // forbidden labels
                            label_on_which_the_lost_feature_edge_is,
                            lg.boundaries[b].other_label(lg.charts,label_on_which_the_lost_feature_edge_is),
                            lg.boundaries[boundary_of_the_second_turning_point].other_label(lg.charts,label_on_which_the_lost_feature_edge_is)
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
                            if(lg.facet2chart[adjacent_facet] == chart_on_which_the_lost_feature_edge_is) {
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

bool auto_fix_monotonicity(Mesh& mesh, const char* attribute_name, LabelingGraph& lg, std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<vec3>& normals) {
    
    size_t nb_processed = 0;

    if (lg.non_monotone_boundaries.empty()) {
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
        nb_processed = (size_t) join_turning_points_pair_with_new_chart(mesh,attribute_name,lg,normals,feature_edges,adj_facets);
        lg.fill_from(mesh,attribute_name,feature_edges);
    } while (nb_processed != 0);

    if (lg.non_monotone_boundaries.empty()) {
        return true;
    }

    size_t nb_iter = 0;
    do {
        nb_processed = (size_t) merge_a_turning_point_and_its_closest_corner(mesh,attribute_name,lg,feature_edges,adj_facets);
        lg.fill_from(mesh,attribute_name,feature_edges);
        nb_iter++;
        if(nb_iter > 20) {
            break;
        }
    } while (nb_processed != 0);

    if (lg.non_monotone_boundaries.empty()) {
        return true;
    }

    move_boundaries_near_turning_points(mesh,attribute_name,lg,feature_edges);
    lg.fill_from(mesh,attribute_name,feature_edges);

    if (lg.non_monotone_boundaries.empty()) {
        return true;
    }

    straighten_boundaries(mesh,attribute_name,lg,adj_facets,feature_edges);
    // `lg` already updated in straighten_boundaries()

    if (lg.non_monotone_boundaries.empty()) {
        return true;
    }

    fmt::println(Logger::out("monotonicity"),"auto fix monotonicity stopped didn't reach all-monotone boundaries"); Logger::out("monotonicity").flush();
    return false;
}