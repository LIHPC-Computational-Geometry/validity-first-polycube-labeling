#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/mesh/mesh_io.h>

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/os.h>

#include <array>            // for std::array
#include <initializer_list> // for std::initializer_list
#include <algorithm>        // for std::max_element(), std::min(), std::max()
#include <iterator>         // for std::distance()
#include <tuple>            // for std::tuple
#include <queue>            // for std::queue
#include <cmath>            // for std::pow(), std::abs()
#include <set>              // for std::set

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
    return (index_t) std::distance(weights.begin(),std::max_element(weights.begin(),weights.end()));
}

void naive_labeling(Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name) {

    // use GEO::Geom::triangle_normal_axis() instead ?

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    for(index_t f: mesh.facets) { // for each facet
        label[f] = nearest_label(normals[f]);
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
    // goal of the fidelity metric : measure how far the assigned direction (label) is from the triangle normal
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

unsigned int remove_surrounded_charts(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);

    // Get charts surronded by only 1 label
    // Broader fix than just looking at charts having 1 boundary

    unsigned int modified_charts_count = 0;
    for(index_t chart_index : slg.invalid_charts) {
        auto boundary_iterator = slg.charts[chart_index].boundaries.begin();
        index_t chart_at_other_side = slg.boundaries[*boundary_iterator].chart_at_other_side(chart_index);
        index_t surrounding_label = slg.charts[chart_at_other_side].label;

        boundary_iterator++; // go to next boundary
        for(;boundary_iterator != slg.charts[chart_index].boundaries.end();++boundary_iterator) {
            chart_at_other_side = slg.boundaries[*boundary_iterator].chart_at_other_side(chart_index);
            if(surrounding_label != slg.charts[chart_at_other_side].label) {
                goto skip_modification; // this chart has several labels around (goto because break in nested loops is a worse idea)
            }
        }

        // if we are here, the goto was not used, so the only label around is surrounding_label
        for(index_t facet_index : slg.charts[chart_index].facets) {
            label[facet_index] = surrounding_label;
        }
        modified_charts_count++;


        skip_modification:
            ; // just end this loop interation
    }

    return modified_charts_count;
}

unsigned int fix_invalid_boundaries(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_half_edges_(mesh); // create an halfedges interface for this mesh

    // For each invalid boundary,
    // go through each vertex
    // and modify the labels in the vertex ring

    unsigned int new_charts_count = 0;
    index_t new_label;
    MeshHalfedges::Halfedge current_halfedge;
    for(index_t boundary_index : slg.invalid_boundaries) { // for each invalid boundary
        const Boundary& current_boundary = slg.boundaries[boundary_index];
        new_label = nearest_label(current_boundary.average_normal);

        // Change the labels around the start/end corner, but only in place of the 2 charts next to the boundary
        for(auto current_corner : std::initializer_list<index_t>({current_boundary.start_corner,current_boundary.end_corner})) {
            geo_assert(current_corner != index_t(-1));
            for(auto& vr : slg.corners[current_corner].vertex_rings_with_boundaries) { // for each vertex ring
                for(auto& be : vr.boundary_edges) { // for each boundary edge
                    if(
                    (slg.facet2chart[be.facet] == current_boundary.left_chart) || 
                    (slg.facet2chart[be.facet] == current_boundary.right_chart) ) {
                        // so be.facet is on one of the charts to shrink
                        label[be.facet] = new_label;
                    }
                }
            }
        }

        // Change the labels around the vertices between the two corners
        auto boundary_halfedge = current_boundary.halfedges.cbegin();
        boundary_halfedge++;
        for(; boundary_halfedge != current_boundary.halfedges.cend(); ++boundary_halfedge) { // walk along the boundary
            current_halfedge = (*boundary_halfedge); // copy the boundary halfedge into a mutable variable
            // modify the labels around this vertex
            do {
                label[current_halfedge.facet] = new_label;
                mesh_half_edges_.move_to_next_around_vertex(current_halfedge,true);
            } while (current_halfedge != *boundary_halfedge); // go around the vertex, until we are back on the initial boundary edge
        }

        new_charts_count++;
    }

    return new_charts_count; // should be == to slg.invalid_boundaries.size()
}

unsigned int fix_invalid_corners(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_half_edges_(mesh); // create an halfedges interface for this mesh

    // Replace the labels around invalid corners
    // by the nearest one from the vertex normal

    unsigned int new_charts_count = 0;
    index_t new_label;
    vec3 vertex_normal = {0,0,0};
    for(index_t corner_index : slg.invalid_corners) { // for each invalid corner

        // compute the normal by adding normals of adjacent facets
        vertex_normal = {0,0,0};
        for(auto& vr : slg.corners[corner_index].vertex_rings_with_boundaries) { // for each vertex ring
            for(auto& be : vr.boundary_edges) { // for each boundary edge
                vertex_normal += normals[be.facet]; // compute normal of associated facet, update vertex_normal
            }
        }

        // compute new label
        new_label = nearest_label(vertex_normal);

        // change the label of adjacent facets
        for(auto& vr : slg.corners[corner_index].vertex_rings_with_boundaries) { // for each vertex ring
            for(auto& be : vr.boundary_edges) { // for each boundary edge
                label[be.facet] = new_label; // change the label of this facet
            }
        }

        new_charts_count++;
    }

    return new_charts_count; // should be == to slg.invalid_corners.size()
}

// from https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/labeling_ops.cpp removeChartMutation()
void remove_invalid_charts(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg) {

    if(slg.invalid_charts.empty()) {
        fmt::println(Logger::out("fix_labeling"),"Warning : operation canceled because there are no invalid charts"); Logger::out("fix_labeling").flush();
        return;
    }

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute

    // Fill invalid charts using a Graph-Cut optimization,
    // preventing existing label from being re-applied

    // compactness = 1
    GraphCutLabeling gcl(mesh,normals);
    gcl.data_cost__set__locked_labels(label); // start by locking all the labels, so valid charts will not be modified
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(1);

    for(index_t chart_index : slg.invalid_charts) { // for each invalid chart

        for(index_t facet_index : slg.charts[chart_index].facets) { // for each facet inside this chart
            gcl.data_cost__change_to__fidelity_based(facet_index,1);
            // if facet next to a boundary, lower the cost of assigning the neighboring label
            FOR(le,3) { // for each local edge
                index_t adjacent_facet = mesh.facets.adjacent(facet_index,le);
                if(label[adjacent_facet] != label[facet_index]) {
                    gcl.data_cost__change_to__scaled(facet_index,label[adjacent_facet],0.5f); // halve the cost
                }
            }
            gcl.data_cost__change_to__forbidden_label(facet_index,label[facet_index]); // prevent the label from staying the same
        }
    }

    gcl.compute_solution(label);
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
            gcl.data_cost__change_to__forbidden_label(facet_index,label[facet_index]); // prevent the label from staying the same
        }
    }

    gcl.compute_solution(label);

}

unsigned int move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {

    unsigned int nb_labels_changed = 0;

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
        for(auto tp : slg.boundaries[b].turning_points) {

            initial_halfedge = slg.boundaries[b].halfedges[tp.outgoing_local_halfedge_index()];
            geo_assert(mesh_halfedges.halfedge_is_border(initial_halfedge));
            geo_assert(MAP_CONTAINS(slg.halfedge2boundary,initial_halfedge));

            // get the vertex index
            index_t current_vertex = mesh.facet_corners.vertex(initial_halfedge.corner);
            geo_assert(slg.is_turning_point(mesh,current_vertex));
            geo_assert(slg.vertex2corner[current_vertex] == index_t(-1)); // a turning point should not be a corner
            // test if the valence of current_vertex is 2
            VertexRingWithBoundaries vr;
            vr.explore(initial_halfedge,mesh_halfedges);
            geo_assert(vr.valence() == 2); // should have only 2 charts

            // new label according to towards which chart the turning point is
            new_label = slg.charts[tp.towards_left() ? slg.boundaries[b].left_chart : slg.boundaries[b].right_chart].label;

            current_halfedge = initial_halfedge;
            // go around the vertex and assign new_label to adjacent facets
            do {
                mesh_halfedges.move_to_next_around_vertex(current_halfedge,true);
                if(label[current_halfedge.facet] != new_label) {
                    label[current_halfedge.facet] = new_label;
                    nb_labels_changed++;
                }
            } while (current_halfedge != initial_halfedge);
        }
    }
    return nb_labels_changed;
}

void straighten_boundary(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg, index_t boundary_index) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    const Boundary& current_boundary = slg.boundaries[boundary_index];
    index_t left_chart_index = current_boundary.left_chart;
    index_t right_chart_index = current_boundary.right_chart;
    const Chart& left_chart = slg.charts[left_chart_index];
    const Chart& right_chart = slg.charts[right_chart_index];
    index_t chart_of_adjacent_facet = index_t(-1);
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
        dump_all_boundaries("boundaries",mesh,mesh_he,slg.boundaries);

        std::map<index_t,bool> facets_of_the_2_charts;
        for(auto f : left_chart.facets) {
            facets_of_the_2_charts[f] = 0; // insert this facet, associated to 0=left
        }
        for(auto f : right_chart.facets) {
            facets_of_the_2_charts[f] = 1; // insert this facet, associated to 1=right
        }
        dump_facets("2_charts","on_wich_side",mesh,facets_of_the_2_charts);
    #endif

    // https://stackoverflow.com/a/72437022
    // Compute the distance to the boundary for all facets in the left and right charts
    std::map<index_t,unsigned int> distance_to_boundary; // map a facet index to a distance. only defined for facet inside the 2 charts
    std::queue<index_t> facets_to_explore; // facets we have to visit to update the distance in their neighborhood. FIFO data structure
    index_t adjacent_facet = index_t(-1);
    // Go through the boundary and set distance of adjacent triangle to 0
    auto boundary_halfedge = current_boundary.halfedges.cbegin(); // get an iterator pointing at the first halfedge
    boundary_halfedge++; // go to the second halfedge (the vertex at the base of the first one has facets of other charts next to it)
    CustomMeshHalfedges::Halfedge current_halfedge;
    for(; boundary_halfedge != current_boundary.halfedges.cend(); ++boundary_halfedge) { // walk along the boundary, from (boundary) halfedge to (boundary) halfedge
        current_halfedge = (*boundary_halfedge); // copy the boundary halfedge into a mutable variable
        do { // for each facet around the vertex at the base of the current boundary halfedge
            distance_to_boundary[current_halfedge.facet] = 0; // set distance to 0 (the current facet touch the boundary by an edge or a vertex)
            facets_to_explore.emplace(current_halfedge.facet);
            mesh_he.move_to_next_around_vertex(current_halfedge,true);
        } while (current_halfedge != *boundary_halfedge); // go around the vertex, until we are back on the initial boundary edge
    }

    #ifndef NDEBUG
        fmt::println("boundary expored, distance of adjacent facet set to 0");
    #endif

    index_t current_facet = index_t(-1);
    unsigned int current_distance = 0; // will store the distance currently computed for the current facet
    while(!facets_to_explore.empty()) { // while there still are facets to explore
        current_facet = facets_to_explore.front(); // pick the facet at the front of the FIFO
        facets_to_explore.pop(); // remove it from the FIFO
        current_distance = distance_to_boundary[current_facet];
        FOR(le,3) { // for each local edge of the current facet
            adjacent_facet = mesh.facets.adjacent(current_facet,le); // get the facet index of the neighbor at the other side of the current local edge
            chart_of_adjacent_facet = slg.facet2chart[adjacent_facet]; // get its chart
            if( (chart_of_adjacent_facet != left_chart_index) && 
                (chart_of_adjacent_facet != right_chart_index) ) {
                    continue; // -> do not explore the adjacent facet, it is not inside the 2 charts
                }
            if(!distance_to_boundary.contains(adjacent_facet)) { // if the adjacent facet is not yet linked to a distance
                distance_to_boundary[adjacent_facet] = current_distance + 1; // set an initial distance of one more than the current facet distance
                facets_to_explore.emplace(adjacent_facet);
                continue;
            }
            if (current_distance + 1 < distance_to_boundary[adjacent_facet]) { // if the distance of the neighbor was too much (passing by the current facet is closer)
                distance_to_boundary[adjacent_facet] = current_distance + 1; // update the distance
                facets_to_explore.emplace(adjacent_facet);
            }
        }
    }

    #ifndef NDEBUG
        fmt::println("distance_to_boundary has {} elements",distance_to_boundary.size());
        fmt::println("combined, the 2 charts have {} facets",left_chart.facets.size()+right_chart.facets.size());
    #endif

    geo_assert(distance_to_boundary.size() == left_chart.facets.size()+right_chart.facets.size());

    #ifndef NDEBUG
        // Export the per-facet distance in a .geogram file
        Mesh per_facet_dist;
        per_facet_dist.copy(mesh,false);
        // keep vertices and facets
        per_facet_dist.edges.clear();
        per_facet_dist.cells.clear();
        Attribute<float> dist(per_facet_dist.facets.attributes(),"dist");
        dist.fill(-1.0f); // set distance of -1 on the whole surface
        for(auto& kv : distance_to_boundary) { // for each facet in the 2 charts
            dist[kv.first] = (float) kv.second; // overwrite with the computed distance
        }
        mesh_save(per_facet_dist,"per_facet_dist.geogram");

        // Export the contour (facets in the perimeter of the union of the 2 charts)
        Mesh contour;
        contour.copy(mesh,false);
        // keep vertices and facets
        contour.edges.clear();
        contour.cells.clear();
        Attribute<bool> on_contour(contour.facets.attributes(),"on_contour");
        dist.fill(false); // init value: not on contour
        for(const auto& kv : distance_to_boundary) { // for each facet in the 2 charts
            if ( !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,0)) || 
                 !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,1)) || 
                 !distance_to_boundary.contains(mesh.facets.adjacent(kv.first,2)) ) { // if one of the neighbors is not inside the 2 charts
                on_contour[kv.first] = true; // on the contour
            }
        }
        mesh_save(contour,"contour.geogram");
    #endif

    auto gcl = GraphCutLabeling(mesh,normals,distance_to_boundary.size()); // graph-cut only on the two charts
    for(const auto& kv : distance_to_boundary) {
        gcl.add_site(kv.first);
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
            gcl.data_cost__change_to__locked_label(kv.first,label[kv.first]);
        }
        else {
            // penalize modification proportionally to the distance to the boundary (only close facets should be modifiable)
            // kv.second is the distance
            // -> multiply the cost of re-assigning by (0.8)^distance
            gcl.data_cost__change_to__scaled(kv.first,label[kv.first],std::pow(0.8,kv.second));
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
                    gcl.neighbors__change_to__shifted(kv.first,adjacent_facet,(1-dot_product)*500); // add a penalty the more ⟂ the edge is (// -> 0, ⟂ -> 300)
                    already_processed.insert({kv.first, adjacent_facet});
                    #ifndef NDEBUG
                        per_edge_dot_product[std::make_pair(vertex1,vertex2)] = dot_product;
                    #endif
                }
            }
        }
        #ifndef NDEBUG
            dump_edges("dot_products","dot_product",mesh,mesh_he,per_edge_dot_product);
        #endif
    }
    gcl.compute_solution(label);
}

void pull_closest_corner(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg, index_t non_monotone_boundary_index) {
    const Boundary& boundary = slg.boundaries[slg.non_monotone_boundaries[non_monotone_boundary_index]];
    if(boundary.turning_points.size() != 1) {
        fmt::println(Logger::out("monotonicity"),"Ignoring boundary {} which has {} turning-points instead of 1 for pull_closest_corner()",slg.non_monotone_boundaries[non_monotone_boundary_index],boundary.turning_points.size()); Logger::out("monotonicity").flush();
        return;
    }
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    const TurningPoint& tp = boundary.turning_points[0];
    index_t halfedge_index_turning_point = tp.outgoing_local_halfedge_index();
    #ifndef NDEBUG
        dump_vertex("current_tp",mesh,boundary.turning_point_vertex(0,mesh));
    #endif
    index_t closest_corner = index_t(-1);
    index_t new_label = index_t(-1); // label to assign between the boundary to remove, and its "parallel" passing by the turning point
    CustomMeshHalfedges mesh_he(mesh);
    mesh_he.set_use_facet_region(attribute_name);
    MeshHalfedges::Halfedge current_halfedge; // first, will store outgoing halfedge of the closest_corner
    // the turning-point is at the base of the halfedge 'halfedge_index_turning_point'
    // compute 2 distances :
    // - from start corner to turning point (= first halfedge to halfedge_index_turning_point-1 included)
    // - from turning point to end corner (= halfedge_index_turning_point to last halfedge included)
    double distance_to_start_corner = 0.0;
    double distance_to_end_corner = 0.0;
    FOR(i,halfedge_index_turning_point-1) {
        distance_to_start_corner += length(Geom::halfedge_vector(mesh,boundary.halfedges[i]));
    }
    FOR(i,halfedge_index_turning_point-1) {
        distance_to_end_corner += length(Geom::halfedge_vector(mesh,boundary.halfedges[i]));
    }
    // TODO check association between closest corner, turning point direction and direction around vertex (clockwise/counterclockwise)
    if(distance_to_start_corner < distance_to_end_corner) {
        closest_corner = boundary.start_corner;
        current_halfedge = boundary.halfedges[0];
        if(tp.towards_left()) {
            do {
                mesh_he.move_to_prev_around_vertex(current_halfedge,true); // next clockwise
            }
            while(!mesh_he.halfedge_is_border(current_halfedge));
            // do while could be replace by move_to_prev_around_vertex(current_halfedge,false); ?
        }
        else {
            do {
                mesh_he.move_to_next_around_vertex(current_halfedge,true); // next counterclockwise
            }
            while(!mesh_he.halfedge_is_border(current_halfedge));
            // do while could be replace by move_to_next_around_vertex(current_halfedge,false); ?
        }
    }
    else {
        closest_corner = boundary.end_corner;
        current_halfedge = *boundary.halfedges.rbegin();
        mesh_he.move_to_opposite(current_halfedge);
        if(tp.towards_right()) {
            do {
                mesh_he.move_to_prev_around_vertex(current_halfedge,true); // next clockwise
            }
            while(!mesh_he.halfedge_is_border(current_halfedge));
            // do while could be replace by move_to_prev_around_vertex(current_halfedge,false); ?
        }
        else {
            do {
                mesh_he.move_to_next_around_vertex(current_halfedge,true); // next counterclockwise
            }
            while(!mesh_he.halfedge_is_border(current_halfedge));
            // do while could be replace by move_to_next_around_vertex(current_halfedge,false); ?
        }
    }
    #ifndef NDEBUG
        dump_vertex("closest_corner",mesh,slg.corners[closest_corner].vertex);
    #endif
    const Boundary& boundary_to_move = slg.boundaries[slg.halfedge2boundary.at(current_halfedge).first];
    #ifndef NDEBUG
        Mesh boundary_to_move_file;
        boundary_to_move_file.copy(mesh,false,MESH_VERTICES);
        boundary_to_move_file.edges.create_edges(boundary_to_move.halfedges.size());
        FOR(he_index,boundary_to_move.halfedges.size()) {
            MeshHalfedges::Halfedge tmp_he = boundary_to_move.halfedges[he_index];
            boundary_to_move_file.edges.set_vertex(he_index,0,mesh.facet_corners.vertex(tmp_he.corner));
            mesh_he.move_to_opposite(tmp_he);
            boundary_to_move_file.edges.set_vertex(he_index,1,mesh.facet_corners.vertex(tmp_he.corner));
        }
        mesh_save(boundary_to_move_file,"boundary_to_move_file.geogram");
    #endif
    // TODO check new_label computation
    if(tp.towards_right()) {
        new_label = slg.charts[boundary_to_move.right_chart].label;
    }
    else {
        new_label = slg.charts[boundary_to_move.left_chart].label;
    }
    // prepare a graph-cut optimization on both sides (2 charts)
    // TODO distance-based : the further a facet is, the bigger is the cost of changing its label. (too restrictive for this operator?)
    // TODO restrict possible labels to the ones of the 2 charts
    std::set<index_t> union_of_2_charts;
    union_of_2_charts.insert(slg.charts[boundary_to_move.left_chart].facets.begin(),slg.charts[boundary_to_move.left_chart].facets.end());
        union_of_2_charts.insert(slg.charts[boundary_to_move.right_chart].facets.begin(),slg.charts[boundary_to_move.right_chart].facets.end());
    #ifndef NDEBUG
        dump_facets("2_charts",mesh,union_of_2_charts);
    #endif
    auto gcl = GraphCutLabeling(mesh,normals,union_of_2_charts.size()); // graph-cut only on the two charts
    for(index_t f : union_of_2_charts) {
        gcl.add_site(f);
    }
    gcl.data_cost__set__fidelity_based(1);
    // Lock labels between the turning point and the corner, to ensure the corner moves where the turning point is
    // Note that MeshHalfedges::Halfedge::facet is the facet at the LEFT of the halfedge
    #ifndef NDEBUG
        std::set<index_t> facets_with_locked_label;
    #endif
    if(closest_corner == boundary.start_corner) {
        // go across the boundary in the opposite way
        for(signed_index_t halfedge_index = halfedge_index_turning_point-1; halfedge_index>= 0; halfedge_index--) {
            current_halfedge = boundary.halfedges[halfedge_index];
            if(tp.towards_right()) {
                mesh_he.move_to_opposite(current_halfedge);
            }
            gcl.data_cost__change_to__locked_label(current_halfedge.facet,new_label);
            #ifndef NDEBUG
                facets_with_locked_label.insert(current_halfedge.facet);
            #endif
        }
    }
    else {
        // go across the boundary in the same direction as the halfedges
        for(index_t halfedge_index = halfedge_index_turning_point; halfedge_index < boundary.halfedges.size(); halfedge_index++) {
            current_halfedge = boundary.halfedges[halfedge_index];
            if(tp.towards_right()) {
                mesh_he.move_to_opposite(current_halfedge);
            }
            gcl.data_cost__change_to__locked_label(current_halfedge.facet,new_label);
            #ifndef NDEBUG
                facets_with_locked_label.insert(current_halfedge.facet);
            #endif
        }
    }
    #ifndef NDEBUG
        dump_facets("facets_with_locked_label",mesh,facets_with_locked_label);
    #endif
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(5); // increase compactness/fidelity ratio. fidelity is less important for this operator
    index_t adjacent_facet = index_t(-1);
    if(boundary_to_move.axis != -1) {
        // penalize boundary tracing through halfedges poorly aligned with the boundary axis (X, Y or Z)
        // parse all undirected edge inside the 2 charts surface and tweak the cost
        #ifndef NDEBUG
            std::map<std::pair<index_t,index_t>,double> per_edge_dot_product;
        #endif
        index_t vertex1 = index_t(-1);
        index_t vertex2 = index_t(-1);
        vec3 edge_vector;
        double dot_product = 0.0; // dot product of the edge vector & boundary axis
        std::set<std::set<index_t>> already_processed; // use a set and not a pair so that the facets are sorted -> unordered pairs
        for(index_t f : union_of_2_charts) {
            FOR(le,3) { // for each local edge of the current facet
                adjacent_facet = mesh.facets.adjacent(f,le);
                if(already_processed.contains({f,adjacent_facet})) {
                    continue; // we already processed this edge (indices flipped)
                }
                if(union_of_2_charts.contains(adjacent_facet)) { // the edge at le is between 2 facets of the 2 charts union
                    // the vertices at the end point of the local edge le are le and (le+1)%3
                    // https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
                    vertex1 = mesh.facet_corners.vertex(mesh.facets.corner(f,le));
                    vertex2 = mesh.facet_corners.vertex(mesh.facets.corner(f,(le+1)%3));
                    edge_vector = normalize(mesh.vertices.point(vertex2)-mesh.vertices.point(vertex1));
                    dot_product = std::abs(dot(edge_vector,label2vector[boundary_to_move.axis*2])); // =1 -> //, =0 -> ⟂
                    gcl.neighbors__change_to__shifted(f,adjacent_facet,(1-dot_product)*500); // add a penalty the more ⟂ the edge is (// -> 0, ⟂ -> 300)
                    already_processed.insert({f, adjacent_facet});
                    #ifndef NDEBUG
                        per_edge_dot_product[std::make_pair(vertex1,vertex2)] = dot_product;
                    #endif
                }
            }
        }
        #ifndef NDEBUG
            dump_edges("dot_products","dot_product",mesh,mesh_he,per_edge_dot_product);
        #endif
    }
    gcl.compute_solution(label);
}