#include <geogram/basic/attributes.h>   // for GEO::Attribute

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/os.h>

#include <array>            // for std::array
#include <initializer_list> // for std::initializer_list
#include <algorithm>        // for std::max_element(), std::min(), std::max()
#include <iterator>         // for std::distance()
#include <tuple>            // for std::tuple

#include "labeling.h"
#include "LabelingGraph.h"
#include "containers.h"
#include "GraphCutLabeling.h"
#include "CustomMeshHalfedges.h"

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

std::tuple<double,double,double> compute_per_facet_fidelity(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* labeling_attribute_name, const char* fidelity_attribute_name) {
    // goal of the fidelity metric : measure how far the assigned direction (label) is from the triangle normal
    //   fidelity=1 (high fidelity) -> label equal to the normal      ex: triangle is oriented towards +X and we assign +X
    //   fidelity=0 (low fidelity)  -> label opposite to the normal   ex: triangle is oriented towards +X and we assign -X
    // 1. compute and normalize the normal
    // 2. compute the vector (ex: 0.0,1.0,0.0) of the label (ex: +Y) with label2vector[] (already normalized)
    // 3. compute the dot product, which is in [-1:1] : 1=equal, -1=opposite
    // 4. add 1 and divide by 2 so that the range is [0:1] : 1=equal, 0=opposite
    Attribute<index_t> label(mesh.facets.attributes(), labeling_attribute_name);
    Attribute<double> per_facet_fidelity(mesh.facets.attributes(), fidelity_attribute_name);
    double min = Numeric::max_float64(),
           max = Numeric::min_float64(),
           avg = 0.0;
    FOR(f,mesh.facets.nb()) {
        per_facet_fidelity[f] = (GEO::dot(normals[f],label2vector[label[f]]) + 1.0)/2.0;
        // stats
        min = std::min(min,per_facet_fidelity[f]);
        max = std::max(max,per_facet_fidelity[f]);
        avg += per_facet_fidelity[f];
    }
    avg /= mesh.facets.nb();
    return {min,max,avg};
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

unsigned int move_boundaries_near_turning_points(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {

    unsigned int nb_labels_changed = 0;

    if(slg.non_monotone_boundaries.empty()) {
        fmt::println(Logger::out("fix_labeling"),"Warning : operation canceled because all boundaries are monotone"); Logger::out("fix_labeling").flush();
        return nb_labels_changed;
    }

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    CustomMeshHalfedges mesh_halfedges(mesh); // create an halfedges interface for this mesh
    mesh_halfedges.set_use_facet_region(attribute_name);
    MeshHalfedges::Halfedge current_halfedge, previous_halfedge;

    // Get turning points, change the label if this doesnt break charts connectivity

    for(auto b : slg.non_monotone_boundaries) {
        for(auto v : slg.boundaries[b].turning_points) { // v is the vertex at the base of halfedge v

            current_halfedge = slg.boundaries[b].halfedges[v];
            geo_assert(mesh_halfedges.halfedge_is_border(current_halfedge));
            geo_assert(MAP_CONTAINS(slg.halfedge2boundary,current_halfedge));

            // get the vertex index
            index_t current_vertex = mesh.facet_corners.vertex(current_halfedge.corner);
            geo_assert(slg.is_turning_point(mesh,current_vertex));
            geo_assert(slg.vertex2corner[current_vertex] == index_t(-1)); // a turning point should not be a corner
            // test if the valence of current_vertex is 2
            VertexRingWithBoundaries vr;
            vr.explore(current_halfedge,mesh_halfedges);
            geo_assert(vr.valence() == 2); // should have only 2 charts

            // Explore clockwise : right side then left side
            // For each side, sum of angles & aggregate facet indices,
            // then replace labels of smallest side (angle-wise) with label of other side
            double left_side_sum_of_angles = 0.0;
            std::vector<index_t> left_side_facets;
            double right_side_sum_of_angles = 0.0;
            std::vector<index_t> right_side_facets;

            previous_halfedge = current_halfedge;
            mesh_halfedges.move_to_next_around_vertex(previous_halfedge,true); // next counterclockwise
            // explore the right side (clockwise)
            do {
                right_side_facets.push_back(current_halfedge.facet);
                right_side_sum_of_angles += angle(halfedge_vector(mesh,previous_halfedge),halfedge_vector(mesh,current_halfedge));
                previous_halfedge = current_halfedge;
                mesh_halfedges.move_to_prev_around_vertex(current_halfedge,true); // previous counterclockwise
            } while (!mesh_halfedges.halfedge_is_border(current_halfedge));
            // explore the right side
            do {
                left_side_facets.push_back(current_halfedge.facet);
                left_side_sum_of_angles += angle(halfedge_vector(mesh,previous_halfedge),halfedge_vector(mesh,current_halfedge));
                previous_halfedge = current_halfedge;
                mesh_halfedges.move_to_prev_around_vertex(current_halfedge,true); // previous counterclockwise
            } while (!mesh_halfedges.halfedge_is_border(current_halfedge));
            
            if(left_side_sum_of_angles < right_side_sum_of_angles) {
                // replace label of left side facet with right side label
                for(auto f : left_side_facets) {
                    label[f] = slg.charts[slg.boundaries[b].left_chart].label; // why label of left_chart and not the one of right_chart ?
                }
                nb_labels_changed += (unsigned int) left_side_facets.size();
            }
            else {
                // replace label of right side facet with left side label
                for(auto f : right_side_facets) {
                    label[f] = slg.charts[slg.boundaries[b].right_chart].label; // why label of right_chart and not the one of left_chart ?
                }
                nb_labels_changed += (unsigned int) left_side_facets.size();
            }

        }
    }
    return nb_labels_changed;
}

void straighten_boundary(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg, index_t boundary_index) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute
    index_t left_chart_index = slg.boundaries[boundary_index].left_chart;
    index_t right_chart_index = slg.boundaries[boundary_index].right_chart;
    const Chart& left_chart = slg.charts[left_chart_index];
    const Chart& right_chart = slg.charts[right_chart_index];
    auto gcl = GraphCutLabeling(mesh,normals,left_chart.facets.size()+right_chart.facets.size()); // graph-cut only on the two charts
    for(index_t f : left_chart.facets)
        gcl.add_site(f);
    for(index_t f : right_chart.facets)
        gcl.add_site(f);
    gcl.data_cost__set__fidelity_based(1);
    // lock label of triangles on the 2 charts boundaries intersection
    index_t chart_of_adjacent_facet = index_t(-1);
    for(index_t f : left_chart.facets) {
        FOR(le,3) {
            chart_of_adjacent_facet = slg.facet2chart[mesh.facets.adjacent(f,le)];
            if( (chart_of_adjacent_facet != left_chart_index) && 
                (chart_of_adjacent_facet != right_chart_index) ) {
                // the facet f is next to a chart different from the 2 next to the boundary
                gcl.data_cost__change_to__locked_label(f,label[f]);
            }
        }
    }
    for(index_t f : right_chart.facets) {
        FOR(le,3) {
            chart_of_adjacent_facet = slg.facet2chart[mesh.facets.adjacent(f,le)];
            if( (chart_of_adjacent_facet != left_chart_index) && 
                (chart_of_adjacent_facet != right_chart_index) ) {
                // the facet f is next to a chart different from the 2 next to the boundary
                gcl.data_cost__change_to__locked_label(f,label[f]);
            }
        }
    }
    gcl.smooth_cost__set__default();
    gcl.neighbors__set__compactness_based(3);
    gcl.compute_solution(label);
}