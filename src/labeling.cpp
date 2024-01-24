#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/numeric.h>      // for min_float64()
#include <geogram/basic/vecg.h>     // for vec3
#include <geogram/basic/matrix.h>   // for mat3

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/os.h>

#include <array>            // for std::array
#include <initializer_list> // for std::initializer_list
#include <algorithm>        // for std::max_element(), std::min(), std::max(), std::sort
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
    return (label1 % 2) != (label2 % 2); // orthogonal if they don't have the same axis
}

void naive_labeling(Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name) {

    // use GEO::Geom::triangle_normal_axis() instead ?

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    for(index_t f: mesh.facets) { // for each facet
        label[f] = nearest_label(normals[f]);
    }
}

void tweaked_naive_labeling(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name) {

    // https://en.wikipedia.org/wiki/Rotation_matrix#Basic_3D_rotations
    // rotation of NAIVE_LABELING_TWEAK_ANGLE around the x, y and z axes
    mat3 rotation; 
    rotation(0,0) = COS_SQUARED_TILT_ANGLE; // <=> cos*cos
    rotation(0,1) = SIN_BY_COS_TILT_ANGLE*(SIN_TILT_ANGLE-1); // <=> sin*sin*cos-cos*sin
    rotation(0,2) = SIN_TILT_ANGLE*(COS_SQUARED_TILT_ANGLE+SIN_TILT_ANGLE); // <=> cos*sin*cos+sin*sin
    rotation(1,0) = SIN_BY_COS_TILT_ANGLE; // <=> cos*sin
    rotation(1,1) = SIN_SQUARED_TILT_ANGLE*SIN_TILT_ANGLE+COS_SQUARED_TILT_ANGLE; // <=> sin*sin*sin+cos*cos
    rotation(1,2) = SIN_BY_COS_TILT_ANGLE*(SIN_TILT_ANGLE-1); // <=> cos*sin*sin-sin*cos
    rotation(2,0) = -SIN_TILT_ANGLE; // <=> -sin
    rotation(2,1) = SIN_BY_COS_TILT_ANGLE; // <=> sin*cos
    rotation(2,2) = COS_SQUARED_TILT_ANGLE; // <=> cos*cos

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    FOR(f,mesh.facets.nb()) { // for each facet
        vec3 normal = normals[f];
        // based on nearest_label()
        std::array<std::pair<double,index_t>,6> per_label_weights = {
            std::make_pair(normal.x < 0.0 ? 0.0 : normal.x,     0), // +X
            std::make_pair(normal.x < 0.0 ? -normal.x : 0.0,    1), // -X
            std::make_pair(normal.y < 0.0 ? 0.0 : normal.y,     2), // +Y
            std::make_pair(normal.y < 0.0 ? -normal.y : 0.0,    3), // -Y
            std::make_pair(normal.z < 0.0 ? 0.0 : normal.z,     4), // +Z
            std::make_pair(normal.z < 0.0 ? -normal.z : 0.0,    5)  // -Z
        };
        std::sort(per_label_weights.begin(),per_label_weights.end());
        if(std::abs(per_label_weights[5].first-per_label_weights[4].first) < NAIVE_LABELING_TWEAK_SENSITIVITY) {
            // the 2 labels with the most weight are too close
            // slightly rotate the the normal to go out of this indecisiveness area
            normal = mult(rotation,normal);

            // find the nearest label of the rotated normal
            label[f] = nearest_label(normal);
        }
        else {
            // assign the closest label
            label[f] = per_label_weights[5].second;
        }
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
                mesh_half_edges_.move_counterclockwise_around_vertex(current_halfedge,true);
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

void pull_closest_corner(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, const StaticLabelingGraph& slg, index_t non_monotone_boundary_index) {
    
    // get the boundary from its index

    const Boundary& boundary = slg.boundaries[slg.non_monotone_boundaries[non_monotone_boundary_index]];
    #ifndef NDEBUG
        dump_boundary_with_halfedges_indices("non-monotone_boundary",mesh,boundary);
    #endif

    // get its only turning point

    if(boundary.turning_points.size() != 1) {
        fmt::println(Logger::out("monotonicity"),"Ignoring boundary {} which has {} turning-points instead of 1 for pull_closest_corner()",slg.non_monotone_boundaries[non_monotone_boundary_index],boundary.turning_points.size()); Logger::out("monotonicity").flush();
        return;
    }
    const TurningPoint& tp = boundary.turning_points[0];
    
    // get labeling attribute and instanciate the halfedges API

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);
    CustomMeshHalfedges mesh_he(mesh);
    mesh_he.set_use_facet_region(attribute_name);
    #ifndef NDEBUG
        dump_vertex("current_tp",mesh_he.mesh(),tp.vertex(boundary,mesh_he.mesh()));
    #endif

    // find which corner of `boundary` is the closest of `tp`
    
    index_t closest_corner = tp.get_closest_corner(boundary,mesh_he);
    #ifndef NDEBUG
        dump_vertex("closest_corner",mesh,slg.corners[closest_corner].vertex);
    #endif

    // find which boundary around `closest_corner` is closest to `tp`, excluding `boundary`

    const Boundary& boundary_to_move = slg.boundaries[boundary.get_closest_boundary_of_turning_point(tp,closest_corner,mesh_he,slg.halfedge2boundary,slg.corners)];
    #ifndef NDEBUG
        Mesh boundary_to_move_file;
        boundary_to_move_file.copy(mesh,false,MESH_VERTICES);
        boundary_to_move_file.edges.create_edges((index_t) boundary_to_move.halfedges.size());
        FOR(he_index,boundary_to_move.halfedges.size()) {
            MeshHalfedges::Halfedge tmp_he = boundary_to_move.halfedges[he_index];
            boundary_to_move_file.edges.set_vertex(he_index,0,mesh.facet_corners.vertex(tmp_he.corner));
            mesh_he.move_to_opposite(tmp_he);
            boundary_to_move_file.edges.set_vertex(he_index,1,mesh.facet_corners.vertex(tmp_he.corner));
        }
        mesh_save(boundary_to_move_file,"boundary_to_move.geogram");
    #endif

    // find which label to assign between the boundary to remove and its "parallel" passing by the turning point
    
    index_t new_label = index_t(-1);
    // TODO check new_label computation
    if(tp.is_towards_right()) {
        new_label = slg.charts[boundary_to_move.right_chart].label;
    }
    else {
        new_label = slg.charts[boundary_to_move.left_chart].label;
    }

    // aggregate the facets of the left and right chart into a single set

    std::set<index_t> union_of_2_charts;
    union_of_2_charts.insert(slg.charts[boundary_to_move.left_chart].facets.begin(),slg.charts[boundary_to_move.left_chart].facets.end());
    union_of_2_charts.insert(slg.charts[boundary_to_move.right_chart].facets.begin(),slg.charts[boundary_to_move.right_chart].facets.end());
    #ifndef NDEBUG
        dump_facets("2_charts",mesh,union_of_2_charts);
    #endif

    // prepare a graph-cut optimization on both sides (2 charts)

    // TODO distance-based : the further a facet is, the bigger is the cost of changing its label. (too restrictive for this operator?)
    // TODO restrict the possible labels
    auto gcl = GraphCutLabeling(mesh,normals,(index_t) union_of_2_charts.size(),{0,1,2,3,4,5});
    for(index_t f : union_of_2_charts) {
        gcl.add_facet(f);
    }
    gcl.data_cost__set__fidelity_based(1);

    // Lock labels between the turning point and the corner, to ensure the corner moves where the turning point is
    // Note that MeshHalfedges::Halfedge::facet is the facet at the LEFT of the halfedge
    #ifndef NDEBUG
        std::set<index_t> facets_with_locked_label;
    #endif
    MeshHalfedges::Halfedge current_halfedge;
    if(closest_corner == boundary.start_corner) {
        // go across the boundary in the opposite way
        for(signed_index_t halfedge_index = (signed_index_t) (tp.outgoing_local_halfedge_index()-1); halfedge_index>= 0; halfedge_index--) {
            current_halfedge = boundary.halfedges[(std::vector<GEO::MeshHalfedges::Halfedge>::size_type) halfedge_index];
            if(tp.is_towards_right()) {
                mesh_he.move_to_opposite(current_halfedge);
            }
            if(!union_of_2_charts.contains(current_halfedge.facet)) {
                fmt::println(Logger::warn("monotonicity"),"pull_closest_corner() : facet {} ignored in GCoptimization because not a part of the 2 charts",current_halfedge.facet); Logger::warn("monotonicity").flush();
                continue;
            }
            gcl.data_cost__change_to__locked_polycube_label(current_halfedge.facet,new_label);
            #ifndef NDEBUG
                facets_with_locked_label.insert(current_halfedge.facet);
            #endif
        }
    }
    else {
        // go across the boundary in the same direction as the halfedges
        for(index_t halfedge_index = tp.outgoing_local_halfedge_index(); halfedge_index < boundary.halfedges.size(); halfedge_index++) {
            current_halfedge = boundary.halfedges[halfedge_index];
            if(tp.is_towards_right()) {
                mesh_he.move_to_opposite(current_halfedge);
            }
            if(!union_of_2_charts.contains(current_halfedge.facet)) {
                fmt::println(Logger::warn("monotonicity"),"pull_closest_corner() : facet {} ignored in GCoptimization because not a part of the 2 charts",current_halfedge.facet); Logger::warn("monotonicity").flush();
                continue;
            }
            gcl.data_cost__change_to__locked_polycube_label(current_halfedge.facet,new_label);
            #ifndef NDEBUG
                facets_with_locked_label.insert(current_halfedge.facet);
            #endif
        }
    }
    #ifndef NDEBUG
        fmt::println(Logger::out("monotonicity"),"facets_with_locked_label contains {} facets",facets_with_locked_label.size()); Logger::out("monotonicity").flush();
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
                    gcl.neighbors__change_to__shifted(f,adjacent_facet,(float) (1-dot_product)*500); // add a penalty the more ⟂ the edge is (// -> 0, ⟂ -> 300)
                    already_processed.insert({f, adjacent_facet});
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

void trace_contour(GEO::Mesh& mesh, const std::vector<vec3>& normals, const char* attribute_name, StaticLabelingGraph& slg, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);
    Attribute<double> per_facet_fidelity(mesh.facets.attributes(), "fidelity");
    // there are no iterators for GEO::Attribute, so I can't use std::min_element()...
    double min_fidelity = 1.0;
    index_t facet_of_min_fidelity = index_t(-1);
    FOR(f,mesh.facets.nb()) {
        if(per_facet_fidelity[f] < min_fidelity) {
            min_fidelity = per_facet_fidelity[f];
            facet_of_min_fidelity = f;
        }
    }
    const Chart& chart_to_refine = slg.charts[slg.facet2chart[facet_of_min_fidelity]];
    index_t prefered_label = nearest_label(normals[facet_of_min_fidelity]);
    geo_assert(are_orthogonal_labels(chart_to_refine.label,prefered_label));
    label[facet_of_min_fidelity] = prefered_label;
    // propagate prefered label to adjacent triangles, if they also prefer this label over the current one
    index_t current_facet = facet_of_min_fidelity,
            adjacent_facet = index_t(-1);
    std::set<index_t> created_chart_facets;
    std::queue<index_t> facets_to_explore; // FIDO of facet indices, only contains facets that prefers the new label (test before insertion)
    facets_to_explore.push(current_facet);
    while(!facets_to_explore.empty()) {
        current_facet = facets_to_explore.front();
        facets_to_explore.pop();
        if(created_chart_facets.contains(current_facet)) {
            continue; // skip this facet, already explored since insertion in FIFO
        }
        label[current_facet] = prefered_label;
        created_chart_facets.insert(current_facet);
        FOR(le,3) {
            adjacent_facet = mesh.facets.adjacent(current_facet,le);
            if(chart_to_refine.facets.contains(adjacent_facet) && is_better_label(normals[adjacent_facet],label[adjacent_facet],prefered_label) && !created_chart_facets.contains(adjacent_facet)) {
                facets_to_explore.push(adjacent_facet);
            }
        }
    }
    #ifndef NDEBUG
        fmt::println(Logger::out("refinement"),"created chart has {} facets",created_chart_facets.size()); Logger::out("refinement").flush();
        dump_facets("created_chart",mesh,created_chart_facets);
    #endif
    slg.fill_from(mesh,attribute_name,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);
    index_t created_chart_index = slg.facet2chart[facet_of_min_fidelity];
    geo_assert(slg.charts[created_chart_index].facets == created_chart_facets);
    // fill holes inside the created chart
    unsigned int nb_removed_charts = 0;
    do {
        nb_removed_charts = remove_surrounded_charts(mesh,attribute_name,slg);
        slg.fill_from(mesh,attribute_name,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);
    } while(nb_removed_charts!=0);
    created_chart_index = slg.facet2chart[facet_of_min_fidelity]; // created_chart_index may have changed -> update it
    index_t total_nb_halfedges = 0,
            he_counter = 0;
    for(index_t b : slg.charts[created_chart_index].boundaries) {
        total_nb_halfedges += (index_t) slg.boundaries[b].halfedges.size();
    }
    CustomMeshHalfedges mesh_he(mesh);
    std::vector<std::pair<MeshHalfedges::Halfedge,index_t>> all_boundary_halfedges(total_nb_halfedges); // associate a direction in {0,1,2}={X,Y,Z}
    std::array<double,3> per_axis_dot_product; // in order to find the max dot product -> axis the most aligned with the current halfedge
    for(index_t b : slg.charts[created_chart_index].boundaries) {
        for(const auto& he : slg.boundaries[b].halfedges) { // for each halfedge of the current boundary
            FOR(axis,3) { // X, Y and Z
                per_axis_dot_product[axis] = std::max(
                    dot(normalize(halfedge_vector(mesh,he)),label2vector[axis*2]),  // dot prodoct with X/Y/Z (according to axis)
                    dot(normalize(halfedge_vector(mesh,he)),label2vector[axis*2+1]) // dot prodoct with -X/-Y/-Z (according to axis)
                    );
            }
            all_boundary_halfedges[he_counter] = std::make_pair(he,VECTOR_MAX_INDEX(per_axis_dot_product));
            he_counter++;
        }
    }
    #ifndef NDEBUG
        // dump all_boundary_halfedges with closest axis as attribute
        std::map<std::pair<GEO::index_t, GEO::index_t>, index_t> edges_and_attributes;
        index_t vertex_1 = index_t(-1),
                vertex_2 = index_t(-1);
        for(const auto& kv : all_boundary_halfedges) {
            vertex_1 = mesh.facet_corners.vertex(kv.first.corner);
            vertex_2 = mesh.facet_corners.vertex(mesh.facets.next_corner_around_facet(kv.first.facet,kv.first.corner));
            edges_and_attributes[std::make_pair(vertex_1,vertex_2)] = kv.second;
        }
        dump_edges("all_boundary_edges","axis",mesh,edges_and_attributes);
    #endif
    // count number of same-axis groups (by counting transitions between different axes)
    index_t nb_same_axis_groups = 0; // 0 transition encountered for now
    FOR(he_index,all_boundary_halfedges.size()-2) {
        if(all_boundary_halfedges[he_index].second != all_boundary_halfedges[he_index+1].second) {
            nb_same_axis_groups++;
        }
    }
    if (all_boundary_halfedges[0].second != all_boundary_halfedges[all_boundary_halfedges.size()-1].second) {
        // first and last edges (which are adjacent) have a different closest axis
        nb_same_axis_groups++;
    }
    nb_same_axis_groups = (nb_same_axis_groups == 0) ? 1 : nb_same_axis_groups; // if no transition -> 1 group (not possible in practice I think)
    fmt::println(Logger::out("refinement"),"boundary of the created chart contains {} group(s) of closest axis",nb_same_axis_groups); Logger::out("refinement").flush();
    // fix validity of the created chart
    if(nb_same_axis_groups == 2) {
        /*
         *  |\               |\ 
         *  | \              | \X
         *  |  \             |  \/ 
         * Z|   \Z    =>    Z|  /\ 
         *  |    \           |    \Z 
         *  |     \          |     \ 
         *  |______\         |______\ 
         *     X                X
         *
         * split a group in two and add an orthogonal separator group
         * start by finding the biggest angle between 2 edges having the same nearest axis
         */
        double biggest_angle = Numeric::min_float64(),
               current_angle = 0.0;
        index_t boundary_halfedge_index_before_biggest_angle = index_t(-1),
                vertex_at_biggest_angle = index_t(-1),
                init_axis = index_t(-1),
                axis_to_insert = index_t(-1);
        FOR(bhe,total_nb_halfedges) { // == all_boundary_halfedges.size()
            if(all_boundary_halfedges[bhe].second != all_boundary_halfedges[(bhe+1)%total_nb_halfedges].second) {
                // ignore this pair of adjacent edges, they are not assigned to the same axis
                continue;
            }
            // compute the angle between this boundary halfedge and the next one
            current_angle = angle(
                Geom::halfedge_vector(mesh,all_boundary_halfedges[bhe].first),
                Geom::halfedge_vector(mesh,all_boundary_halfedges[(bhe+1)%total_nb_halfedges].first)
            );
            if(current_angle > biggest_angle) {
                biggest_angle = current_angle;
                boundary_halfedge_index_before_biggest_angle = bhe;
                vertex_at_biggest_angle = Geom::halfedge_vertex_index_to(mesh,all_boundary_halfedges[bhe].first);
                init_axis = all_boundary_halfedges[bhe].second; // axis associated to the 2 edges
                // find the axis to insert between the 2 part of this group
                // among the 2 others axes, find the one that is best suited for the current and next edge (= the second choice)
                std::array<double,3> per_axis_dot_product;
                per_axis_dot_product.fill(Numeric::min_float64());
                FOR(other__axis,3) {
                    if(other__axis == init_axis) {
                        continue; // leave min_float64() value
                    }
                    // find the max dot product for this axis,
                    // for both the current edge and the next one,
                    // for both the positive and negative axis direction
                    per_axis_dot_product[other__axis] = std::max({
                        dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[bhe].first)),label2vector[other__axis*2]),
                        dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[bhe].first)),label2vector[other__axis*2+1]),
                        dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[(bhe+1)%total_nb_halfedges].first)),label2vector[other__axis*2]),
                        dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[(bhe+1)%total_nb_halfedges].first)),label2vector[other__axis*2+1])
                    });
                }
                axis_to_insert = (index_t) VECTOR_MAX_INDEX(per_axis_dot_product);
            }
        }
        geo_assert(vertex_at_biggest_angle != index_t(-1));
        dump_vertex("vertex_at_biggest_angle",mesh,vertex_at_biggest_angle);
        // axis_to_insert is orthogonal to init_axis because they are not equal
        // but axis_to_insert may not be equal to the axis of the 2nd group...
        fmt::println(Logger::out("refinement"),"axis_to_insert={}",axis_to_insert); Logger::out("refinement").flush();
        // find wich side of the vertex_at_biggest_angle should be assigned to axis_to_insert
        // -> compute the average dot product when incrementing indices and when decrementing indices, inside this group of edges currently assigned to the same axis
        double score_of_assigning_new_axis_upward = 0.0,
               score_of_assigning_new_axis_downward = 0.0;
        index_t count = 0;
        for(index_t bhe = boundary_halfedge_index_before_biggest_angle+1; bhe < total_nb_halfedges; bhe++) {
            if(all_boundary_halfedges[bhe].second != init_axis) {
                break;
            }
            score_of_assigning_new_axis_upward += std::max(
                dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[bhe].first)),label2vector[axis_to_insert*2]),
                dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[bhe].first)),label2vector[axis_to_insert*2+1])
            );
            count++;
        }
        score_of_assigning_new_axis_upward /= (double) count;
        count = 0;
        for(signed_index_t bhe = (signed_index_t) boundary_halfedge_index_before_biggest_angle; bhe >= 0; bhe--) {
            if(all_boundary_halfedges[(size_t) bhe].second != init_axis) {
                break;
            }
            score_of_assigning_new_axis_downward += std::max(
                dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[(size_t) bhe].first)),label2vector[axis_to_insert*2]),
                dot(normalize(Geom::halfedge_vector(mesh,all_boundary_halfedges[(size_t) bhe].first)),label2vector[axis_to_insert*2+1])
                );
            count++;
        }
        score_of_assigning_new_axis_downward /= (double) count;
        if(score_of_assigning_new_axis_upward > score_of_assigning_new_axis_downward) {
            fmt::println(Logger::out("refinement"),"better to insert the axis upward the vertex at biggest angle"); Logger::out("refinement").flush();
            dump_edge("where_to_insert_the_axis",mesh,all_boundary_halfedges[boundary_halfedge_index_before_biggest_angle+1].first);
        }
        else {
            fmt::println(Logger::out("refinement"),"better to insert the axis downward the vertex at biggest angle"); Logger::out("refinement").flush();
            dump_edge("where_to_insert_the_axis",mesh,all_boundary_halfedges[boundary_halfedge_index_before_biggest_angle].first);
        }
    }
    else if(nb_same_axis_groups == 3) {
        fmt::println(Logger::err("refinement"),"trace_contour operator cannot handle boundaries with 3 groups of closest axis (yet)",nb_same_axis_groups); Logger::err("refinement").flush();
        // TODO cancel modifications
        return;
    }
    else if(nb_same_axis_groups != 4) {
        fmt::println(Logger::err("refinement"),"trace_contour operator cannot handle boundaries with {} groups of closest axis",nb_same_axis_groups); Logger::err("refinement").flush();
        // TODO cancel modifications
        return;
        // Future work: kind of graph cut with increasing compactness until there are 4 groups
    }
    // created chart = front
    // counterclockwise : right, back then left
    // the top chart must have the same label as the initial chart
}