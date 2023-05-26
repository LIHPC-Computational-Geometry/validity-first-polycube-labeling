#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <GCoptimization.h>

#include <array>            // for std::array
#include <initializer_list> // for std::initializer_list
#include <algorithm>        // for std::max_element()
#include <iterator>         // for std::distance()

#include "labeling.h"
#include "LabelingGraph.h"
#include "containers.h"

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

void naive_labeling(Mesh& mesh, const char* attribute_name) {

    // use GEO::Geom::triangle_normal_axis() instead ?

    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    GEO::vec3 normal;
    for(index_t f: mesh.facets) { // for each facet
        normal = Geom::mesh_facet_normal(mesh,f);
        label[f] = nearest_label(normal);
    }
}

// from https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/graphcut_labeling.cpp
// which is based on ext/GraphCutOptimization/src/example.cpp GridGraph_DArraySArraySpatVarying()
void graphcut_labeling(Mesh& mesh, const char* labeling_attribute_name, const char* locked_labels_attribute_name, const char* forbidden_labels_attribute_name, int compact_coeff, int fidelity_coeff) {
    
    // facet attribute for the output labeling
    Attribute<index_t> label(mesh.facets.attributes(), labeling_attribute_name);

    // store number of facets
    index_t nb_facets = mesh.facets.nb();

    // fill data cost with fidelity-based value (label further away from normal -> higher cost)
    int *data = new int[nb_facets*6]; // because 6 labels
    FOR(facet_index,nb_facets) {
        vec3 normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
        FOR(label,6) {
            double dot = (GEO::dot(normal,label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data[facet_index*6+label] = (int) (fidelity_coeff*100*cost);
        }
    }

    // use data costs to enforce locked labels

    if(mesh.facets.attributes().is_defined(locked_labels_attribute_name)) {
        Attribute<index_t> locked_labels(mesh.facets.attributes(), locked_labels_attribute_name);
        if(min(locked_labels)!=LabelingGraph::UNDEFINED) {
            FOR(facet_index,nb_facets) { // for each facet
                if(locked_labels[facet_index]==LabelingGraph::UNDEFINED) continue; // nothing specified for this facet
                FOR(label,6) {
                    // zero-cost for the locked label, high cost for other labels
                    data[facet_index*6+label] = (label==locked_labels[facet_index]) ? 0 : 10e4;
                }
            }
        }
        else {
            fmt::println(Logger::out("fix_labeling"),"Warning : no locked labels given to graphcut_labeling() - the attribute is defined but is filled with UNDEFINED"); Logger::out("fix_labeling").flush();
        }
    }
    else {
        fmt::println(Logger::out("fix_labeling"),"Warning : no locked labels given to graphcut_labeling() - the attribute is not defined"); Logger::out("fix_labeling").flush();
    }

    // use data costs to enforce forbidden labels

    if(mesh.facets.attributes().is_defined(forbidden_labels_attribute_name)) {
        Attribute<index_t> forbidden_labels(mesh.facets.attributes(), forbidden_labels_attribute_name);
        if(min(forbidden_labels)!=LabelingGraph::UNDEFINED) {
            FOR(facet_index,nb_facets) { // for each facet
                if(forbidden_labels[facet_index]==LabelingGraph::UNDEFINED) continue; // nothing specified for this facet
                FOR(label,6) {
                    if (label==forbidden_labels[facet_index]) {
                        // high cost for this label
                        data[facet_index*6+label] = 10e4;
                    }
                }
            }
        }
        else {
            fmt::println(Logger::out("fix_labeling"),"Warning : no forbidden labels given to graphcut_labeling() - the attribute is defined but is filled with UNDEFINED"); Logger::out("fix_labeling").flush();
        }
    }
    else {
        fmt::println(Logger::out("fix_labeling"),"Warning : no forbidden labels given to graphcut_labeling() - the attribute is not defined"); Logger::out("fix_labeling").flush();
    }

    // smooth costs
    // TODO re-implement 'prevent_opposite_neighbors' mode

    int *smooth = new int[6*6]; // 6 labels x 6 labels
    FOR(label1,6) {
        FOR(label2,6) {
            // same label = very smooth edge, different label = less smooth
            smooth[label1+label2*6] = (label1==label2) ? 0 : 1;
        }
    }

    try{
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph( (GCoptimization::SiteID) nb_facets,6);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);
		
        FOR(facet_index,nb_facets) {
            vec3 facet_normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
            FOR(le,3) { // for each local edge of the current facet
                index_t neighbor_index = mesh.facets.adjacent(facet_index,le);
                vec3 neighbor_normal = normalize(Geom::mesh_facet_normal(mesh,neighbor_index)); // TODO precompute
                double dot = (GEO::dot(facet_normal,neighbor_normal)-1)/0.25;
                double cost = std::exp(-(1./2.)*std::pow(dot,2));
                gc->setNeighbors( (GCoptimization::SiteID) facet_index, (GCoptimization::SiteID) neighbor_index, (int) (compact_coeff*100*cost));
            }
        }

		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);

        // get results
        FOR(facet,nb_facets)
			label[facet] = (index_t) gc->whatLabel( (GCoptimization::SiteID) facet);

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

	delete [] smooth;
	delete [] data;

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
    CustomMeshHalfedges::Halfedge current_halfedge;
    for(index_t boundary_index : slg.invalid_boundaries) { // for each invalid boundary
        const Boundary& current_boundary = slg.boundaries[boundary_index];
        new_label = nearest_label(current_boundary.average_normal);

        // Change the labels around the start/end corner, but only in place of the 2 charts next to the boundary
        for(auto current_corner : std::initializer_list<index_t>({current_boundary.start_corner,current_boundary.end_corner})) {
            geo_assert(current_corner != LabelingGraph::UNDEFINED);
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

unsigned int fix_invalid_corners(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {
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
                vertex_normal += Geom::mesh_facet_normal(mesh,be.facet); // compute normal of associated facet, update vertex_normal
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
void remove_invalid_charts(GEO::Mesh& mesh, const char* attribute_name, const StaticLabelingGraph& slg) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // get labeling attribute

    // Fill invalid charts using a Graph-Cut optimization,
    // preventing existing label from being re-applied

    Attribute<index_t> locked_labels(mesh.facets.attributes(), "locked_label");
    Attribute<index_t> forbidden_labels(mesh.facets.attributes(), "forbidden_label");

    locked_labels.copy(label); // starting with all facets locked to the current label
    forbidden_labels.fill(LabelingGraph::UNDEFINED); // starting with no forbidden labels

    for(index_t chart_index : slg.invalid_charts) { // for each invalid chart

        for(index_t facet_index : slg.charts[chart_index].facets) { // for each facet inside this chart
            locked_labels[facet_index] = LabelingGraph::UNDEFINED; // unlock the label that was previously locked
            forbidden_labels[facet_index] = label[facet_index]; // prevent the label from staying the same
        }

        graphcut_labeling(mesh,attribute_name,"locked_label","forbidden_label",1,1);
    }
}