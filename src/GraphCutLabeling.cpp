#include <geogram/mesh/mesh.h> // for GEO::Mesh
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <fmt/os.h> // for fmt::output_file

#include <GCoptimization.h>

#include <algorithm> // for std::max()

#include "labeling.h" // for label2vector

#include "GraphCutLabeling.h"

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh) : mesh_(mesh), gco_( (GCoptimization::SiteID) mesh.facets.nb(),6) {
    data_cost_.resize(mesh_.facets.nb()*6);
    normals_.resize(mesh_.facets.nb());
    FOR(f,mesh_.facets.nb()) {
        normals_[f] = normalize(Geom::mesh_facet_normal(mesh,f));
    }
    data_cost_set_ = false;
    smooth_cost_set_ = false;
    neighbors_set_ = false;
}

void GraphCutLabeling::data_cost__set__fidelity_based(int fidelity) {
    if(data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"data cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    // cost of assigning a facet to a label, weight based on fidelity coeff & dot product between normal & label direction
    FOR(f,mesh_.facets.nb()) {
        FOR(label,6) {
            double dot = (GEO::dot(normals_[f],label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost_[f*6+label] = (int) (fidelity*100*cost);
        }
    }
    data_cost_set_ = true;
}

void GraphCutLabeling::data_cost__set__locked_labels(const Attribute<index_t>& per_facet_label) {
    if(data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"data cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    FOR(facet_index,mesh_.facets.nb()) {
        // zero-cost for the locked label, high cost for other labels
        FOR(label,6) {
            data_cost_[facet_index*6+label] = (label==per_facet_label[facet_index]) ? 0 : HIGH_COST;
        }
    }
    data_cost_set_ = true;
}

void GraphCutLabeling::data_cost__change_to__fidelity_based(index_t facet_index, int fidelity) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    FOR(label,6) {
        double dot = (GEO::dot(normals_[facet_index],label2vector[label]) - 1.0)/0.2;
        double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
        data_cost_[facet_index*6+label] = (int) (fidelity*100*cost);
    }
}

void GraphCutLabeling::data_cost__change_to__locked_label(index_t facet_index, index_t locked_label) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    FOR(label,6) {
        data_cost_[facet_index*6+label] = (label==locked_label) ? 0 : HIGH_COST;
    }
}

void GraphCutLabeling::data_cost__change_to__forbidden_label(index_t facet_index, index_t forbidden_label) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    // FOR(label,6) {
    //     data_cost_[facet_index*6+label] = (label==forbidden_label) ? HIGH_COST : 0;
    // }
    data_cost_[facet_index*6+forbidden_label] = HIGH_COST; // do not edit weights of other labels
}

void GraphCutLabeling::data_cost__change_to__scaled(index_t facet_index, index_t label, float factor) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    data_cost_[facet_index*6+label] = (int) (((float) data_cost_[facet_index*6+label]) * factor);
}

void GraphCutLabeling::data_cost__change_to__shifted(index_t facet_index, index_t label, float delta) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    data_cost_[facet_index*6+label] = (int) std::max(0.0f,(((float) data_cost_[facet_index*6+label]) + delta)); // forbid negative cost, min set to 0
}

void GraphCutLabeling::data_cost__change_to__per_label_weights(index_t facet_index, const vec6i& per_label_weights) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    memcpy(data_cost_.data()+(facet_index*6),per_label_weights.data(),sizeof(int)*6); // memcpy <3
}

void GraphCutLabeling::smooth_cost__set__default() {
    if(smooth_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    // cost of assigning two labels to adjacent facets
    FOR(label1,6) {
        FOR(label2,6) {
            // same label = very smooth edge, different label = less smooth
            smooth_cost_[label1+label2*6] = (label1==label2) ? 0 : 1;
        }
    }
    smooth_cost_set_ = true;
}

void GraphCutLabeling::smooth_cost__set__prevent_opposite_neighbors() {
    if(smooth_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    FOR(label1,6) {
        FOR(label2,6) {
            smooth_cost_[label1+label2*6] = (label1==label2) ? 0 : ( // if samel label -> null cost
                                            (label2axis(label1)==label2axis(label2)) ? HIGH_COST : 1 // if same axis but different label (= opposite labels) -> high cost, else small cost (1)
                                            );
        }
    }
    smooth_cost_set_ = true;
}

void GraphCutLabeling::smooth_cost__set__custom(const std::array<int,6*6>& smooth_cost) {
    if(smooth_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    smooth_cost_ = smooth_cost;
    smooth_cost_set_ = true;
}

void GraphCutLabeling::neighbors__set__compactness_based(int compactness) {
    if(neighbors_set_) {
        fmt::println(Logger::err("graph-cut"),"neighbors already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    FOR(facet_index,mesh_.facets.nb()) {
        // define facet adjacency on the graph, weight based on compactness coeff & dot product of the normals
        FOR(le,3) { // for each local edge of the current facet
            index_t neighbor_index = mesh_.facets.adjacent(facet_index,le);
            double dot = (GEO::dot(normals_[facet_index],normals_[neighbor_index])-1)/0.25;
            double cost = std::exp(-(1./2.)*std::pow(dot,2));
            gco_.setNeighbors( (GCoptimization::SiteID) facet_index, (GCoptimization::SiteID) neighbor_index, (int) (compactness*100*cost));
            // can only be set once, see GCoptimizationGeneralGraph::setNeighbors()
        }
    }
    neighbors_set_ = true;
}

vec6i GraphCutLabeling::data_cost__get__for_facet(index_t facet_index) const {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    vec6i result;
    memcpy(result.data(),data_cost_.data()+(facet_index*6),sizeof(int)*6); // what? you don't like memcpy?
    return result;
}

int GraphCutLabeling::data_cost__get__for_facet_and_label(index_t facet_index, index_t label) const {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    return data_cost_[facet_index*6+label];
}

void GraphCutLabeling::dump_costs() const {

    // write data costs
    auto ofs_data = fmt::output_file("data.csv");
    ofs_data.print("facet,+X,-X,+Y,-X,+Z,-Z\n");
    FOR(facet_index,mesh_.facets.nb()) {
        ofs_data.print("{},{},{},{},{},{},{}\n",
            facet_index,
            data_cost_[facet_index*6+0],
            data_cost_[facet_index*6+1],
            data_cost_[facet_index*6+2],
            data_cost_[facet_index*6+3],
            data_cost_[facet_index*6+4],
            data_cost_[facet_index*6+5]
            );
    }
    ofs_data.flush();
    ofs_data.close();

    // write smooth costs
    auto ofs_smooth = fmt::output_file("smooth.csv");
    ofs_smooth.print("label1,0,1,2,3,4,5\n");
    FOR(label1,6) {
        ofs_smooth.print("{},{},{},{},{},{},{}\n",
            label1,
            smooth_cost_[label1+0*6],
            smooth_cost_[label1+1*6],
            smooth_cost_[label1+2*6],
            smooth_cost_[label1+3*6],
            smooth_cost_[label1+4*6],
            smooth_cost_[label1+5*6]
            );
    }
    ofs_smooth.flush();
    ofs_smooth.close();
}

void GraphCutLabeling::compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations) {

    geo_assert(output_labeling.size()==mesh_.facets.nb());

    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"compute_solution() called before setting of the data cost"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(!smooth_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"compute_solution() called before setting of the smooth cost"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(!neighbors_set_) {
        fmt::println(Logger::err("graph-cut"),"compute_solution() called before setting of the neighbors"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }

    try {
        gco_.setDataCost(data_cost_.data());
        gco_.setSmoothCost(smooth_cost_.data());
        gco_.expansion(max_nb_iterations);// run expansion
        // gco.swap(num_iterations) instead ?

        // get results
        FOR(facet,mesh_.facets.nb())
            output_labeling[facet] = (index_t) gco_.whatLabel( (GCoptimization::SiteID) facet);
    }
    catch (GCException e) {
		e.Report();
	}
}