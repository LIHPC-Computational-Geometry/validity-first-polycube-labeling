#include <geogram/mesh/mesh.h> // for GEO::Mesh
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <fmt/os.h> // for fmt::output_file

#include <GCoptimization.h>

#include "labeling.h" // for label2vector

#include "GraphCutLabeling.h"

inline void _facet_neighbors__set_compactness_based(index_t facet_index, GCoptimizationGeneralGraph& gco, const vec3& facet_normal, const Mesh& mesh, int compact_coeff) {
    // define facet adjacency on the graph, weight based on compactness coeff & dot product of the normals
    FOR(le,3) { // for each local edge of the current facet
        index_t neighbor_index = mesh.facets.adjacent(facet_index,le);
        vec3 neighbor_normal = normalize(Geom::mesh_facet_normal(mesh,neighbor_index)); // TODO precompute
        double dot = (GEO::dot(facet_normal,neighbor_normal)-1)/0.25;
        double cost = std::exp(-(1./2.)*std::pow(dot,2));
        gco.setNeighbors( (GCoptimization::SiteID) facet_index, (GCoptimization::SiteID) neighbor_index, (int) (compact_coeff*100*cost));
    }
}

inline void _facet_data_cost__set_fidelity_based(index_t facet_index, std::vector<int>& data_cost, const vec3& facet_normal, int fidelity_coeff) {
    geo_assert(facet_index*6+5 < data_cost.size());
     // cost of assigning a facet to a label, weight based on fidelity coeff & dot product between normal & label direction
    FOR(label,6) {
        double dot = (GEO::dot(facet_normal,label2vector[label]) - 1.0)/0.2;
        double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
        data_cost[facet_index*6+label] = (int) (fidelity_coeff*100*cost);
    }
}

inline void _facet_data_cost__lock_label(index_t facet_index, index_t locked_label, std::vector<int>& data_cost) {
    geo_assert(facet_index*6+5 < data_cost.size());
    // zero-cost for the locked label, high cost for other labels
    FOR(label,6) {
        data_cost[facet_index*6+label] = (label==locked_label) ? 0 : HIGH_COST;
    }
}

void _smooth_cost__fill(std::array<int,6*6>& smooth_cost) {
    // cost of assigning two labels to adjacent facets ?
    // TODO re-implement 'prevent_opposite_neighbors' mode
    FOR(label1,6) {
        FOR(label2,6) {
            // same label = very smooth edge, different label = less smooth
            smooth_cost[label1+label2*6] = (label1==label2) ? 0 : 1;
        }
    }
}

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, int compact_coeff, int fidelity_coeff) : NB_FACETS_(mesh.facets.nb()), mesh_(mesh), gco_( (GCoptimization::SiteID) mesh.facets.nb(),6) {
    
    data_cost_.resize(NB_FACETS_*6);
    FOR(facet_index,NB_FACETS_) {
        vec3 normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
        _facet_neighbors__set_compactness_based(facet_index,gco_,normal,mesh,compact_coeff);
        _facet_data_cost__set_fidelity_based(facet_index,data_cost_,normal,fidelity_coeff);
    }
    _smooth_cost__fill(smooth_cost_);
}

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, int compact_coeff, int fidelity_coeff, const std::array<int,6*6>& smooth_cost) : NB_FACETS_(mesh.facets.nb()), mesh_(mesh), gco_( (GCoptimization::SiteID) mesh.facets.nb(),6) {
    data_cost_.resize(NB_FACETS_*6);
    FOR(facet_index,NB_FACETS_) {
        vec3 normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
        _facet_neighbors__set_compactness_based(facet_index,gco_,normal,mesh,compact_coeff);
        _facet_data_cost__set_fidelity_based(facet_index,data_cost_,normal,fidelity_coeff);
    }
    smooth_cost_ = smooth_cost; // user-given smoothness cost
}

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, int compact_coeff, const Attribute<index_t>& per_facet_locked_label) : NB_FACETS_(mesh.facets.nb()), mesh_(mesh), gco_( (GCoptimization::SiteID) mesh.facets.nb(),6) {
    geo_assert(per_facet_locked_label.size()==NB_FACETS_);

    data_cost_.resize(NB_FACETS_*6);
    FOR(facet_index,NB_FACETS_) {
        vec3 normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
        _facet_neighbors__set_compactness_based(facet_index,gco_,normal,mesh,compact_coeff);
        _facet_data_cost__lock_label(facet_index, per_facet_locked_label[facet_index], data_cost_); // difference from the first constructor : data cost is independant of the fidelity coeff, all labels are enforced
    }
    _smooth_cost__fill(smooth_cost_);
}

void GraphCutLabeling::lock_label_on_facet(index_t facet_index, index_t locked_label) {
    _facet_data_cost__lock_label(facet_index, locked_label, data_cost_);
}

void GraphCutLabeling::forbid_label_on_facet(index_t facet_index, index_t forbidden_label) {
    // high cost for this label
    data_cost_[facet_index*6+forbidden_label] = HIGH_COST;
}

void GraphCutLabeling::set_fidelity_based_data_cost(index_t facet_index, int fidelity_coeff) {
    vec3 normal = normalize(Geom::mesh_facet_normal(mesh_,facet_index)); // TODO precompute
    _facet_data_cost__set_fidelity_based(facet_index,data_cost_,normal,fidelity_coeff);
}

int GraphCutLabeling::get_data_cost(index_t facet_index, index_t label) {
    return data_cost_[facet_index*6+label];
}

void GraphCutLabeling::set_data_cost(index_t facet_index, index_t label, int value) {
    data_cost_[facet_index*6+label] = value;
}

void GraphCutLabeling::dump_costs() {

    // write data costs
    auto ofs_data = fmt::output_file("data.csv");
    ofs_data.print("facet,+X,-X,+Y,-X,+Z,-Z\n");
    FOR(facet_index,NB_FACETS_) {
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
    geo_assert(output_labeling.size()==NB_FACETS_);

    try {
        gco_.setDataCost(data_cost_.data());
        gco_.setSmoothCost(smooth_cost_.data());

        gco_.expansion(max_nb_iterations);// run expansion
        // gco.swap(num_iterations) instead ?

        // get results
        FOR(facet,NB_FACETS_)
            output_labeling[facet] = (index_t) gco_.whatLabel( (GCoptimization::SiteID) facet);
    }
    catch (GCException e) {
		e.Report();
	}
}