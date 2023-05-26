#include <geogram/mesh/mesh.h> // for GEO::Mesh
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <GCoptimization.h>

#include "labeling.h" // for label2vector

#include "GraphCutLabeling.h"


GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, int compact_coeff, int fidelity_coeff) : NB_FACETS(mesh.facets.nb()), HIGH_COST(10e4), gco( (GCoptimization::SiteID) mesh.facets.nb(),6) {

    // define facet adjacency on the graph, weight based on compactness coeff & dot product of the normals

    FOR(facet_index,NB_FACETS) {
        vec3 facet_normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
        FOR(le,3) { // for each local edge of the current facet
            index_t neighbor_index = mesh.facets.adjacent(facet_index,le);
            vec3 neighbor_normal = normalize(Geom::mesh_facet_normal(mesh,neighbor_index)); // TODO precompute
            double dot = (GEO::dot(facet_normal,neighbor_normal)-1)/0.25;
            double cost = std::exp(-(1./2.)*std::pow(dot,2));
            gco.setNeighbors( (GCoptimization::SiteID) facet_index, (GCoptimization::SiteID) neighbor_index, (int) (compact_coeff*100*cost));
        }
    }

    // cost of assigning a facet to a label, weight based on fidelity coeff & dot product between normal & label direction

    data_cost.resize(NB_FACETS*6);
    FOR(facet_index,NB_FACETS) {
        vec3 normal = normalize(Geom::mesh_facet_normal(mesh,facet_index));
        FOR(label,6) {
            double dot = (GEO::dot(normal,label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost[facet_index*6+label] = (int) (fidelity_coeff*100*cost);
        }
    }

    // cost of assigning two labels to adjacent facets ?

    smooth_cost.resize(6*6); // 6 labels x 6 labels
    FOR(label1,6) {
        FOR(label2,6) {
            // same label = very smooth edge, different label = less smooth
            smooth_cost[label1+label2*6] = (label1==label2) ? 0 : 1;
        }
    }

    data_cost_initial = data_cost; // copy the initial data cost, if the user wants to recover them (see restore_initial_weights_of_facet())
}

void GraphCutLabeling::lock_label_on_facet(index_t facet_index, index_t locked_label) {
    // zero-cost for the locked label, high cost for other labels
    FOR(label,6) {
        data_cost[facet_index*6+label] = (label==locked_label) ? 0 : HIGH_COST;
    }
}

void GraphCutLabeling::forbid_label_on_facet(index_t facet_index, index_t forbidden_label) {
    // high cost for this label
    data_cost[facet_index*6+forbidden_label] = HIGH_COST;
}

void GraphCutLabeling::lock_all_facets(const Attribute<index_t>& per_facet_locked_label) {
    geo_assert(per_facet_locked_label.size()==NB_FACETS);
    FOR(facet_index,NB_FACETS) {
        lock_label_on_facet(facet_index,per_facet_locked_label[facet_index]);
    }
}

void GraphCutLabeling::restore_initial_weights_of_facet(index_t facet_index) {
    FOR(label,6) {
        data_cost[facet_index*6+label] = data_cost_initial[facet_index*6+label];
    }
}

void GraphCutLabeling::compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations) {
    geo_assert(output_labeling.size()==NB_FACETS);

    try {
        gco.setDataCost(data_cost.data());
        gco.setSmoothCost(smooth_cost.data());

        gco.expansion(max_nb_iterations);// run expansion
        // gco.swap(num_iterations) instead ?

        // get results
        FOR(facet,NB_FACETS)
            output_labeling[facet] = (index_t) gco.whatLabel( (GCoptimization::SiteID) facet);
    }
    catch (GCException e) {
		e.Report();
	}
}