#pragma once

// from https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/graphcut_labeling.cpp
// which is based on ext/GraphCutOptimization/src/example.cpp GridGraph_DArraySArraySpatVarying()

#include <geogram/basic/numeric.h>      // for GEO::index_t
#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/mesh/mesh.h>          // for GEO::Mesh

#include <GCoptimization.h>

#include <vector>
#include <array>
#include <map>

using namespace GEO;

class GraphCutLabeling {

public:

    //// Initialization //////////////////

    /**
     * \brief Prepare a labeling computation with a Graph-Cut optimization
     * \param[in] mesh A surface triangle mesh
     * \param[in] compact_coeff Compactness coefficient
     * \param[in] fidelity_coeff Fidelity coefficient
     */
    GraphCutLabeling(const Mesh& mesh, int compact_coeff, int fidelity_coeff);

    //// Tweak weigths //////////////////

    /**
     * \brief Enforce a label on a facet
     * \param[in] facet_index Which facet
     * \param[in] locked_label Which label
     */
    void lock_label_on_facet(index_t facet_index, index_t locked_label);

    /**
     * \brief Prevent a label from being assigned to a facet
     * \param[in] facet_index Which facet
     * \param[in] forbidden_label Which label
     */
    void forbid_label_on_facet(index_t facet_index, index_t forbidden_label);

    /**
     * \brief Enforce all the labels (useless optimization if restore_initial_weights_of_facet() is not used after)
     * \param[in] per_facet_locked_label Labels to lock on the facets
     */
    void lock_all_facets(const Attribute<index_t>& per_facet_locked_label);

    /**
     * \brief Restore initial weigths, computed in the constructor, for a given facet
     * \param[in] facet_index Which facet
     */
    void restore_initial_weights_of_facet(index_t facet_index);

    //// Optimization //////////////////

    /**
     * \brief Launch optimizer and store result
     * \param[in] output_labeling Facet attribute in which the labeling will be stored
     * \param[in] max_nb_iterations Maximum number of iteration for the GCoptimization::expansion() function
     */
    void compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations = 2);

private:

    const index_t NB_FACETS;
    const int HIGH_COST;
    std::vector<int> data_cost;
    std::vector<int> data_cost_initial; // not memory-friendly...
    std::vector<int> smooth_cost;
    GCoptimizationGeneralGraph gco;
};