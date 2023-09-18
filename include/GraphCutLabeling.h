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

#define HIGH_COST 10e4

using namespace GEO;

// to avoid redundancy with constructors, methods...
inline void _facet_neighbors__set_compactness_based(index_t facet_index, GCoptimizationGeneralGraph& gco, const vec3& facet_normal, const Mesh& mesh, int compact_coeff);
inline void _facet_data_cost__set_fidelity_based(index_t facet_index, std::vector<int>& data_cost, const vec3& facet_normal, int fidelity);
inline void _facet_data_cost__lock_label(index_t facet_index, index_t locked_label, std::vector<int>& data_cost);
inline void _smooth_cost__fill(std::array<int,6*6>& smooth_cost); // for now, no other arguments

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

    /**
     * \brief Prepare a labeling computation with a Graph-Cut optimization, with custom smooth cost
     * \param[in] mesh A surface triangle mesh
     * \param[in] compact_coeff Compactness coefficient
     * \param[in] fidelity_coeff Fidelity coefficient
     * \param[in] smooth_cost Smoothness cost
     */
    GraphCutLabeling(const Mesh& mesh, int compact_coeff, int fidelity_coeff, const std::array<int,6*6>& smooth_cost);

    /**
     * \brief Prepare a labeling computation with a Graph-Cut optimization, and enfore all the labels (useless optimization if set_fidelity_based_data_cost() is not used after)
     * \param[in] mesh A surface triangle mesh
     * \param[in] compact_coeff Compactness coefficient
     * \param[in] per_facet_locked_labels Enforced labels
     */
    GraphCutLabeling(const Mesh& mesh, int compact_coeff, const Attribute<index_t>& per_facet_locked_label);

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
     * \brief Set per label cost based on how far the label direction is from the facet normal
     * \param[in] facet_index Which facet
     * \param[in] fidelity_coeff Fidelity coefficient
     */
    void set_fidelity_based_data_cost(index_t facet_index, int fidelity_coeff);

    /**
     * \brief Get the current data cost of assigning a given label on a given facet
     * \param[in] facet_index Which facet
     * \param[in] label Which label
     * \return The current data cost value
     */
    int get_data_cost(index_t facet_index, index_t label);

    /**
     * \brief Set the data cost of assigning a given label on a given facet
     * \param[in] facet_index Which facet
     * \param[in] label Which label
     * \param[in] value The new data cost value
     */
    void set_data_cost(index_t facet_index, index_t label, int value);

    //// Debug //////////////////

    void dump_costs();

    //// Optimization //////////////////

    /**
     * \brief Launch optimizer and store result
     * \param[out] output_labeling Facet attribute in which the labeling will be stored
     * \param[in] max_nb_iterations Maximum number of iteration for the GCoptimization::expansion() function
     */
    void compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations = 2);

private:

    const index_t NB_FACETS_;
    const Mesh& mesh_; // needed in set_fidelity_based_data_cost()
    std::vector<int> data_cost_;
    std::array<int,6*6> smooth_cost_;
    GCoptimizationGeneralGraph gco_;
};