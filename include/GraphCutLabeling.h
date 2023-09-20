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

class GraphCutLabeling {

public:

    //// Initialization //////////////////

    /**
     * \brief Prepare a labeling computation with a Graph-Cut optimization
     * \param[in] mesh A surface triangle mesh
     * TODO pre-computed normals as optional argument
     */
    GraphCutLabeling(const Mesh& mesh);

    //// Data cost definition //////////////////

    void data_cost__set__fidelity_based(int fidelity);

    void data_cost__set__locked_labels(const Attribute<index_t>& per_facet_label);

    void data_cost__change_to__fidelity_based(index_t facet_index, int fidelity);
    
    void data_cost__change_to__locked_label(index_t facet_index, index_t locked_label);

    void data_cost__change_to__forbidden_label(index_t facet_index, index_t forbidden_label);

    void data_cost__change_to__scaled(index_t facet_index, index_t label, float factor);

    //// Smooth cost definition //////////////////

    void smooth_cost__set__default();

    void smooth_cost__set__prevent_opposite_neighbors();

    void smooth_cost__set__custom(const std::array<int,6*6>& smooth_cost);

    //// Neighbors definition //////////////////

    void neighbors__set__compactness_based(int compactness);

    //// Getters & debug //////////////////

    int data_cost__get__for_facet_and_label(index_t facet_index, index_t label) const;

    void dump_costs() const;

    //// Optimization //////////////////

    /**
     * \brief Launch optimizer and store result
     * \param[out] output_labeling Facet attribute in which the labeling will be stored
     * \param[in] max_nb_iterations Maximum number of iteration for the GCoptimization::expansion() function
     */
    void compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations = 2);

private:

    const Mesh& mesh_; // needed for the facet number (everywhere) and facet adjacency (neighborhood weights)
    std::vector<int> data_cost_; // #facet*6, cost of assigning a facet to a label
    std::array<int,6*6> smooth_cost_; // label adjacency cost
    GCoptimizationGeneralGraph gco_; // underlying optimizer
    std::vector<vec3> normals_; // facet normals, computed in the constructor
    bool data_cost_set_; // if the data cost is already set or not, required for compute_solution()
    bool smooth_cost_set_; // if the smooth cost is already set or not, required for compute_solution()
    bool neighbors_set_; // if the neighbors are already set or not, required for compute_solution()
};