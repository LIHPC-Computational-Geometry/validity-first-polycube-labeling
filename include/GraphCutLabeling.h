#pragma once

// from https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/graphcut_labeling.cpp
// which is based on ext/GraphCutOptimization/src/example.cpp GridGraph_DArraySArraySpatVarying()

#include <geogram/basic/numeric.h>      // for GEO::index_t
#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/mesh/mesh.h>          // for GEO::Mesh
#include <geogram/basic/vecg.h>         // for GEO::vecng

#include <GCoptimization.h>

#include <vector>
#include <array>
#include <map>
#include <utility>      // for std::pair

using namespace GEO;

typedef vecng<6, float> vec6f; // useful to store per label weights
typedef vecng<6, int> vec6i; // useful to store per label weights

// GCoptimizationGeneralGraph has 2 methods for neighbors weigths : setNeighbors() and setAllNeighbors()
// -> setNeighbors() will require redundancy: weights in GCoptimizationGeneralGraph given by value, weights in GraphCutLabeling, maybe also in GraphCutLabelingApp
// -> setAllNeighbors() read weights by address
// The class NeighborsCosts is meant to hide the complexity of the GCoptimizationGeneralGraph interface (set_neighbors() is like GCoptimization::setNeighbors()),
// while allowing to pass neighbors weights by address -> carefully respect the types
class NeighborsCosts {

public:
    NeighborsCosts();
    ~NeighborsCosts();

    index_t nb_sites() const;
    void set_nb_sites(index_t nb_sites);
    void set_neighbors(index_t facet1, index_t facet2, int cost);

    // see GCoptimizationGeneralGraph::setAllNeighbors() for the following 3 arrays
    // public because we need the pointers in compute_solution()
    GCoptimization::SiteID* per_facet_neighbors_; // for each site, gives the number of neighbors
    GCoptimization::SiteID** per_facet_neighbor_indices_; // for each site neighbors, gives the SiteID of the neighbor
    GCoptimization::EnergyTermType** per_facet_neighbor_weight_; // for each site neighbors, gives the weight between the site and the neighbor

private:
    index_t nb_sites_;
};

class GraphCutLabeling {

public:

    static const int HIGH_COST = 10e4;

    //// Initialization //////////////////

    /**
     * \brief Prepare a labeling computation with a Graph-Cut optimization
     * \param[in] mesh A surface triangle mesh
     */
    GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals);

    //// Data cost definition //////////////////

    void data_cost__set__fidelity_based(int fidelity);

    void data_cost__set__locked_labels(const Attribute<index_t>& per_facet_label);

    void data_cost__set__all_at_once(const std::vector<int>& data_cost);

    void data_cost__change_to__fidelity_based(index_t facet_index, int fidelity);
    
    void data_cost__change_to__locked_label(index_t facet_index, index_t locked_label);

    void data_cost__change_to__forbidden_label(index_t facet_index, index_t forbidden_label);

    void data_cost__change_to__scaled(index_t facet_index, index_t label, float factor);

    void data_cost__change_to__shifted(index_t facet_index, index_t label, float delta);

    void data_cost__change_to__per_label_weights(index_t facet_index, const vec6i& per_label_weights);

    //// Smooth cost definition //////////////////

    void smooth_cost__set__default();

    void smooth_cost__set__prevent_opposite_neighbors();

    void smooth_cost__set__custom(const std::array<int,6*6>& smooth_cost);

    //// Neighbors definition //////////////////

    void neighbors__set__compactness_based(int compactness);

    void neighbors__set__all_at_once(const NeighborsCosts& neighbors_costs);

    //// Getters & debug //////////////////

    vec6i data_cost__get__for_facet(index_t facet_index) const;

    int data_cost__get__for_facet_and_label(index_t facet_index, index_t label) const;

    void dump_costs() const;

    //// Optimization //////////////////

    /**
     * \brief Launch optimizer and store result
     * \param[out] output_labeling Facet attribute in which the labeling will be stored
     * \param[in] max_nb_iterations Maximum number of iteration for the GCoptimization::expansion() function
     */
    void compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations = 2);

    //// Static functions (for class methods & for data cost managed externally) //////////////////

    static void fill_data_cost__fidelity_based(const Mesh& mesh, const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity);

    static void shift_data_cost(std::vector<int>& data_cost, index_t facet_index, index_t label, float delta);

    static vec6i per_facet_data_cost_as_vector(const std::vector<int>& data_cost, index_t facet_index);

    static void fill_neighbors_cost__compactness_based(const Mesh& mesh, const std::vector<vec3>& normals, int compactness, NeighborsCosts& neighbors_costs);

private:

    const Mesh& mesh_; // needed for the facet number (everywhere) and facet adjacency (neighborhood weights)
    std::vector<int> data_cost_; // #facet*6, cost of assigning a facet to a label
    std::array<int,6*6> smooth_cost_; // label adjacency cost
    NeighborsCosts neighbors_costs_;
    GCoptimizationGeneralGraph gco_; // underlying optimizer
    const std::vector<vec3>& normals_; // facet normals
    bool data_cost_set_; // if the data cost is already set or not, required for compute_solution()
    bool smooth_cost_set_; // if the smooth cost is already set or not, required for compute_solution()
    bool neighbors_set_; // if the neighbors are already set or not, required for compute_solution()
};