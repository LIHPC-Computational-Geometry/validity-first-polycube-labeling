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

typedef vecng<6, float> vec6f; // useful to store per label weights in `graphcut_labeling` app (in case of an optimization with 6 possible labels per facet)
typedef vecng<6, int> vec6i; // useful to store per label weights (in case of an optimization with 6 possible labels per facet)

// GCoptimizationGeneralGraph has 2 methods for neighbors weigths : setNeighbors() and setAllNeighbors()
// -> setNeighbors() will require redundancy: weights in GCoptimizationGeneralGraph given by value, weights in GraphCutLabeling, maybe also in GraphCutLabelingApp
// -> setAllNeighbors() read weights by address
// The class NeighborsCosts is meant to hide the complexity of the GCoptimizationGeneralGraph interface (set_neighbors() is like GCoptimization::setNeighbors()),
// while allowing to pass neighbors weights by address -> carefully respect the types
class NeighborsCosts {

public:
    NeighborsCosts();
    ~NeighborsCosts();

    GCoptimization::SiteID nb_sites() const;
    void set_nb_sites(GCoptimization::SiteID nb_sites);
    bool are_neighbors(GCoptimization::SiteID site1, GCoptimization::SiteID site2, GCoptimization::SiteID& site1_local_index); // check if site2 is among the neighbors of site1. Do not check for reciprocity
    void set_neighbors(GCoptimization::SiteID site1, GCoptimization::SiteID site2, GCoptimization::EnergyTermType cost); // ensure reciprocity
    GCoptimization::EnergyTermType get_neighbors_cost(GCoptimization::SiteID site1, GCoptimization::SiteID site2); // ensure reciprocity

    // see GCoptimizationGeneralGraph::setAllNeighbors() for the following 3 arrays
    // public because we need the pointers in compute_solution()
    GCoptimization::SiteID* per_site_neighbors_; // for each site, gives the number of neighbors
    GCoptimization::SiteID** per_site_neighbor_indices_; // for each site neighbors, gives the SiteID of the neighbor
    GCoptimization::EnergyTermType** per_site_neighbor_weight_; // for each site neighbors, gives the weight between the site and the neighbor

private:
    GCoptimization::SiteID nb_sites_;
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

    GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals, index_t nb_facets, std::vector<index_t> labels); // if graph-cut on a subset of the surface and/or a subset of the labels

    //// Facets definition //////////////////

    void add_facet(index_t facet);

    //// Data cost definition //////////////////

    void data_cost__set__fidelity_based(int fidelity);

    void data_cost__set__locked_labels(const Attribute<index_t>& per_facet_polycube_label);

    void data_cost__set__all_at_once(const std::vector<int>& data_cost);

    void data_cost__change_to__fidelity_based(index_t facet_index, int fidelity);

    void data_cost__change_to__locked_polycube_label(index_t facet_index, index_t locked_polycube_label);

    void data_cost__change_to__forbidden_polycube_label(index_t facet_index, index_t forbidden_polycube_label);

    void data_cost__change_to__scaled(index_t facet_index, index_t polycube_label, float factor);

    void data_cost__change_to__shifted(index_t facet_index, index_t polycube_label, float delta);

    void data_cost__change_to__per_label_weights(index_t facet_index, const vec6i& per_polycube_label_weights); // only usable when nb_labels_ == 6

    //// Smooth cost definition //////////////////

    void smooth_cost__set__default();

    void smooth_cost__set__prevent_opposite_neighbors();

    void smooth_cost__set__custom(const std::vector<int>& smooth_cost);

    //// Neighbors definition //////////////////

    void neighbors__set__compactness_based(int compactness);

    void neighbors__set__all_at_once(const NeighborsCosts& neighbors_costs);

    void neighbors__change_to__scaled(index_t facet1, index_t facet2, float factor);

    void neighbors__change_to__shifted(index_t facet1, index_t facet2, float delta);

    //// Getters & debug //////////////////

    vec6i data_cost__get__for_facet(index_t facet_index) const; // only usable when nb_labels_ == 6

    int data_cost__get__for_facet_and_polycube_label(index_t facet_index, index_t polycube_label) const;

    inline bool facets_defined() const;

    void dump_costs() const;

    //// Optimization //////////////////

    /**
     * \brief Launch optimizer and store result
     * \param[out] output_labeling Facet attribute in which the labeling will be stored
     * \param[in] max_nb_iterations Maximum number of iteration for the GCoptimization::expansion() function
     */
    void compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations = 2);

    //// Static functions (for class methods & for data cost managed externally) //////////////////

    static void fill_data_cost__fidelity_based(const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity, const Mesh& mesh); // only usable when optimization on the whole mesh

    static void fill_data_cost__fidelity_based(const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity, const std::map<index_t,GCoptimization::SiteID>& facet2siteID, const std::map<index_t,GCoptimization::LabelID>& polycubeLabel2labelID);

    static void shift_data_cost(std::vector<int>& data_cost, GCoptimization::SiteID siteID, GCoptimization::LabelID labelID, float delta);

    static vec6i per_siteID_data_cost_as_vector(const std::vector<int>& data_cost, GCoptimization::SiteID siteID, index_t nb_labels, index_t nb_siteID);

    static void fill_neighbors_cost__compactness_based(const Mesh& mesh, const std::vector<vec3>& normals, int compactness, NeighborsCosts& neighbors_costs, const std::map<index_t,GCoptimization::SiteID>& facet2siteID);

private:

    // mesh

    const Mesh& mesh_; // needed for the facet number (everywhere) and facet adjacency (neighborhood weights)
    const std::vector<vec3>& normals_; // facet normals

    // facets on which the optimization will operate

    index_t nb_facets_;
    index_t count_facets_; // count how many facets have been defined
    std::map<index_t,GCoptimization::SiteID> facet2siteID_; // map between facet index (what the user manage) and siteID (internal to GCO)

    // possible labels for facets

    index_t nb_labels_;
    std::map<index_t,GCoptimization::LabelID> polycubeLabel2labelID_; // map between polycube label (what the user manage) and labelID (internal to GCO)
    std::vector<index_t> labelID2polycubeLabel_; // map between labelID and polycube label

    // underlying optimizer
    
    GCoptimizationGeneralGraph gco_;

    // optimizer settings

    std::vector<int> data_cost_; // nb_facets_*nb_labels_ items, cost of assigning a facet to a label
    bool data_cost_defined_; // if the data cost is already set or not, required for compute_solution()

    std::vector<int> smooth_cost_; // nb_labels_*nb_labels_ items, label adjacency cost
    bool smooth_cost_defined_; // if the smooth cost is already set or not, required for compute_solution()

    NeighborsCosts neighbors_costs_;
    bool neighbors_defined_; // if the neighbors are already set or not, required for compute_solution()
};