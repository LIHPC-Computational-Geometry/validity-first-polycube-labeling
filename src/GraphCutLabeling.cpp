#include <geogram/mesh/mesh.h> // for GEO::Mesh
#include <geogram/basic/geometry.h>     // for GEO::vec3

#include <fmt/os.h> // for fmt::output_file

#include <GCoptimization.h>

#include <nlohmann/json.hpp>

#include <utility>      // for std::pair
#include <map>          // for std::map
#include <algorithm>    // for std::max()

#include "labeling.h"   // for label2vector

#include "GraphCutLabeling.h"

NeighborsCosts::NeighborsCosts() :
    per_site_neighbors_(nullptr),
    per_site_neighbor_indices_(nullptr),
    per_site_neighbor_weight_(nullptr),
    nb_sites_(0) {}

NeighborsCosts::~NeighborsCosts() {
    if(nb_sites_ == 0) {
        return; // nothing to dealloc
    }
    FOR(s,nb_sites_) {
        delete per_site_neighbor_indices_[s];
        delete per_site_neighbor_weight_[s];
    }
    delete per_site_neighbor_indices_;
    delete per_site_neighbor_weight_;
    delete per_site_neighbors_;
}

GCoptimization::SiteID NeighborsCosts::nb_sites() const {
    return nb_sites_;
}

void NeighborsCosts::set_nb_sites(GCoptimization::SiteID nb_sites) {
    geo_assert(nb_sites_ == 0); // dont allow re-sizing, the memory allocation must be done once
    nb_sites_ = nb_sites;
    per_site_neighbors_ = new GCoptimization::SiteID[nb_sites];
    per_site_neighbor_indices_ = new GCoptimization::SiteID*[nb_sites];
    per_site_neighbor_weight_ = new GCoptimization::EnergyTermType*[nb_sites];
    FOR(s,nb_sites) {
        per_site_neighbors_[s] = (GCoptimization::SiteID) 0;
        per_site_neighbor_indices_[s] = new GCoptimization::SiteID[3]; // a facet has at most 3 neighbors
        per_site_neighbor_weight_[s] = new GCoptimization::EnergyTermType[3]; // a facet has at most 3 neighbors
    }
}

bool NeighborsCosts::are_neighbors(GCoptimization::SiteID site1, GCoptimization::SiteID site2, GCoptimization::SiteID& site1_neighbor_index) {
    FOR(current_site1_neighbor_index,per_site_neighbors_[site1]) {
        if(per_site_neighbor_indices_[site1][current_site1_neighbor_index] == site2) {
            site1_neighbor_index = (GCoptimization::SiteID) current_site1_neighbor_index; // = on which neighbor index of site1 is site2
            return true;
        }
    }
    return false;
}

void NeighborsCosts::set_neighbors(GCoptimization::SiteID site1, GCoptimization::SiteID site2, GCoptimization::EnergyTermType cost) {
    geo_assert(nb_sites_ != 0); // set_nb_sites() must heve been called before
    geo_assert(site1 < nb_sites_);
    geo_assert(site2 < nb_sites_);
    geo_assert(cost >= 0);
    GCoptimization::SiteID site1_neighbor_index = GCoptimization::SiteID(-1);
    GCoptimization::SiteID site2_neighbor_index = GCoptimization::SiteID(-1);
    if(are_neighbors(site1,site2,site1_neighbor_index)) {
        geo_assert(are_neighbors(site2,site1,site2_neighbor_index)); // if site1 is a neighbor of site2, site2 must be a neighbor of site1
        geo_assert(per_site_neighbor_weight_[site1][site1_neighbor_index] == per_site_neighbor_weight_[site2][site2_neighbor_index]); // assert the cost is consistent
        // overwrite the cost
        per_site_neighbor_weight_[site1][site1_neighbor_index] = cost;
        per_site_neighbor_weight_[site2][site2_neighbor_index] = cost;
    }
    else {
        // site1 and site2 are not yet neighbors
        site1_neighbor_index = per_site_neighbors_[site1];
        site2_neighbor_index = per_site_neighbors_[site2];
        geo_assert(site1_neighbor_index < 3); // else array overflow
        geo_assert(site2_neighbor_index < 3); // else array overflow
        per_site_neighbor_indices_[site1][site1_neighbor_index] = site2;
        per_site_neighbor_indices_[site2][site2_neighbor_index] = site1;
        per_site_neighbor_weight_[site1][site1_neighbor_index] = cost;
        per_site_neighbor_weight_[site2][site2_neighbor_index] = cost;
        per_site_neighbors_[site1]++;
        per_site_neighbors_[site2]++;
    }    
}

GCoptimization::EnergyTermType NeighborsCosts::get_neighbors_cost(GCoptimization::SiteID site1, GCoptimization::SiteID site2) {
    geo_assert(site1 < nb_sites_);
    geo_assert(site2 < nb_sites_);
    GCoptimization::SiteID site1_neighbor_index = GCoptimization::SiteID(-1);
    GCoptimization::SiteID site2_neighbor_index = GCoptimization::SiteID(-1);
    if(are_neighbors(site1,site2,site1_neighbor_index)) {
        geo_assert(are_neighbors(site2,site1,site2_neighbor_index)); // if site1 is a neighbor of site2, site2 must be a neighbor of site1
        geo_assert(per_site_neighbor_weight_[site1][site1_neighbor_index] == per_site_neighbor_weight_[site2][site2_neighbor_index]); // assert the cost is consistent
        return per_site_neighbor_weight_[site1][site1_neighbor_index]; // return the neighbor cost
    }
    // else: site1 and site2 are not neighbors, no cost to return
    geo_assert_not_reached; 
}

// graph-cut
//  - on the whole surface -> nb_facets = mesh.facets.nb()
//  and
//  - on all labels -> labels = {0,1,2,3,4,5}
GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals) : GraphCutLabeling(mesh, normals, mesh.facets.nb(), {0,1,2,3,4,5}) {}

// graph-cut
// - on a subset of the surface (only `nb_facets` triangles)
// and/or
// - on a subset of the labels (the ones in `labels`)
//
// But GCO always need contiguous values
// so for example, a an optimization on facet 46, 12, 16, 49 and 8
// with only the labels 2=+Y and 4=+Z would result in:
//
//   facet | siteID    polycube label | labelID
//  -------|--------   ---------------|---------
//      46 | 0                      2 | 0
//      12 | 1                      4 | 1
//      16 | 2
//      49 | 3                     X --> X = polycubeLabel2labelID_ (which is a map)
//       8 | 4                     X <-- X = labelID2polycubeLabel_ (which is a vector)
//
//      X --> X = facet2siteID_ (which is a map)
//      X <-- X = not needed

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals, index_t nb_facets, std::vector<index_t> labels)
      : mesh_(mesh),
        normals_(normals),
        nb_facets_(nb_facets),
        count_facets_(0), // the user has to call add_facet() nb_facets_ times
        // default initialization for facet2siteID_ (empty map)
        nb_labels_((index_t) labels.size()),
        // default initialization for polycubeLabel2labelID_ (empty map)
        labelID2polycubeLabel_(labels), // fill vector from initializer_list
        gco_( (GCoptimization::SiteID) nb_facets_, (GCoptimization::LabelID) labels.size()),
        data_cost_((std::vector<int>::size_type) nb_facets_*labels.size(),0),
        data_cost_defined_(false),
        smooth_cost_((std::vector<int>::size_type) nb_labels_*nb_labels_,1),
        smooth_cost_defined_(false),
        neighbors_costs_(),
        neighbors_defined_(false) {

    // manage definition of sites
    if(nb_facets_ == mesh.facets.nb()) {
        // graph-cut on the whole surface :
        FOR(f,mesh.facets.nb()) {
            facet2siteID_[f] = (GCoptimization::SiteID) f;
        }
        count_facets_ = nb_facets_; // makes facets_defined() true
    }

    // manage definition of labels
    // basically we need the reverse map of labelID2polycubeLabel_, for compute_solution()
    geo_assert(nb_labels_ > 1); // an optimization with 0 or 1 label makes no sense
    FOR(i,nb_labels_) {
        geo_assert(labelID2polycubeLabel_[i] <= 5); // assert is a valid polycube label
        geo_assert(!polycubeLabel2labelID_.contains(labelID2polycubeLabel_[i])); // prevent duplication of a label
        polycubeLabel2labelID_[labelID2polycubeLabel_[i]] = (GCoptimization::SiteID) i;
    }
}

void GraphCutLabeling::add_facet(index_t facet) {
    if(facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"all facets already defined, add_facet() cannot be called after"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(!facet2siteID_.contains(facet)); // prevent re-insertion of a facet
    // register the facet and store the association between facet index and site ID
    facet2siteID_[facet] = (GCoptimization::SiteID) count_facets_;
    count_facets_++;
}

void GraphCutLabeling::data_cost__set__fidelity_based(int fidelity) {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"data cost already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    fill_data_cost__fidelity_based(normals_,data_cost_,fidelity,facet2siteID_,polycubeLabel2labelID_);
    data_cost_defined_ = true;
}

void GraphCutLabeling::data_cost__set__locked_labels(const Attribute<index_t>& per_facet_polycube_label) {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"data cost already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    // this method must only be called when graph-cut is applied on the whole mesh
    geo_assert(nb_facets_ == mesh_.facets.nb());
    geo_assert(nb_facets_ == per_facet_polycube_label.size());
    FOR(siteID,mesh_.facets.nb()) {
        geo_assert(polycubeLabel2labelID_.contains(per_facet_polycube_label[siteID])); // assert per_facet_polycube_label[siteID] is among the possible polycube labels
        // per_facet_label link a facet index (here equal to siteID) to a polycube label
        // but we need a LabelID -> convert with polycubeLabel2labelID_
        FOR(labelID,nb_labels_) {
            // zero-cost for the locked label, high cost for other labels
            data_cost_[siteID*nb_labels_+labelID] = ( (GCoptimization::LabelID) labelID == polycubeLabel2labelID_[per_facet_polycube_label[siteID]]) ? 0 : HIGH_COST;
        }
    }
    data_cost_defined_ = true;
}

void GraphCutLabeling::data_cost__set__all_at_once(const std::vector<int>& data_cost) {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"data cost already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(data_cost.size()==data_cost_.size());
    data_cost_ = data_cost;
    data_cost_defined_ = true;
}

void GraphCutLabeling::data_cost__change_to__fidelity_based(index_t facet_index, int fidelity) {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being defined yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_[facet_index] < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    FOR(labelID,nb_labels_) { // for each labelID
        double dot = (GEO::dot(normals_[facet_index],label2vector[labelID2polycubeLabel_[labelID]]) - 1.0)/0.2;
        double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
        data_cost_[
            (std::vector<int>::size_type) facet2siteID_[facet_index] * (std::vector<int>::size_type) nb_labels_ +
            (std::vector<int>::size_type) labelID
        ] = (int) (fidelity*100*cost);
    }
}

void GraphCutLabeling::data_cost__change_to__locked_polycube_label(index_t facet_index, index_t locked_polycube_label) {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being defined yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_[facet_index] < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(polycubeLabel2labelID_.contains(locked_polycube_label)); // assert locked_polycube_label is among the possible polycube labels
    FOR(labelID,nb_labels_) { // for each labelID
        data_cost_[
            (std::vector<int>::size_type) facet2siteID_[facet_index] * (std::vector<int>::size_type) nb_labels_ +
            (std::vector<int>::size_type) labelID
        ] = ( (GCoptimization::SiteID) labelID == polycubeLabel2labelID_[locked_polycube_label]) ? 0 : HIGH_COST;
    }
}

void GraphCutLabeling::data_cost__change_to__forbidden_polycube_label(index_t facet_index, index_t forbidden_polycube_label) {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being defined yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_[facet_index] < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(polycubeLabel2labelID_.contains(forbidden_polycube_label)); // assert forbidden_polycube_label is among the possible polycube labels
    data_cost_[
        (std::vector<int>::size_type) facet2siteID_[facet_index] * (std::vector<int>::size_type) nb_labels_ +
        (std::vector<int>::size_type) polycubeLabel2labelID_[forbidden_polycube_label]
    ] = HIGH_COST;
    // do not edit weights of other labels
}

void GraphCutLabeling::data_cost__change_to__scaled(index_t facet_index, index_t polycube_label, float factor) {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being defined yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_[facet_index] < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(polycubeLabel2labelID_.contains(polycube_label)); // assert polycube_label is among the possible polycube labels
    std::vector<int>::size_type index_in_data_cost =
        (std::vector<int>::size_type) facet2siteID_[facet_index] * (std::vector<int>::size_type) nb_labels_ +
        (std::vector<int>::size_type) polycubeLabel2labelID_[polycube_label];
    data_cost_[index_in_data_cost] = (int) (((float) data_cost_[index_in_data_cost]) * factor);
}

void GraphCutLabeling::data_cost__change_to__shifted(index_t facet_index, index_t polycube_label, float delta) {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being defined yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_[facet_index] < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(polycubeLabel2labelID_.contains(polycube_label)); // assert polycube_label is among the possible polycube labels
    shift_data_cost(data_cost_,facet2siteID_[facet_index],polycubeLabel2labelID_[polycube_label],delta);
}

void GraphCutLabeling::data_cost__change_to__per_label_weights(index_t facet_index, const vec6i& per_polycube_label_weights) {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being defined yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_[facet_index] < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(nb_labels_ == 6); // this method expect GCO on all labels (returns a vec6i)
    memcpy(data_cost_.data()+(facet2siteID_[facet_index]*6),per_polycube_label_weights.data(),sizeof(int)*6);
}

void GraphCutLabeling::smooth_cost__set__default() {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(smooth_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    // cost of assigning two labels to adjacent facets
    FOR(labelID1,nb_labels_) {
        FOR(labelID2,nb_labels_) {
            // same labelID = very smooth edge = zero cost, different labelID = less smooth = small cost
            smooth_cost_[(std::vector<int>::size_type) labelID1 + (std::vector<int>::size_type) labelID2 * (std::vector<int>::size_type) nb_labels_] = (labelID1==labelID2) ? 0 : 1;
        }
    }
    smooth_cost_defined_ = true;
}

void GraphCutLabeling::smooth_cost__set__prevent_opposite_neighbors() {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(smooth_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    FOR(labelID1,nb_labels_) {
        FOR(labelID2,nb_labels_) {
            smooth_cost_[(std::vector<int>::size_type) labelID1 + (std::vector<int>::size_type) labelID2 * (std::vector<int>::size_type) nb_labels_] =
                (labelID1==labelID2) ? 0 : ( // if samel label -> null cost
                (label2axis(labelID2polycubeLabel_[labelID1])==label2axis(labelID2polycubeLabel_[labelID2])) ? HIGH_COST : 1 // if same axis but different label (= opposite labels) -> high cost, else small cost (1)
                );
        }
    }
    smooth_cost_defined_ = true;
}

void GraphCutLabeling::smooth_cost__set__custom(const std::vector<int>& smooth_cost) {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(smooth_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(smooth_cost.size()==smooth_cost_.size());
    smooth_cost_ = smooth_cost;
    smooth_cost_defined_ = true;
}

void GraphCutLabeling::neighbors__set__compactness_based(int compactness) {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(neighbors_defined_) {
        fmt::println(Logger::err("graph-cut"),"neighbors already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    fill_neighbors_cost__compactness_based(mesh_,normals_,compactness,neighbors_costs_,facet2siteID_);
    neighbors_defined_ = true;
}

void GraphCutLabeling::neighbors__set__all_at_once(const NeighborsCosts& neighbors_costs) {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(neighbors_defined_) {
        fmt::println(Logger::err("graph-cut"),"neighbors already defined, and can only be defined once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    neighbors_costs_ = neighbors_costs;
    neighbors_defined_ = true;
}

void GraphCutLabeling::neighbors__change_to__scaled(index_t facet1, index_t facet2, float factor) {
    if(!neighbors_defined_) {
        fmt::println(Logger::err("graph-cut"),"the neighbors cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    float previous_cost = (float) neighbors_costs_.get_neighbors_cost(facet2siteID_[facet1],facet2siteID_[facet2]);
    neighbors_costs_.set_neighbors(facet2siteID_[facet1],facet2siteID_[facet2], (GCoptimization::EnergyTermType) (previous_cost*factor));
}

void GraphCutLabeling::neighbors__change_to__shifted(index_t facet1, index_t facet2, float delta) {
    if(!neighbors_defined_) {
        fmt::println(Logger::err("graph-cut"),"the neighbors cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    float previous_cost = (float) neighbors_costs_.get_neighbors_cost(facet2siteID_[facet1],facet2siteID_[facet2]);
    neighbors_costs_.set_neighbors(facet2siteID_[facet1],facet2siteID_[facet2], (GCoptimization::EnergyTermType) std::max(0.0f,previous_cost+delta)); // forbid negative cost, min set to 0
}

vec6i GraphCutLabeling::data_cost__get__for_facet(index_t facet_index) const {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_.at(facet_index) < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(nb_labels_ == 6); // this method expect GCO on all labels (returns a vec6i)
    return per_siteID_data_cost_as_vector(data_cost_,facet2siteID_.at(facet_index),nb_labels_,nb_facets_);
}

int GraphCutLabeling::data_cost__get__for_facet_and_polycube_label(index_t facet_index, index_t polycube_label) const {
    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index)); // `facet_index` should be among the facets defined for the optimization
    geo_assert(facet2siteID_.at(facet_index) < (GCoptimization::SiteID) nb_facets_); // the corresponding siteID should be in the range [0:nb_facets_[ = [0:#siteID[
    geo_assert(polycubeLabel2labelID_.contains(polycube_label)); // assert polycube_label is among the possible polycube labels
    geo_assert(polycubeLabel2labelID_.at(polycube_label) < (GCoptimization::LabelID) nb_labels_);
    return data_cost_[
        (std::vector<int>::size_type) facet2siteID_.at(facet_index) * (std::vector<int>::size_type) nb_labels_ +
        (std::vector<int>::size_type) polycubeLabel2labelID_.at(polycube_label)
    ];
}

bool GraphCutLabeling::facets_defined() const {
    return count_facets_ == nb_facets_;
}

void GraphCutLabeling::dump_costs() const {
    if(!facets_defined()) {
        fmt::println(Logger::err("graph-cut"),"not all facets are defined, you need to call add_facet() until all facets are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }

    //////////////////////////////////////////////////
	// Write data costs as CSV
	//////////////////////////////////////////////////

    auto ofs_data = fmt::output_file("data.csv");
    // print header. for ex: "facet,siteID,+X,+Y,-Z". siteID is the internal index for GCO
    ofs_data.print("facet,siteID");
    FOR(labelID,nb_labels_) {
        ofs_data.print(",{}",LABEL2STR(labelID2polycubeLabel_[labelID]));
    }
    ofs_data.print("\n");
    // print data costs values
    for(const auto& kv : facet2siteID_) { // first = key = facet index, second = value = siteID
        ofs_data.print("{},{}",kv.first,kv.second);
        FOR(labelID,nb_labels_) {
            ofs_data.print(",{}",
                data_cost_[
                    (std::vector<int>::size_type) kv.second * (std::vector<int>::size_type) nb_labels_ +
                    (std::vector<int>::size_type) labelID
                ]
            );
        }
        ofs_data.print("\n");
    }
    ofs_data.flush();
    ofs_data.close();

    //////////////////////////////////////////////////
	// Write smooth costs as CSV
	//////////////////////////////////////////////////

    auto ofs_smooth = fmt::output_file("smooth.csv");
    // print header. for ex: "label1,+X,+Y,-Z"
    ofs_smooth.print("label1");
    FOR(labelID_2,nb_labels_) {
        ofs_smooth.print(",{}",LABEL2STR(labelID2polycubeLabel_[labelID_2]));
    }
    ofs_smooth.print("\n");
    // print smooth costs values
    FOR(labelID_1,nb_labels_) {
        ofs_data.print("{}",LABEL2STR(labelID2polycubeLabel_[labelID_1]));
        FOR(labelID_2,nb_labels_) {
            ofs_smooth.print(",{}",smooth_cost_[labelID_1+labelID_2*nb_labels_]);
        }
        ofs_smooth.print("\n");
    }
    ofs_smooth.flush();
    ofs_smooth.close();

    //////////////////////////////////////////////////
	// Write neighbors and their weights as CSV
	//////////////////////////////////////////////////

    auto ofs_neighbors_1 = fmt::output_file("per_facet_neighbors.csv");
    ofs_neighbors_1.print("site,#neighbors\n");
    FOR(f,neighbors_costs_.nb_sites()) {
        // f and the number of neighbors of f
        ofs_neighbors_1.print("{},{}\n",
            f,
            neighbors_costs_.per_site_neighbors_[f]
        );
    }
    ofs_neighbors_1.flush();
    ofs_neighbors_1.close();

    auto ofs_neighbors_2 = fmt::output_file("per_facet_neighbor_indices.csv");
    ofs_neighbors_2.print("site,neighbors0,neighbors1,neighbors2,neighbors3,neighbors4,neighbors5\n");
    FOR(f,neighbors_costs_.nb_sites()) {
        ofs_neighbors_2.print("{},",f);
        // f and the SiteID (facet index) of neighbors of f
        FOR(n,neighbors_costs_.per_site_neighbors_[f]) {
            ofs_neighbors_2.print("{}",neighbors_costs_.per_site_neighbor_indices_[f][n]);
            if(n!=5) {
                ofs_neighbors_2.print(",");
            }
        }
        // fill the remaining of the columns
        for(index_t n = (index_t) neighbors_costs_.per_site_neighbors_[f]; n < 6; ++n) {
            ofs_neighbors_2.print(" "); // out of bound -> space
            if(n!=5) {
                ofs_neighbors_2.print(",");
            }
        }
        ofs_neighbors_2.print("\n");
    }
    ofs_neighbors_2.flush();
    ofs_neighbors_2.close();

    auto ofs_neighbors_3 = fmt::output_file("per_facet_neighbor_weights.csv");
    ofs_neighbors_3.print("site,neighbors0,neighbors1,neighbors2,neighbors3,neighbors4,neighbors5\n");
    FOR(f,neighbors_costs_.nb_sites()) {
        ofs_neighbors_3.print("{},",f);
        // f and the weight associated to each neighbor of f
        FOR(n,neighbors_costs_.per_site_neighbors_[f]) {
            ofs_neighbors_3.print("{}",neighbors_costs_.per_site_neighbor_weight_[f][n]);
            if(n!=5) {
                ofs_neighbors_3.print(",");
            }
        }
        // fill the remaining of the columns
        for(index_t n = (index_t) neighbors_costs_.per_site_neighbors_[f]; n < 6; ++n) {
            ofs_neighbors_3.print(" "); // out of bound -> space
            if(n!=5) {
                ofs_neighbors_3.print(",");
            }
        }
        ofs_neighbors_3.print("\n");
    }
    ofs_neighbors_3.flush();
    ofs_neighbors_3.close();

    //////////////////////////////////////////////////
	// Write neighbors and their weights as CSV
	//////////////////////////////////////////////////
    
    // also write neighbors as JSON
    nlohmann::json debug_json;
    std::string key, subkey;
    GCoptimizationGeneralGraph::SiteID nb_neighbors = 0;
    FOR(f,neighbors_costs_.nb_sites()) {
        key = fmt::format("{}",f);
        nb_neighbors = neighbors_costs_.per_site_neighbors_[f];
        FOR(n,nb_neighbors) {
            subkey = fmt::format("{}",neighbors_costs_.per_site_neighbor_indices_[f][n]);
            if(debug_json[key].contains(subkey)) { // if this neighbor was already encountered
                if(debug_json[key][subkey] != neighbors_costs_.per_site_neighbor_weight_[f][n]) {
                    // hmmm the weight is not consistent -> add another entry by changing the key
                    fmt::println(
                        Logger::err("graph-cut"),
                        "In GraphCutLabeling::dump_costs(), inconsistent neighbor weight : {} and {} for neighbor n°{} of site {}",
                        debug_json[key][subkey],
                        neighbors_costs_.per_site_neighbor_weight_[f][n],
                        n,
                        f
                    );
                    subkey += "_";
                }
                else {
                    // no need to re-write the same value in debug_json
                    break;
                }
            }
            debug_json[key][subkey] = neighbors_costs_.per_site_neighbor_weight_[f][n];
        }
    }
    std::ofstream ofs_json("neighbors.json");
    ofs_json << std::setw(4) << debug_json << std::endl;
    ofs_json.close();
}

void GraphCutLabeling::compute_solution(Attribute<index_t>& output_labeling, int max_nb_iterations) {

    geo_assert(output_labeling.size()==mesh_.facets.nb());

    if(!data_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"compute_solution() called before definition of the data cost"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(!smooth_cost_defined_) {
        fmt::println(Logger::err("graph-cut"),"compute_solution() called before definition of the smooth cost"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(!neighbors_defined_) {
        fmt::println(Logger::err("graph-cut"),"compute_solution() called before definition of the neighbors"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }

    #ifndef NDEBUG
        fmt::println(Logger::out("graph-cut"),"Launching GCO optimizer with {} facet and {} labels",nb_facets_,nb_labels_); Logger::out("graph-cut").flush();
    #endif

    try {
        gco_.setDataCost(data_cost_.data());
        gco_.setSmoothCost(smooth_cost_.data());
        gco_.setAllNeighbors(
            neighbors_costs_.per_site_neighbors_,
            neighbors_costs_.per_site_neighbor_indices_,
            neighbors_costs_.per_site_neighbor_weight_
        );
        gco_.expansion(max_nb_iterations);// run expansion
        // gco.swap(num_iterations) instead ?

        // get results
        for(auto kv : facet2siteID_)
            output_labeling[kv.first] = labelID2polycubeLabel_[ (std::vector<index_t>::size_type) gco_.whatLabel(kv.second)];
    }
    catch (GCException e) {
		e.Report();
	}
}

void GraphCutLabeling::fill_data_cost__fidelity_based(const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity, const Mesh& mesh) {
    // TODO avoid redundancy with the other fill_data_cost__fidelity_based() ?
    geo_assert(data_cost.size()==mesh.facets.nb()*6);
    FOR(f,mesh.facets.nb()) {
        FOR(label,6) {
            double dot = (GEO::dot(normals[f],label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost[f*6+label] = (int) (fidelity*100*cost);
        }
    }
}

void GraphCutLabeling::fill_data_cost__fidelity_based(const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity, const std::map<index_t,GCoptimization::SiteID>& facet2siteID, const std::map<index_t,GCoptimization::LabelID>& polycubeLabel2labelID) {
    // cost of assigning a facet to a label, weight based on fidelity coeff & dot product between normal & label direction
    index_t nb_labels = (index_t) polycubeLabel2labelID.size();
    geo_assert(data_cost.size()==facet2siteID.size()*nb_labels); // data_cost should have #facets * #labels elements
    for(const auto& sites : facet2siteID) { // first = facet index, second = siteID
        for(const auto& labels : polycubeLabel2labelID) { // first = polycube label, second = labelID
            double dot = (GEO::dot(normals[sites.first],label2vector[labels.first]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost[
                (std::vector<int>::size_type) sites.second * (std::vector<int>::size_type) nb_labels +
                (std::vector<int>::size_type) labels.second
            ] = (int) (fidelity*100*cost);
        }
    }
}

void GraphCutLabeling::shift_data_cost(std::vector<int>& data_cost, GCoptimization::SiteID siteID, GCoptimization::LabelID labelID, float delta) {
    std::vector<int>::size_type index_in_data_cost = (std::vector<int>::size_type) siteID * 6 + (std::vector<int>::size_type) labelID;
    data_cost[index_in_data_cost] = (int) std::max(0.0f,(((float) data_cost[index_in_data_cost]) + delta));
}

vec6i GraphCutLabeling::per_siteID_data_cost_as_vector(const std::vector<int>& data_cost, GCoptimization::SiteID siteID, index_t nb_labels, index_t nb_siteID) {
    geo_assert(nb_labels == 6); // this method expect GCO on all labels (returns a vec6i)
    geo_assert(siteID < (GCoptimization::SiteID) nb_siteID);
    vec6i result;
    memcpy(result.data(),data_cost.data()+(siteID*6),sizeof(int)*6);
    return result;
}

void GraphCutLabeling::fill_neighbors_cost__compactness_based(const Mesh& mesh, const std::vector<vec3>& normals, int compactness, NeighborsCosts& neighbors_costs, const std::map<index_t,GCoptimization::SiteID>& facet2siteID) {
    neighbors_costs.set_nb_sites((GCoptimization::SiteID) facet2siteID.size());
    for(auto kv : facet2siteID) {
        GCoptimization::SiteID current_site = kv.second;
        FOR(le,3) { // for each local edge of the current facet
            index_t neighbor_index = mesh.facets.adjacent(kv.first,le);
            if(!facet2siteID.contains(neighbor_index)) { // the neighbor is not a part of the sites
                continue;
            }
            GCoptimization::SiteID neighbor_site = facet2siteID.at(neighbor_index);
            geo_assert(neighbor_site != GCoptimization::SiteID(-1));
            // the normals are in [-1.0:1.0]^3
            // the dot product of the normals is in [-1:1]
            // applying (x-1)/0.25 makes it in [-4:0]
            double dot = (GEO::dot(normals[kv.first],normals[neighbor_index])-1)/0.25;
            // apply a Gaussian function of mean 0
            double cost = std::exp(-(1./2.)*std::pow(dot,2));
            neighbors_costs.set_neighbors(current_site, neighbor_site, (GCoptimization::EnergyTermType) (compactness*100*cost));
        }
    }
}