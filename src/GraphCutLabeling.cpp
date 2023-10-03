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
        per_site_neighbor_indices_[s] = new GCoptimization::SiteID[6]; // a facet has at most 3 neighbors -> 6 directionnal weights
        per_site_neighbor_weight_[s] = new GCoptimization::EnergyTermType[6]; // a facet has at most 3 neighbors -> 6 directionnal weights
    }
}

void NeighborsCosts::set_neighbors(GCoptimization::SiteID site1, GCoptimization::SiteID site2, GCoptimization::EnergyTermType cost) {
    geo_assert(nb_sites_ != 0); // set_nb_sites() must heve been called before
    // TODO check if a weight was already set for the 2 -> overwrite instead of expand the arrays
    geo_assert(site1 < nb_sites_);
    geo_assert(site2 < nb_sites_);
    geo_assert(cost >= 0);
    GCoptimization::SiteID site1_neighbor_index = per_site_neighbors_[site1];
    GCoptimization::SiteID site2_neighbor_index = per_site_neighbors_[site2];
    geo_assert(site1_neighbor_index <= 6); // else array overflow
    geo_assert(site2_neighbor_index <= 6); // else array overflow
    per_site_neighbor_indices_[site1][site1_neighbor_index] = site2;
    per_site_neighbor_indices_[site2][site2_neighbor_index] = site1;
    per_site_neighbor_weight_[site1][site1_neighbor_index] = cost;
    per_site_neighbor_weight_[site2][site2_neighbor_index] = cost;
    per_site_neighbors_[site1]++;
    per_site_neighbors_[site2]++;
}

// graph-cut on the whole surface -> nb_sites = mesh.facets.nb()
GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals) : GraphCutLabeling(mesh,normals, (GCoptimization::SiteID) mesh.facets.nb()) {}

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals, GCoptimization::SiteID nb_sites)
      : mesh_(mesh),
        nb_sites_(nb_sites),
        count_sites_(nb_sites),
        data_cost_(nb_sites*6,0),
        neighbors_costs_(),
        gco_( (GCoptimization::SiteID) nb_sites,6),
        normals_(normals),
        data_cost_set_(false),
        smooth_cost_set_(false),
        neighbors_set_(false) {
    if(nb_sites_ != (GCoptimization::SiteID) mesh.facets.nb()) {
        count_sites_ = 0; // now the user has to call add_site() nb_sites times
    }
    else {
        // graph-cut on the whole surface :
        FOR(f,mesh.facets.nb()) {
            facet2siteID_[f] = f;
        }
        // leave count_sites_ == nb_sites_ (making sites_set_() true)
    }
}

void GraphCutLabeling::add_site(index_t facet) {
    if(sites_set()) {
        fmt::println(Logger::err("graph-cut"),"sites already set, add_site() cannot be called after"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(!facet2siteID_.contains(facet)); // prevent re-insertion of a facet
    // register the facet and store the association between facet index and site ID
    facet2siteID_[facet] = count_sites_;
    count_sites_++;
}

void GraphCutLabeling::data_cost__set__fidelity_based(int fidelity) {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"data cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    fill_data_cost__fidelity_based(mesh_,normals_,data_cost_,fidelity,facet2siteID_);
    data_cost_set_ = true;
}

void GraphCutLabeling::data_cost__set__locked_labels(const Attribute<index_t>& per_facet_label) {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"data cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(nb_sites_ == per_facet_label.size()); // this method must only be called when graph-cut is applied on the whole mesh
    FOR(siteID,mesh_.facets.nb()) {
        // zero-cost for the locked label, high cost for other labels
        FOR(label,6) {
            data_cost_[siteID*6+label] = (label==per_facet_label[siteID]) ? 0 : HIGH_COST;
        }
    }
    data_cost_set_ = true;
}

void GraphCutLabeling::data_cost__set__all_at_once(const std::vector<int>& data_cost) {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"data cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(data_cost.size()==data_cost_.size());
    data_cost_ = data_cost;
    data_cost_set_ = true;
}

void GraphCutLabeling::data_cost__change_to__fidelity_based(index_t facet_index, int fidelity) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(facet2siteID_.contains(facet_index));
    geo_assert(facet2siteID_[facet_index] < nb_sites_);
    FOR(label,6) {
        double dot = (GEO::dot(normals_[facet_index],label2vector[label]) - 1.0)/0.2;
        double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
        data_cost_[facet2siteID_[facet_index]*6+label] = (int) (fidelity*100*cost);
    }
}

void GraphCutLabeling::data_cost__change_to__locked_label(index_t facet_index, index_t locked_label) {
    data_cost__change_to__locked_label(facet2siteID_[facet_index],locked_label);
}

void GraphCutLabeling::data_cost__change_to__locked_label(GCoptimization::SiteID siteID, index_t locked_label) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(siteID < nb_sites_);
    FOR(label,6) {
        data_cost_[siteID*6+label] = (label==locked_label) ? 0 : HIGH_COST;
    }
}

void GraphCutLabeling::data_cost__change_to__forbidden_label(GCoptimization::SiteID siteID, index_t forbidden_label) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(siteID < nb_sites_);
    data_cost_[siteID*6+forbidden_label] = HIGH_COST; // do not edit weights of other labels
}

void GraphCutLabeling::data_cost__change_to__scaled(GCoptimization::SiteID siteID, index_t label, float factor) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(siteID < nb_sites_);
    data_cost_[siteID*6+label] = (int) (((float) data_cost_[siteID*6+label]) * factor);
}

void GraphCutLabeling::data_cost__change_to__shifted(GCoptimization::SiteID siteID, index_t label, float delta) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(siteID < nb_sites_);
    shift_data_cost(data_cost_,siteID,label,delta);
}

void GraphCutLabeling::data_cost__change_to__per_label_weights(GCoptimization::SiteID siteID, const vec6i& per_label_weights) {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"the data cost cannot be change because it has not being set yet"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(siteID < nb_sites_);
    memcpy(data_cost_.data()+(siteID*6),per_label_weights.data(),sizeof(int)*6); // memcpy <3
}

void GraphCutLabeling::smooth_cost__set__default() {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
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
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
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
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(smooth_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"smooth cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    smooth_cost_ = smooth_cost;
    smooth_cost_set_ = true;
}

void GraphCutLabeling::neighbors__set__compactness_based(int compactness) {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(neighbors_set_) {
        fmt::println(Logger::err("graph-cut"),"neighbors already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    fill_neighbors_cost__compactness_based(mesh_,normals_,compactness,neighbors_costs_,facet2siteID_);
    neighbors_set_ = true;
}

void GraphCutLabeling::neighbors__set__all_at_once(const NeighborsCosts& neighbors_costs) {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    if(neighbors_set_) {
        fmt::println(Logger::err("graph-cut"),"neighbors already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    neighbors_costs_ = neighbors_costs;
    neighbors_set_ = true;
}

vec6i GraphCutLabeling::data_cost__get__for_site(GCoptimization::SiteID siteID) const {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    return per_site_data_cost_as_vector(data_cost_,siteID);
}

int GraphCutLabeling::data_cost__get__for_site_and_label(GCoptimization::SiteID siteID, index_t label) const {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    geo_assert(siteID < nb_sites_);
    return data_cost_[siteID*6+label];
}

bool GraphCutLabeling::sites_set() const {
    return count_sites_ == nb_sites_;
}

void GraphCutLabeling::dump_costs() const {
    if(!sites_set()) {
        fmt::println(Logger::err("graph-cut"),"not all sites are set, you need to call add_site() until all sites are defined"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }

    // TODO if nb_sites_ != #facets, also dump facet2siteID_

    // write data costs
    auto ofs_data = fmt::output_file("data.csv");
    ofs_data.print("site,+X,-X,+Y,-X,+Z,-Z\n");
    FOR(s,nb_sites_) {
        ofs_data.print("{},{},{},{},{},{},{}\n",
            s,
            data_cost_[s*6+0],
            data_cost_[s*6+1],
            data_cost_[s*6+2],
            data_cost_[s*6+3],
            data_cost_[s*6+4],
            data_cost_[s*6+5]
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

    // write neighbors and their weights
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
        for(index_t n = neighbors_costs_.per_site_neighbors_[f]; n < 6; ++n) {
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
        for(index_t n = neighbors_costs_.per_site_neighbors_[f]; n < 6; ++n) {
            ofs_neighbors_3.print(" "); // out of bound -> space
            if(n!=5) {
                ofs_neighbors_3.print(",");
            }
        }
        ofs_neighbors_3.print("\n");
    }
    ofs_neighbors_3.flush();
    ofs_neighbors_3.close();
    
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
        gco_.setAllNeighbors(
            neighbors_costs_.per_site_neighbors_,
            neighbors_costs_.per_site_neighbor_indices_,
            neighbors_costs_.per_site_neighbor_weight_
        );
        gco_.expansion(max_nb_iterations);// run expansion
        // gco.swap(num_iterations) instead ?

        // get results
        for(auto kv : facet2siteID_)
            output_labeling[kv.first] = (index_t) gco_.whatLabel(kv.second);
    }
    catch (GCException e) {
		e.Report();
	}
}

void GraphCutLabeling::fill_data_cost__fidelity_based(const Mesh& mesh, const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity) {
    // TODO avoid redundancy with the other fill_data_cost__fidelity_based()
    geo_assert(data_cost.size()==mesh.facets.nb()*6);
    FOR(f,mesh.facets.nb()) {
        FOR(label,6) {
            double dot = (GEO::dot(normals[f],label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost[f*6+label] = (int) (fidelity*100*cost);
        }
    }
}

void GraphCutLabeling::fill_data_cost__fidelity_based(const Mesh& mesh, const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity, const std::map<index_t,GCoptimization::SiteID>& facet2siteID) {
    // cost of assigning a facet to a label, weight based on fidelity coeff & dot product between normal & label direction
    geo_assert(data_cost.size()==facet2siteID.size()*6);
    for(auto kv : facet2siteID) {
        FOR(label,6) {
            double dot = (GEO::dot(normals[kv.first],label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost[kv.second*6+label] = (int) (fidelity*100*cost);
        }
    }
}

void GraphCutLabeling::shift_data_cost(std::vector<int>& data_cost, GCoptimization::SiteID site, index_t label, float delta) {
    geo_assert(site*6+label < data_cost.size());
    data_cost[site*6+label] = (int) std::max(0.0f,(((float) data_cost[site*6+label]) + delta)); // forbid negative cost, min set to 0
}

vec6i GraphCutLabeling::per_site_data_cost_as_vector(const std::vector<int>& data_cost, GCoptimization::SiteID site) {
    vec6i result;
    memcpy(result.data(),data_cost.data()+(site*6),sizeof(int)*6); // what? you don't like memcpy?
    return result;
}

void GraphCutLabeling::fill_neighbors_cost__compactness_based(const Mesh& mesh, const std::vector<vec3>& normals, int compactness, NeighborsCosts& neighbors_costs, const std::map<index_t,GCoptimization::SiteID>& facet2siteID) {
    neighbors_costs.set_nb_sites(facet2siteID.size());
    for(auto kv : facet2siteID) {
        GCoptimization::SiteID current_site = kv.second;
        FOR(le,3) { // for each local edge of the current facet
            index_t neighbor_index = mesh.facets.adjacent(kv.first,le);
            if(!facet2siteID.contains(neighbor_index)) { // the neighbor is not a part of the sites
                continue;
            }
            GCoptimization::SiteID neighbor_site = facet2siteID.at(neighbor_index);
            geo_assert(neighbor_site != GCoptimization::SiteID(-1));
            double dot = (GEO::dot(normals[kv.first],normals[neighbor_index])-1)/0.25;
            double cost = std::exp(-(1./2.)*std::pow(dot,2));
            neighbors_costs.set_neighbors(current_site, neighbor_site, (GCoptimization::EnergyTermType) (compactness*100*cost));
        }
    }
}