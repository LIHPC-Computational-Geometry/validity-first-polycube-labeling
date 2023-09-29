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
    per_facet_neighbors_(nullptr),
    per_facet_neighbor_indices_(nullptr),
    per_facet_neighbor_weight_(nullptr),
    nb_sites_(0) {}

NeighborsCosts::~NeighborsCosts() {
    if(nb_sites_ == 0) {
        return; // nothing to dealloc
    }
    FOR(f,nb_sites_) {
        delete per_facet_neighbor_indices_[f];
        delete per_facet_neighbor_weight_[f];
    }
    delete per_facet_neighbor_indices_;
    delete per_facet_neighbor_weight_;
    delete per_facet_neighbors_;
}

index_t NeighborsCosts::nb_sites() const {
    return nb_sites_;
}

void NeighborsCosts::set_nb_sites(index_t nb_sites) {
    geo_assert(nb_sites_ == 0); // dont allow re-sizing, the memory allocation must be done once
    nb_sites_ = nb_sites;
    per_facet_neighbors_ = new GCoptimization::SiteID[nb_sites];
    per_facet_neighbor_indices_ = new GCoptimization::SiteID*[nb_sites];
    per_facet_neighbor_weight_ = new GCoptimization::EnergyTermType*[nb_sites];
    FOR(f,nb_sites) {
        per_facet_neighbors_[f] = (GCoptimization::SiteID) 0;
        per_facet_neighbor_indices_[f] = new GCoptimization::SiteID[6]; // a facet has at most 3 neighbors -> 6 directionnal weights
        per_facet_neighbor_weight_[f] = new GCoptimization::EnergyTermType[6]; // a facet has at most 3 neighbors -> 6 directionnal weights
    }
}

void NeighborsCosts::set_neighbors(index_t facet1, index_t facet2, int cost) {
    geo_assert(nb_sites_ != 0); // set_nb_sites() must heve been called before
    // TODO check if a weight was already set for the 2 -> overwrite instead of expand the arrays
    geo_assert(facet1 < nb_sites_);
    geo_assert(facet2 < nb_sites_);
    geo_assert(cost >= 0);
    GCoptimization::SiteID facet1_neighbor_index = per_facet_neighbors_[facet1];
    GCoptimization::SiteID facet2_neighbor_index = per_facet_neighbors_[facet2];
    geo_assert(facet1_neighbor_index <= 6); // else array overflow
    geo_assert(facet2_neighbor_index <= 6); // else array overflow
    per_facet_neighbor_indices_[facet1][facet1_neighbor_index] = (GCoptimization::SiteID) facet2;
    per_facet_neighbor_indices_[facet2][facet2_neighbor_index] = (GCoptimization::SiteID) facet1;
    per_facet_neighbor_weight_[facet1][facet1_neighbor_index] = cost;
    per_facet_neighbor_weight_[facet2][facet2_neighbor_index] = cost;
    per_facet_neighbors_[facet1]++;
    per_facet_neighbors_[facet2]++;
}

GraphCutLabeling::GraphCutLabeling(const Mesh& mesh, const std::vector<vec3>& normals)
      : mesh_(mesh),
        data_cost_(mesh.facets.nb()*6,0),
        neighbors_costs_(),
        gco_( (GCoptimization::SiteID) mesh.facets.nb(),6),
        normals_(normals),
        data_cost_set_(false),
        smooth_cost_set_(false),
        neighbors_set_(false) {
}

void GraphCutLabeling::data_cost__set__fidelity_based(int fidelity) {
    if(data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"data cost already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    fill_data_cost__fidelity_based(mesh_,normals_,data_cost_,fidelity);
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

void GraphCutLabeling::data_cost__set__all_at_once(const std::vector<int>& data_cost) {
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
    shift_data_cost(data_cost_,facet_index,label,delta);
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
    fill_neighbors_cost__compactness_based(mesh_,normals_,compactness,neighbors_costs_);
    neighbors_set_ = true;
}

void GraphCutLabeling::neighbors__set__all_at_once(const NeighborsCosts& neighbors_costs) {
    if(neighbors_set_) {
        fmt::println(Logger::err("graph-cut"),"neighbors already set, and can only be set once"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    neighbors_costs_ = neighbors_costs;
    neighbors_set_ = true;
}

vec6i GraphCutLabeling::data_cost__get__for_facet(index_t facet_index) const {
    if(!data_cost_set_) {
        fmt::println(Logger::err("graph-cut"),"getter of data cost called before setter"); Logger::err("graph-cut").flush();
        geo_assert_not_reached;
    }
    return per_facet_data_cost_as_vector(data_cost_,facet_index);
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

    // write neighbors and their weights
    auto ofs_neighbors_1 = fmt::output_file("per_facet_neighbors.csv");
    ofs_neighbors_1.print("facet,#neighbors\n");
    FOR(f,neighbors_costs_.nb_sites()) {
        // f and the number of neighbors of f
        ofs_neighbors_1.print("{},{}\n",
            f,
            neighbors_costs_.per_facet_neighbors_[f]
        );
    }
    ofs_neighbors_1.flush();
    ofs_neighbors_1.close();

    auto ofs_neighbors_2 = fmt::output_file("per_facet_neighbor_indices.csv");
    ofs_neighbors_2.print("facet,neighbors0,neighbors1,neighbors2,neighbors3,neighbors4,neighbors5\n");
    FOR(f,neighbors_costs_.nb_sites()) {
        ofs_neighbors_2.print("{},",f);
        // f and the SiteID (facet index) of neighbors of f
        FOR(n,neighbors_costs_.per_facet_neighbors_[f]) {
            ofs_neighbors_2.print("{}",neighbors_costs_.per_facet_neighbor_indices_[f][n]);
            if(n!=5) {
                ofs_neighbors_2.print(",");
            }
        }
        // fill the remaining of the columns
        for(index_t n = neighbors_costs_.per_facet_neighbors_[f]; n < 6; ++n) {
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
    ofs_neighbors_3.print("facet,neighbors0,neighbors1,neighbors2,neighbors3,neighbors4,neighbors5\n");
    FOR(f,neighbors_costs_.nb_sites()) {
        ofs_neighbors_3.print("{},",f);
        // f and the weight associated to each neighbor of f
        FOR(n,neighbors_costs_.per_facet_neighbors_[f]) {
            ofs_neighbors_3.print("{}",neighbors_costs_.per_facet_neighbor_weight_[f][n]);
            if(n!=5) {
                ofs_neighbors_3.print(",");
            }
        }
        // fill the remaining of the columns
        for(index_t n = neighbors_costs_.per_facet_neighbors_[f]; n < 6; ++n) {
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
        nb_neighbors = neighbors_costs_.per_facet_neighbors_[f];
        FOR(n,nb_neighbors) {
            subkey = fmt::format("{}",neighbors_costs_.per_facet_neighbor_indices_[f][n]);
            if(debug_json[key].contains(subkey)) { // if this neighbor was already encountered
                if(debug_json[key][subkey] != neighbors_costs_.per_facet_neighbor_weight_[f][n]) {
                    // hmmm the weight is not consistent -> add another entry by changing the key
                    fmt::println(
                        Logger::err("graph-cut"),
                        "In GraphCutLabeling::dump_costs(), inconsistent neighbor weight : {} and {} for neighbor n°{} of site {}",
                        debug_json[key][subkey],
                        neighbors_costs_.per_facet_neighbor_weight_[f][n],
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
            debug_json[key][subkey] = neighbors_costs_.per_facet_neighbor_weight_[f][n];
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
            neighbors_costs_.per_facet_neighbors_,
            neighbors_costs_.per_facet_neighbor_indices_,
            neighbors_costs_.per_facet_neighbor_weight_
        );
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

void GraphCutLabeling::fill_data_cost__fidelity_based(const Mesh& mesh, const std::vector<vec3>& normals, std::vector<int>& data_cost, int fidelity) {
    // cost of assigning a facet to a label, weight based on fidelity coeff & dot product between normal & label direction
    FOR(f,mesh.facets.nb()) {
        FOR(label,6) {
            double dot = (GEO::dot(normals[f],label2vector[label]) - 1.0)/0.2;
            double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
            data_cost[f*6+label] = (int) (fidelity*100*cost);
        }
    }
}

void GraphCutLabeling::shift_data_cost(std::vector<int>& data_cost, index_t facet_index, index_t label, float delta) {
    data_cost[facet_index*6+label] = (int) std::max(0.0f,(((float) data_cost[facet_index*6+label]) + delta)); // forbid negative cost, min set to 0
}

vec6i GraphCutLabeling::per_facet_data_cost_as_vector(const std::vector<int>& data_cost, index_t facet_index) {
    vec6i result;
    memcpy(result.data(),data_cost.data()+(facet_index*6),sizeof(int)*6); // what? you don't like memcpy?
    return result;
}

void GraphCutLabeling::fill_neighbors_cost__compactness_based(const Mesh& mesh, const std::vector<vec3>& normals, int compactness, NeighborsCosts& neighbors_costs) {
    neighbors_costs.set_nb_sites(mesh.facets.nb());
    FOR(facet_index,mesh.facets.nb()) {
        // define facet adjacency on the graph, weight based on compactness coeff & dot product of the normals
        FOR(le,3) { // for each local edge of the current facet
            index_t neighbor_index = mesh.facets.adjacent(facet_index,le);
            double dot = (GEO::dot(normals[facet_index],normals[neighbor_index])-1)/0.25;
            double cost = std::exp(-(1./2.)*std::pow(dot,2));
            neighbors_costs.set_neighbors(facet_index,neighbor_index,(int) (compactness*100*cost));
        }
    }
}