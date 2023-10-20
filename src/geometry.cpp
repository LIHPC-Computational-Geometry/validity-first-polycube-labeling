#include <queue>

#include "geometry.h"

void per_facet_distance(const Mesh& mesh, std::map<index_t,unsigned int>& distance) {
    // thank you Trizalio https://stackoverflow.com/a/72437022
    // Compute the distance to facets for which distance[facet]==0, for all facets in the keys of 'distance'
    std::queue<index_t> facets_to_explore; // facets we have to visit to update the distance in their neighborhood. FIFO data structure
    index_t current_facet = index_t(-1),
            adjacent_facet = index_t(-1);
    unsigned int current_distance = 0; // will store the distance currently computed for the current facet
    for(const auto& kv : distance) {
        if(kv.second == 0) {
            facets_to_explore.push(kv.first);
        }
    }
    while(!facets_to_explore.empty()) { // while there still are facets to explore
        current_facet = facets_to_explore.front(); // pick the facet at the front of the FIFO
        facets_to_explore.pop(); // remove it from the FIFO
        current_distance = distance[current_facet];
        FOR(le,3) { // for each local edge of the current facet
            adjacent_facet = mesh.facets.adjacent(current_facet,le); // get the facet index of the neighbor at the other side of the current local edge
            if(!distance.contains(adjacent_facet)) {
                continue; // -> do not explore the adjacent facet, it is not inside the set of facets
            }
            if (current_distance + 1 < distance[adjacent_facet]) { // if the distance of the neighbor was too much (passing by the current facet is closer)
                distance[adjacent_facet] = current_distance + 1; // update the distance
                facets_to_explore.emplace(adjacent_facet);
            }
        }
    }
}