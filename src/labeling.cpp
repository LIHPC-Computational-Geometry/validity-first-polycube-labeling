#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <array>     // for std::array
#include <algorithm> // for std::max_element()
#include <iterator>  // for std::distance()

#include "labeling.h"

void naive_labeling(Mesh& mesh, const char* attribute_name) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // create a facet attribute in this mesh
    GEO::vec3 normal;
    std::array<double,6> weights;
    for(index_t f: mesh.facets) { // for each facet
        normal = Geom::mesh_facet_normal(mesh,f); // get 3 signed components {x,y,z} of the normal
        // from 3 signed to 6 unsigned components:
        weights = {
            normal.x < 0.0 ? 0.0 : normal.x,  // +X
            normal.x < 0.0 ? -normal.x : 0.0, // -X
            normal.y < 0.0 ? 0.0 : normal.y,  // +Y
            normal.y < 0.0 ? -normal.y : 0.0, // -Y
            normal.z < 0.0 ? 0.0 : normal.z,  // +Z
            normal.z < 0.0 ? -normal.z : 0.0  // -Z
        };
        label[f] = (index_t) std::distance(weights.begin(),std::max_element(weights.begin(),weights.end())); // get index of max
    }
}