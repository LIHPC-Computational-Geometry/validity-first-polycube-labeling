#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <array>     // for std::array
#include <algorithm> // for std::max_element()
#include <iterator>  // for std::distance()

#include "labeling.h"

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name) {

    //open the file
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        fmt::println(Logger::out("I/O"),"Could not open labeling file '{}'",filename); 
        return false;
    }

    fmt::println(Logger::out("I/O"),"Loading file {}...",filename); Logger::out("I/O") << std::flush;

    std::string current_line;
    unsigned long current_label; // necessary type for string to unsigned int conversion (stoul)
    GEO::index_t current_line_number = 0;
    GEO::index_t facets_number = mesh.facets.nb(); // expected line number (one line per facet)

    // Add an attribute on facets. see https://github.com/BrunoLevy/geogram/wiki/Mesh#attributes
    // The type must be Attribute<index_t> to be used with geogram/mesh/mesh_halfedges.h . See MeshHalfedges::facet_region_
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);

    //fill the attribute
    while (ifs >> current_line) {
        try {
            current_label = std::stoul(current_line);
            if(current_label > 5) {
                fmt::println(Logger::err("I/O"),"In load_labeling(), each line must be a label in [0,5]\nBut found {}",current_label);
                return false;
            }
        }
        catch (const std::exception& e) { // if stoul() failed
            fmt::println(Logger::err("I/O"),"In load_labeling(), each line must be an unsigned integer\nBut found '{}'\nException message : {}",current_line,e.what());
            return false;
        }
        if(current_line_number >= facets_number) {
            fmt::println(Logger::err("I/O"),"In load_labeling(), the number of labels is greater than the number of facets");
            fmt::println(Logger::err("I/O"),"Number of labels so far = {}",current_line_number+1);
            fmt::println(Logger::err("I/O"),"Number of facets = {}",facets_number);
            return false;
        }
        label[current_line_number] = (index_t) current_label;
        current_line_number++;
    }

    ifs.close();

    //compare with expected size
    if (current_line_number != facets_number){
        fmt::println(Logger::err("I/O"),"In load_labeling(), the number of labels is lesser than the number of facets");
        fmt::println(Logger::err("I/O"),"Number of labels = {}",current_line_number+1);
        fmt::println(Logger::err("I/O"),"Number of facets = {}",facets_number);
        return false;
    }

    return true;
}

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