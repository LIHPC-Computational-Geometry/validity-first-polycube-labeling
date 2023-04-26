#include <geogram/basic/attributes.h>   // for GEO::Attribute
#include <geogram/basic/geometry.h>     // for GEO::vec3
#include <geogram/mesh/mesh_geometry.h> // for Geom::mesh_facet_normal()

#include <array>     // for std::array
#include <algorithm> // for std::max_element()
#include <iterator>  // for std::distance()

#include "labeling.h"

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name) {

    //open the file
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        Logger::out("I/O") << "Could not open labeling file '" << filename << "'" << std::endl;
        return false;
    }

    Logger::out("I/O") << "Loading file " << filename << "..." << std::endl;

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
                Logger::err("I/O") << "In load_labeling(), each line must be a label in [0,5]\n"
                                    << "But found " << current_label << std::endl;
                label.unbind();
                mesh.facets.attributes().delete_attribute_store(attribute_name);
                return false;
            }
        }
        catch (const std::exception& e) { // if stoul() failed
            Logger::err("I/O") << "In load_labeling(), each line must be an unsigned integer\n"
                                << "But found '" << current_line << "'\n"
                                << "Exception message : " << e.what() << std::endl;
            label.unbind();
            mesh.facets.attributes().delete_attribute_store(attribute_name);
            return false;
        }
        if(current_line_number >= facets_number) {
            Logger::err("I/O") << "In load_labeling(), the number of labels is greater than the number of facets\n"
                                << "Number of labels so far = " << current_line_number+1 << "\n"
                                << "Number of facets = " << facets_number << "\n";
            label.unbind();
            mesh.facets.attributes().delete_attribute_store(attribute_name);
            Logger::err("I/O") << "Labeling removed" << std::endl;
            return false;
        }
        label[current_line_number] = (index_t) current_label;
        current_line_number++;
    }

    //compare with expected size
    if (current_line_number != facets_number){
        Logger::err("I/O") << "In load_labeling(), the number of labels is lesser than the number of facets\n"
                            << "Number of labels = " << current_line_number+1 << "\n"
                            << "Number of facets = " << facets_number << "\n";
        label.unbind();
        mesh.facets.attributes().delete_attribute_store(attribute_name);
        Logger::err("I/O") << "Labeling removed" << std::endl;
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