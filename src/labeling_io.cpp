#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include <string>

#include "labeling_io.h"
#include "labeling.h"

bool load_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name) {

    //open the file
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        fmt::println(Logger::out("I/O"),"Could not open labeling file '{}'",filename); 
        return false;
    }

    fmt::println(Logger::out("I/O"),"Loading file {}...",filename); Logger::out("I/O").flush();

    std::string current_line;
    unsigned long current_label; // necessary type for string to unsigned int conversion (stoul)
    index_t current_line_number = 0;
    index_t facets_number = mesh.facets.nb(); // expected line number (one line per facet)

    // Add an attribute on facets. see https://github.com/BrunoLevy/geogram/wiki/Mesh#attributes
    // The type must be Attribute<index_t> to be used with geogram/mesh/mesh_halfedges.h . See MeshHalfedges::facet_region_
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name);

    //fill the attribute
    while (ifs >> current_line) {
        try {
            current_label = std::stoul(current_line);
            if(current_label > 5) {
                fmt::println(Logger::err("I/O"),"In load_labeling(), each line must be a label in [0,5]\nBut found {}",current_label); Logger::err("I/O").flush();
                return false;
            }
        }
        catch (const std::exception& e) { // if stoul() failed
            fmt::println(Logger::err("I/O"),"In load_labeling(), each line must be an unsigned integer\nBut found '{}'\nException message : {}",current_line,e.what()); Logger::err("I/O").flush();
            return false;
        }
        if(current_line_number >= facets_number) {
            fmt::println(Logger::err("I/O"),"In load_labeling(), the number of labels is greater than the number of facets");
            fmt::println(Logger::err("I/O"),"Number of labels so far = {}",current_line_number+1);
            fmt::println(Logger::err("I/O"),"Number of facets = {}",facets_number); Logger::err("I/O").flush();
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
        fmt::println(Logger::err("I/O"),"Number of facets = {}",facets_number); Logger::err("I/O").flush();
        return false;
    }

    return true;
}

bool save_labeling(const std::string& filename, Mesh& mesh, const char* attribute_name) {
    Attribute<index_t> label(mesh.facets.attributes(), attribute_name); // fetch the labeling
    auto out = fmt::output_file(filename);
    FOR(f,mesh.facets.nb()) {
        out.print("{}",label[f]);
        if(f != mesh.facets.nb()-1)
            out.print("\n");
    }
    return true;
}