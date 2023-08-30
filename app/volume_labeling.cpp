// Rewriting of https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/flagging_utils.cpp#L111
// using Geogram instead of libigl.
// 
// Given a surface labeling (one label per surface triangle) and the surface map (on which cell facets are surface triangles),
// write the volume labeling (on label per cell facet).

#include <geogram/basic/file_system.h>  // for FileSystem::is_file()
#include <geogram/basic/logger.h>       // for Logger::out() and Logger::err()

#include <fmt/core.h>
#include <fmt/os.h>         // for fmt::output_file()
#include <fmt/ostream.h>    // to use fmt::print() on ostreams

#include <vector>   // for std::vector
#include <fstream>  // for std::ifstream
#include <string>   // for std::stoi(), std::getline()

using namespace GEO;

int main(int argc, const char** argv) {

    if (argc<4) {
        fmt::println("Missing arguments");
        fmt::println("./volume_labeling surface_labeling.txt surface_map.txt tetra_labeling.txt");
        return 1;
    }

    GEO::initialize();

    geo_assert(FileSystem::is_file(argv[1]));
    geo_assert(FileSystem::is_file(argv[2]));

    std::string line;
    std::size_t separator_pos = std::string::npos;
    index_t nb_facets = index_t(-1);
    index_t nb_cells = index_t(-1);
    index_t line_count = 0;

    // first, read the surface map (contains nb facets and nb cells)

    std::ifstream ifs(argv[2]);
	if (!ifs.is_open()) {
        fmt::println(Logger::err("I/O"),"Unable to open {}",argv[2]); Logger::err("I/O").flush();
		geo_abort();
	}
    fmt::println(Logger::out("I/O"),"Reading {}...",argv[2]); Logger::out("I/O").flush();
    // read first line of header
    std::getline(ifs, line);
    separator_pos = line.find(" ");
    if(separator_pos != std::string::npos) { // if line contains ' '
        line = line.substr(0,separator_pos); // remove what is after ' '
    }
    nb_facets = (index_t) std::stoi(line); // integer conversion
    fmt::println(Logger::out("I/O"),"{} triangles",nb_facets); Logger::out("I/O").flush();
    // read second line of header
    std::getline(ifs, line);
    separator_pos = line.find(" ");
    if(separator_pos != std::string::npos) { // if line contains ' '
        line = line.substr(0,separator_pos); // remove what is after ' '
    }
    nb_cells = (index_t) std::stoi(line); // integer conversion
    fmt::println(Logger::out("I/O"),"{} tetrahedra",nb_cells); Logger::out("I/O").flush();
    std::vector<index_t> surface_map(nb_facets);
    line_count = 2;
	FOR(f, nb_facets) {
		if (ifs.eof()) break;
        std::getline(ifs, line);
		surface_map[f] = (index_t) std::stoi(line);
        line_count++;
	}
	ifs.close();
    geo_assert(line_count == nb_facets+2); // check if the file is long enough

    // read the surface labeling

    std::vector<index_t> surface_labeling(nb_facets);
    ifs.open(argv[1]);
    if (!ifs.is_open()) {
        fmt::println(Logger::err("I/O"),"Unable to open {}",argv[1]); Logger::err("I/O").flush();
		geo_abort();
	}
    fmt::println(Logger::out("I/O"),"Reading {}...",argv[1]); Logger::out("I/O").flush();
    line_count = 0;
    FOR(f, nb_facets) {
		if (ifs.eof()) break;
        ifs >> line;
		surface_labeling[f] = (index_t) std::stoi(line);
        line_count++;
	}
    ifs.close();
    geo_assert(line_count == nb_facets); // check if the file is long enough

    // fill volume labeling

    std::vector<signed_index_t> volume_labeling(4 * nb_cells,-1);
    FOR(f, nb_facets) {
        volume_labeling[surface_map[f]] = (signed_index_t) surface_labeling[f];
	}

    // write volume labeling

    fmt::println(Logger::out("I/O"),"Writing {}...",argv[3]); Logger::out("I/O").flush();
    auto out = fmt::output_file(argv[3]);
    FOR(cf, 4 * nb_cells) { // for each cell facet
        out.print("{}\n",volume_labeling[cf]);
	}
    out.close();    

    return 0;
}