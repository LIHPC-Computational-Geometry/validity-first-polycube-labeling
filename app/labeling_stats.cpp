// Usage :
//
//   ./bin/labeling_stats input_mesh.obj labeling.txt
//     -> write to stdout
//
//   ./bin/labeling_stats input_mesh.obj labeling.txt output.json
//     -> write to output.json

#include <geogram/mesh/mesh.h>              // for Mesh
#include <geogram/mesh/mesh_io.h>           // for mesh_load(), mesh_save()
#include <geogram/basic/file_system.h>      // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <nlohmann/json.hpp>

#include <iostream>
#include <iomanip>
#include <vector>

#include "labeling.h"
#include "labeling_io.h"
#include "labeling_graph.h"
#include "stats.h"
#include "geometry_halfedges.h"

#define LABELING_ATTRIBUTE_NAME "label"

using namespace GEO;

int main(int argc, char** argv) {
    
    GEO::initialize();

    CmdLine::declare_arg("allow-opposite-labels",true,"Allow boundaries between opposite labels");

    std::vector<std::string> filenames;
    if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh input_surface_labeling <output_JSON>" // third filename is facultative
		))
	{
        fmt::println(Logger::err("I/O"),"Usage should be\n{} input_mesh input_surface_labeling",argv[0]); Logger::err("I/O").flush();
		return 1;
	}

    if(!FileSystem::is_file(filenames[0])) {
        fmt::println(Logger::err("I/O"),"{} does not exist",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }

    if(!FileSystem::is_file(filenames[1])) {
        fmt::println(Logger::err("I/O"),"{} does not exist",filenames[1]); Logger::err("I/O").flush();
        return 1;
    }

    Mesh input_mesh;
    if(!mesh_load(filenames[0],input_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to load the mesh from {}",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }

    MeshExt mesh_ext(input_mesh); // compute normals, vertex-facets adjacency, feature edges

    Attribute<index_t> labeling(input_mesh.facets.attributes(),LABELING_ATTRIBUTE_NAME);
    if (!load_labeling(filenames[1],input_mesh,labeling)) {
        exit(1);
    }

    LabelingGraph lg;
    lg.fill_from(input_mesh,labeling,CmdLine::get_arg_bool("allow-opposite-labels"));
    mesh_ext.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);

    nlohmann::json output_JSON;

    ////////////////////////////////
    // Parameter
    ////////////////////////////////
    
    output_JSON["is_allowing_boundaries_between_opposite_labels"] = lg.is_allowing_boundaries_between_opposite_labels();

    ////////////////////////////////
    // Charts
    ////////////////////////////////

    output_JSON["charts"]["nb"] = lg.nb_charts();
    output_JSON["charts"]["invalid"] = lg.nb_invalid_charts();

    ////////////////////////////////
    // Boundaries
    ////////////////////////////////

    output_JSON["boundaries"]["nb"] = lg.nb_boundaries();
    output_JSON["boundaries"]["invalid"] = lg.nb_invalid_boundaries();
    output_JSON["boundaries"]["non-monotone"] = lg.nb_non_monotone_boundaries();

    ////////////////////////////////
    // Corners
    ////////////////////////////////

    output_JSON["corners"]["nb"] = lg.nb_corners();
    output_JSON["corners"]["invalid"] = lg.nb_invalid_corners();

    ////////////////////////////////
    // Turning-points
    ////////////////////////////////

    output_JSON["turning-points"]["nb"] = lg.nb_turning_points();

    ////////////////////////////////
    // Fidelity
    ////////////////////////////////

    // compute normals
    std::vector<vec3> normals(input_mesh.facets.nb());
    FOR(f,input_mesh.facets.nb()) {
        normals[f] = normalize(Geom::mesh_facet_normal(input_mesh,f));
    }
    // compute fidelity
    IncrementalStats fidelity_stats;
    compute_per_facet_fidelity(mesh_ext, labeling,"fidelity",fidelity_stats);
    // fill JSON
    output_JSON["fidelity"]["min"] = fidelity_stats.min();
    output_JSON["fidelity"]["max"] = fidelity_stats.max();
    output_JSON["fidelity"]["avg"] = fidelity_stats.avg();
    output_JSON["fidelity"]["sd"] = fidelity_stats.sd();

    ///////////////////////////////
    // Feature edges preservation
    ///////////////////////////////

    unsigned int nb_lost = count_lost_feature_edges(mesh_ext);
    output_JSON["feature-edges"]["removed"] = input_mesh.edges.nb() - mesh_ext.feature_edges.nb();
    output_JSON["feature-edges"]["lost"] = nb_lost;
    output_JSON["feature-edges"]["preserved"] = mesh_ext.feature_edges.nb() - nb_lost;

    // write to stdout or file, according to the number of arguments

    if (filenames.size() <= 2) {
        // no output JSON filename provided
        std::cout << std::setw(4) << output_JSON << std::endl;
    }
    else {
        // write values in a JSON file
        std::fstream ofs(filenames[2],std::ios_base::out);
        if(ofs.good()) {
            fmt::println(Logger::out("I/O"),"Saving file {}...",filenames[2]); Logger::out("I/O").flush();
            ofs << std::setw(4) << output_JSON << std::endl;
        }
        else {
            fmt::println(Logger::err("I/O"),"Cannot write into {}",filenames[2]); Logger::err("I/O").flush();
        }
    }

    return 0;
}