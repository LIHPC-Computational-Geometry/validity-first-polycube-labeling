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
#include "LabelingGraph.h"
#include "basic_stats.h"
#include "CustomMeshHalfedges.h"

#define LABELING_ATTRIBUTE_NAME "label"

using namespace GEO;

int main(int argc, char** argv) {
    
    // inside of the GEO::initialize() function, modified to have the logger in the "minimal" mode
    Environment* env = Environment::instance();
    env->set_value("version", "inaccessible"); // some code after expects the "version" environment variable to exist
    env->set_value("release_date", "inaccessible"); // idem
    env->set_value("SVN revision", "inaccessible"); // idem
    FileSystem::initialize();
    Logger::initialize();
    Logger::instance()->set_minimal(true);
    CmdLine::initialize();
    CmdLine::import_arg_group("sys"); // declares sys:compression_level, needed by mesh_save() for .geogram files

    CmdLine::declare_arg("allow-opposite-labels",true,"Allow boundaries between opposite labels");

    std::vector<std::string> filenames;
    if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh input_surface_labeling" // second filename is facultative
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

    mesh_io_initialize();
    Mesh input_mesh;
    if(!mesh_load(filenames[0],input_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to load the mesh from {}",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }

    // remove feature edges on edge with small angle
    unsigned int init_nb_feature_edges = input_mesh.edges.nb();
	std::vector<std::vector<index_t>> adj_facets; // for each vertex, store adjacent facets. no ordering
	remove_feature_edges_with_low_dihedral_angle(input_mesh,adj_facets);
	// store them as a set
	std::set<std::pair<index_t,index_t>> feature_edges;
	transfer_feature_edges(input_mesh,feature_edges);

    if (!load_labeling(filenames[1],input_mesh,LABELING_ATTRIBUTE_NAME)) {
        exit(1);
    }

    StaticLabelingGraph slg;
    slg.fill_from(input_mesh,LABELING_ATTRIBUTE_NAME,CmdLine::get_arg_bool("allow-opposite-labels"),feature_edges);
    CustomMeshHalfedges mesh_he(input_mesh);
    mesh_he.set_use_facet_region(LABELING_ATTRIBUTE_NAME);

    nlohmann::json output_JSON;

    ////////////////////////////////
    // Parameter
    ////////////////////////////////
    
    output_JSON["is_allowing_boundaries_between_opposite_labels"] = slg.is_allowing_boundaries_between_opposite_labels();

    ////////////////////////////////
    // Charts
    ////////////////////////////////

    output_JSON["charts"]["nb"] = slg.nb_charts();
    output_JSON["charts"]["invalid"] = slg.nb_invalid_charts();

    ////////////////////////////////
    // Boundaries
    ////////////////////////////////

    output_JSON["boundaries"]["nb"] = slg.nb_boundaries();
    output_JSON["boundaries"]["invalid"] = slg.nb_invalid_boundaries();
    output_JSON["boundaries"]["non-monotone"] = slg.nb_non_monotone_boundaries();

    ////////////////////////////////
    // Corners
    ////////////////////////////////

    output_JSON["corners"]["nb"] = slg.nb_corners();
    output_JSON["corners"]["invalid"] = slg.nb_invalid_corners();

    ////////////////////////////////
    // Turning-points
    ////////////////////////////////

    output_JSON["turning-points"]["nb"] = slg.nb_turning_points();

    ////////////////////////////////
    // Fidelity
    ////////////////////////////////

    // compute normals
    std::vector<vec3> normals(input_mesh.facets.nb());
    FOR(f,input_mesh.facets.nb()) {
        normals[f] = normalize(Geom::mesh_facet_normal(input_mesh,f));
    }
    // compute fidelity
    BasicStats fidelity_stats;
    compute_per_facet_fidelity(input_mesh, normals, LABELING_ATTRIBUTE_NAME,"fidelity",fidelity_stats);
    // fill JSON
    output_JSON["fidelity"]["min"] = fidelity_stats.min();
    output_JSON["fidelity"]["max"] = fidelity_stats.max();
    output_JSON["fidelity"]["avg"] = fidelity_stats.avg();
    output_JSON["fidelity"]["sd"] = fidelity_stats.sd();

    ///////////////////////////////
    // Feature edges preservation
    ///////////////////////////////

    unsigned int nb_lost = count_lost_feature_edges(mesh_he,feature_edges);
    output_JSON["feature-edges"]["removed"] = init_nb_feature_edges - feature_edges.size();
    output_JSON["feature-edges"]["lost"] = nb_lost;
    output_JSON["feature-edges"]["preserved"] = feature_edges.size() - nb_lost;

    std::cout << std::setw(4) << output_JSON << std::endl;

    return 0;
}