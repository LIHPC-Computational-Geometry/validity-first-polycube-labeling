#include <geogram/mesh/mesh.h>          // for Mesh
#include <geogram/mesh/mesh_io.h>       // for mesh_load()
#include <geogram/basic/file_system.h>  // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <nlohmann/json.hpp>

#include <iostream>
#include <iomanip>

#include "labeling.h"
#include "basic_stats.h"

using namespace GEO;

int main(int argc, const char** argv) {

    if (argc<1) {
        fmt::println("Missing arguments");
        fmt::println("./volume_labeling input_mesh");
        return 1;
    }
    FileSystem::initialize();
    Logger::initialize();
    Logger::instance()->set_minimal(true);
    CmdLine::initialize();
    mesh_io_initialize();

    if(!FileSystem::is_file(argv[1])) {
        fmt::println(Logger::err("I/O"),"{} does not exist",argv[1]); Logger::err("I/O").flush();
        return 1;
    }

    Mesh input_mesh;
    if(!mesh_load(argv[1],input_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to load the mesh from {}",argv[1]); Logger::err("I/O").flush();
        return 1;
    }

    nlohmann::json output_JSON;
    
    BasicStats dim_x_stats,
               dim_y_stats,
               dim_z_stats;
    output_JSON["vertices"]["nb"] = input_mesh.vertices.nb();
    vec3 coordinates;
    FOR(v,input_mesh.vertices.nb()) {
        coordinates = input_mesh.vertices.point(v);
        dim_x_stats.insert(coordinates.x);
        dim_y_stats.insert(coordinates.y);
        dim_z_stats.insert(coordinates.z);
    }

    // x axis
    output_JSON["vertices"]["x"]["min"] = dim_x_stats.min();
    output_JSON["vertices"]["x"]["max"] = dim_x_stats.max();
    output_JSON["vertices"]["x"]["avg"] = dim_x_stats.avg();
    output_JSON["vertices"]["x"]["sd"]  = dim_x_stats.sd();
    // y axis
    output_JSON["vertices"]["y"]["min"] = dim_y_stats.min();
    output_JSON["vertices"]["y"]["max"] = dim_y_stats.max();
    output_JSON["vertices"]["y"]["avg"] = dim_y_stats.avg();
    output_JSON["vertices"]["y"]["sd"]  = dim_y_stats.sd();
    // z axis
    output_JSON["vertices"]["z"]["min"] = dim_z_stats.min();
    output_JSON["vertices"]["z"]["max"] = dim_z_stats.max();
    output_JSON["vertices"]["z"]["avg"] = dim_z_stats.avg();
    output_JSON["vertices"]["z"]["sd"]  = dim_z_stats.sd();

    std::cout << std::setw(4) << output_JSON << std::endl;

    return 0;
}