#include <geogram/mesh/mesh.h>          // for Mesh, cell types
#include <geogram/mesh/mesh_io.h>       // for mesh_load()
#include <geogram/mesh/mesh_geometry.h> // for mesh_facet_area(), mesh_cell_volume()
#include <geogram/basic/file_system.h>  // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>         // for length()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <nlohmann/json.hpp>

#include <iostream>
#include <iomanip>
#include <array>

#include "labeling.h"
#include "basic_stats.h"

#define CELL_TYPE_KEY_STR(cell_type)  ( ((cell_type) == MESH_TET)       ? "tetrahedra"  : (\
                                        ((cell_type) == MESH_HEX)       ? "hexahedra"   : (\
                                        ((cell_type) == MESH_PRISM)     ? "prisms"      : (\
                                        ((cell_type) == MESH_PYRAMID)   ? "pyramids"    : (\
                                        ((cell_type) == MESH_CONNECTOR) ? "connectors"  : (\
                                        "invalid cell type" ))))))

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

    ////////////////////////////////
    // Vertices
    ////////////////////////////////
    
    BasicStats dim_x_stats,
               dim_y_stats,
               dim_z_stats;
    output_JSON["vertices"]["nb"] = input_mesh.vertices.nb();
    if(input_mesh.vertices.nb()!=0) {
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
    }

    ////////////////////////////////
    // Edges
    ////////////////////////////////

    output_JSON["edges"]["nb"] = input_mesh.edges.nb();
    if(input_mesh.edges.nb()!=0) {
        BasicStats edges_length_stats;
        FOR(e,input_mesh.edges.nb()) {
            edges_length_stats.insert(length(
                input_mesh.vertices.point(input_mesh.edges.vertex(e,1)) -
                input_mesh.vertices.point(input_mesh.edges.vertex(e,0))
            ));
        }
        output_JSON["edges"]["length"]["min"] = edges_length_stats.min();
        output_JSON["edges"]["length"]["max"] = edges_length_stats.max();
        output_JSON["edges"]["length"]["avg"] = edges_length_stats.avg();
        output_JSON["edges"]["length"]["sd"]  = edges_length_stats.sd();
    }

    ////////////////////////////////
    // Facets
    ////////////////////////////////

    output_JSON["facets"]["nb"] = input_mesh.facets.nb();
    if(input_mesh.facets.nb()!=0) {
        BasicStats facets_area_stats;
        FOR(f,input_mesh.facets.nb()) {
            facets_area_stats.insert(mesh_facet_area(input_mesh,f));
        }
        output_JSON["facets"]["area"]["min"] = facets_area_stats.min();
        output_JSON["facets"]["area"]["max"] = facets_area_stats.max();
        output_JSON["facets"]["area"]["avg"] = facets_area_stats.avg();
        output_JSON["facets"]["area"]["sd"]  = facets_area_stats.sd();
        output_JSON["facets"]["area"]["sum"] = facets_area_stats.sum();
    }

    ////////////////////////////////
    // Cells
    ////////////////////////////////

    output_JSON["cells"]["nb"] = input_mesh.cells.nb();
    if(input_mesh.cells.nb()!=0) {
        std::array<unsigned int,MESH_NB_CELL_TYPES> per_type_count;
        per_type_count.fill(0);
        BasicStats cells_volume_stats;
        FOR(c,input_mesh.cells.nb()) {
            cells_volume_stats.insert(mesh_cell_volume(input_mesh,c));
            per_type_count[input_mesh.cells.type(c)]++;
        }
        output_JSON["cells"]["volume"]["min"] = cells_volume_stats.min();
        output_JSON["cells"]["volume"]["max"] = cells_volume_stats.max();
        output_JSON["cells"]["volume"]["avg"] = cells_volume_stats.avg();
        output_JSON["cells"]["volume"]["sd"]  = cells_volume_stats.sd();
        output_JSON["cells"]["volume"]["sum"] = cells_volume_stats.sum();
        FOR(cell_type,MESH_NB_CELL_TYPES) {
            if(per_type_count[cell_type]!=0) {
                output_JSON["cells"]["by_type"][CELL_TYPE_KEY_STR(cell_type)] = per_type_count[cell_type];
            }
        }
    }

    std::cout << std::setw(4) << output_JSON << std::endl;

    return 0;
}