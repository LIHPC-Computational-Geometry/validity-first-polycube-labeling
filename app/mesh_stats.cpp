// Usage :
//
//   ./bin/mesh_stats input.mesh
//     -> write to stdout
//
//   ./bin/mesh_stats input.mesh output.json
//     -> write to output.json
//
//   ./bin/mesh_stats input.mesh principal_axes=principal_axes.obj
//     -> export a mesh with the center and the 3 vectors (as edges) of the principal axes

#include <geogram/mesh/mesh.h>              // for Mesh, cell types
#include <geogram/mesh/mesh_io.h>           // for mesh_load(), mesh_save()
#include <geogram/mesh/mesh_geometry.h>     // for mesh_facet_area(), mesh_cell_volume()
#include <geogram/basic/file_system.h>      // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>             // for length()
#include <geogram/points/principal_axes.h>  // for PrincipalAxes3d

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <nlohmann/json.hpp>

#include <iostream>
#include <iomanip>
#include <array>
#include <set>
#include <vector>

#include "labeling.h"
#include "stats.h"
#include "geometry_hexahedra.h" // for compute_scaled_jacobian()
#include "geometry.h" // for facet_normals_are_inward()

#define CELL_TYPE_KEY_STR(cell_type)  ( ((cell_type) == MESH_TET)       ? "tetrahedra"  : (\
                                        ((cell_type) == MESH_HEX)       ? "hexahedra"   : (\
                                        ((cell_type) == MESH_PRISM)     ? "prisms"      : (\
                                        ((cell_type) == MESH_PYRAMID)   ? "pyramids"    : (\
                                        ((cell_type) == MESH_CONNECTOR) ? "connectors"  : (\
                                        "invalid cell type" ))))))

#define VERTEX_ORIGIN       0
#define VERTEX_TIP_1ST_AXIS 1
#define VERTEX_TIP_2ND_AXIS 2
#define VERTEX_TIP_3RD_AXIS 3

using namespace GEO;

int main(int argc, char** argv) {
    
    GEO::initialize();

    CmdLine::declare_arg(
        "principal_axes",
        "",
        "path to the output principal axes mesh"
    );

    std::vector<std::string> filenames;
    if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh <output_JSON>" // second filename is facultative
		))
	{
		return 1;
	}

    std::string principal_axes_filename = GEO::CmdLine::get_arg("principal_axes");

    if(!FileSystem::is_file(filenames[0])) {
        fmt::println(Logger::err("I/O"),"{} does not exist",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }

    mesh_io_initialize();
    Mesh input_mesh;
    if(!mesh_load(filenames[0],input_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to load the mesh from {}",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }

    nlohmann::json output_JSON;

    ////////////////////////////////
    // Vertices
    ////////////////////////////////
    
    output_JSON["vertices"]["nb"] = input_mesh.vertices.nb();
    if(input_mesh.vertices.nb()!=0) {
        IncrementalStats dim_x_stats,
                         dim_y_stats,
                         dim_z_stats;
        PrincipalAxes3d principal_axes;
        principal_axes.begin();
        vec3 coordinates;
        FOR(v,input_mesh.vertices.nb()) {
            coordinates = input_mesh.vertices.point(v);
            dim_x_stats.insert(coordinates.x);
            dim_y_stats.insert(coordinates.y);
            dim_z_stats.insert(coordinates.z);
            principal_axes.add_point(coordinates);
        }
        principal_axes.end();
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
        // principal axes, as a list (3 items) of objects (2 keys for each)
        output_JSON["vertices"]["principal_axes"] = {
            { nlohmann::json::object({
                {"axis", {principal_axes.axis(0)[0], principal_axes.axis(0)[1], principal_axes.axis(0)[2]} },
                {"eigenvalue", principal_axes.eigen_value(0) }})
            },
            { nlohmann::json::object({
                {"axis", {principal_axes.axis(1)[0], principal_axes.axis(1)[1], principal_axes.axis(1)[2]} },
                {"eigenvalue", principal_axes.eigen_value(1) }})
            },
            { nlohmann::json::object({
                {"axis", {principal_axes.axis(2)[0], principal_axes.axis(2)[1], principal_axes.axis(2)[2]} },
                {"eigenvalue", principal_axes.eigen_value(2) }})
            }
        };
        if(!principal_axes_filename.empty()) {
            Mesh principal_axes_mesh;
            principal_axes_mesh.vertices.create_vertices(4); // {center, first axis tip, second axis tip, third axis tip}
            principal_axes_mesh.vertices.point(VERTEX_ORIGIN) = principal_axes.center();
            principal_axes_mesh.vertices.point(VERTEX_TIP_1ST_AXIS) = principal_axes.center() + principal_axes.axis(0)*principal_axes.eigen_value(0);
            principal_axes_mesh.vertices.point(VERTEX_TIP_2ND_AXIS) = principal_axes.center() + principal_axes.axis(1)*principal_axes.eigen_value(1);
            principal_axes_mesh.vertices.point(VERTEX_TIP_3RD_AXIS) = principal_axes.center() + principal_axes.axis(2)*principal_axes.eigen_value(2);
            principal_axes_mesh.edges.create_edges(3);
            // edge n°0 : 1st principal axis
            principal_axes_mesh.edges.set_vertex(0,0,VERTEX_ORIGIN);
            principal_axes_mesh.edges.set_vertex(0,1,VERTEX_TIP_1ST_AXIS);
            // edge n°1 : 2nd principal axis
            principal_axes_mesh.edges.set_vertex(1,0,VERTEX_ORIGIN);
            principal_axes_mesh.edges.set_vertex(1,1,VERTEX_TIP_2ND_AXIS);
            // edge n°2 : 3rd principal axis
            principal_axes_mesh.edges.set_vertex(2,0,VERTEX_ORIGIN);
            principal_axes_mesh.edges.set_vertex(2,1,VERTEX_TIP_3RD_AXIS);
            mesh_save(principal_axes_mesh,principal_axes_filename);
        }
    }

    ////////////////////////////////
    // Edges
    ////////////////////////////////

    // In Geogram, edges are unrelated to the surfacic and volumetric parts of the mesh
    // https://github.com/BrunoLevy/geogram/wiki/Mesh#mesh-edges

    output_JSON["edges"]["nb"] = input_mesh.edges.nb();
    if(input_mesh.edges.nb()!=0) {
        IncrementalStats edges_length_stats;
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
        IncrementalStats facets_area_stats;
        FOR(f,input_mesh.facets.nb()) {
            facets_area_stats.insert(mesh_facet_area(input_mesh,f));
        }
        output_JSON["facets"]["area"]["min"] = facets_area_stats.min();
        output_JSON["facets"]["area"]["max"] = facets_area_stats.max();
        output_JSON["facets"]["area"]["avg"] = facets_area_stats.avg();
        output_JSON["facets"]["area"]["sd"]  = facets_area_stats.sd();
        output_JSON["facets"]["area"]["sum"] = facets_area_stats.sum();
        if(input_mesh.cells.nb()==0) { // if there is no cells in the mesh (surface only)
            output_JSON["facets"]["normals_outward"] = (facet_normals_are_inward(input_mesh) == false);
        }
    }

    ////////////////////////////////
    // Facet edges
    ////////////////////////////////

    std::set<std::set<index_t>> unique_facet_edges;
    index_t nb_corners = index_t(-1);
    index_t v0 = index_t(-1);
    index_t v1 = index_t(-1);
    if(input_mesh.facets.nb()!=0) {
        FOR(f,input_mesh.facets.nb()) { // for each facet
            nb_corners = input_mesh.facets.nb_corners(f);
            FOR(lv,nb_corners) { // for each local vertex of f
                v0 = input_mesh.facet_corners.vertex(input_mesh.facets.corner(f,lv));
                v1 = input_mesh.facet_corners.vertex(input_mesh.facets.corner(f,(lv+1) % nb_corners));
                unique_facet_edges.insert({v0,v1}); // auto-sorted by the set
            }
        }
        output_JSON["facets"]["edges"]["nb"] = unique_facet_edges.size();
        IncrementalStats facet_edges_length_stats;
        for(const std::set<index_t>& edge : unique_facet_edges) {
            geo_debug_assert(edge.size()==2); // only 2 vertex indices
            facet_edges_length_stats.insert(length(
                mesh_vertex(input_mesh,*edge.rbegin()) - mesh_vertex(input_mesh,*edge.begin())
            ));
        }
        output_JSON["facets"]["edges"]["length"]["min"] = facet_edges_length_stats.min();
        output_JSON["facets"]["edges"]["length"]["max"] = facet_edges_length_stats.max();
        output_JSON["facets"]["edges"]["length"]["avg"] = facet_edges_length_stats.avg();
        output_JSON["facets"]["edges"]["length"]["sd"] = facet_edges_length_stats.sd();
        fmt::println(
            "Euler characteristic for surface polyhedra = V - E + F = {} - {} + {} = {}",
            input_mesh.vertices.nb(),
            unique_facet_edges.size(),
            input_mesh.facets.nb(),
            input_mesh.vertices.nb() - unique_facet_edges.size() + input_mesh.facets.nb()
        );
    }

    ////////////////////////////////
    // Cells
    ////////////////////////////////

    output_JSON["cells"]["nb"] = input_mesh.cells.nb();
    if(input_mesh.cells.nb()!=0) {
        std::array<unsigned int,MESH_NB_CELL_TYPES> per_type_count;
        per_type_count.fill(0);
        IncrementalStats cells_volume_stats;
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
        if(per_type_count[MESH_HEX]!=0) {
            IncrementalStats scaled_jacobian_stats;
            compute_scaled_jacobian(input_mesh,scaled_jacobian_stats);
            output_JSON["cells"]["quality"]["hex_SJ"]["min"] = scaled_jacobian_stats.min();
            output_JSON["cells"]["quality"]["hex_SJ"]["max"] = scaled_jacobian_stats.max();
            output_JSON["cells"]["quality"]["hex_SJ"]["avg"] = scaled_jacobian_stats.avg();
            output_JSON["cells"]["quality"]["hex_SJ"]["sd"]  = scaled_jacobian_stats.sd();
        }
    }

    ////////////////////////////////
    // Cells facets
    ////////////////////////////////

    std::map<std::set<index_t>,std::vector<index_t>> unique_cell_facets; // set as key -> ensure uniqueness, vector as value -> keep vertices ordering for later
    std::vector<index_t> facet_vertices;
    if(input_mesh.cells.nb()!=0) {
        FOR(c,input_mesh.cells.nb()) { // for each cell
            index_t nb_facets = input_mesh.cells.nb_facets(c);
            FOR(lf,nb_facets) { // for each local facet of c
                index_t nb_vertices = input_mesh.cells.facet_nb_vertices(c,lf);
                facet_vertices.clear();
                FOR(lv,nb_vertices) { // for each local vertex of this facet
                    facet_vertices.push_back(input_mesh.cells.facet_vertex(c,lf,lv));
                }
                unique_cell_facets.insert_or_assign(
                    std::set(facet_vertices.begin(),facet_vertices.end()),// auto-sorted by the set
                    facet_vertices
                );
            }
        }
        output_JSON["cells"]["facets"]["nb"] = unique_cell_facets.size();
        // TODO compute min/max/avg/sd area like GEO::Geom::mesh_facet_area() does
    }


    ////////////////////////////////
    // Cells edges
    ////////////////////////////////

    std::set<std::set<index_t>> unique_cell_edges;
    if(input_mesh.cells.nb()!=0) {
        FOR(c,input_mesh.cells.nb()) { // for each cell
            FOR(ce,input_mesh.cells.nb_edges(c)) { // for each edge of c
                v0 = input_mesh.cells.edge_vertex(c,ce,0);
                v1 = input_mesh.cells.edge_vertex(c,ce,1);
                unique_cell_edges.insert({v0,v1}); // auto-sorted by the set
            }
        }
        output_JSON["cells"]["edges"]["nb"] = unique_cell_edges.size();
        IncrementalStats cell_edges_length_stats;
        for(const std::set<index_t>& edge : unique_cell_edges) {
            geo_debug_assert(edge.size()==2); // only 2 vertex indices
            cell_edges_length_stats.insert(length(
                mesh_vertex(input_mesh,*edge.rbegin()) - mesh_vertex(input_mesh,*edge.begin())
            ));
        }
        output_JSON["cells"]["edges"]["length"]["min"] = cell_edges_length_stats.min();
        output_JSON["cells"]["edges"]["length"]["max"] = cell_edges_length_stats.max();
        output_JSON["cells"]["edges"]["length"]["avg"] = cell_edges_length_stats.avg();
        output_JSON["cells"]["edges"]["length"]["sd"] = cell_edges_length_stats.sd();
        fmt::println(
            "Euler characteristic for volume polyhedra = V - E + F - C = {} - {} + {} - {} = {}",
            input_mesh.vertices.nb(),
            unique_cell_edges.size(),
            unique_cell_facets.size(),
            input_mesh.cells.nb(),
            input_mesh.vertices.nb() - unique_cell_edges.size() + unique_cell_facets.size() - input_mesh.cells.nb()
        );
    }

    // write to stdout or file, according to the number of arguments

    if(filenames.size() <= 1) {
        // no output JSON filename provided
        std::cout << std::setw(4) << output_JSON << std::endl;
    }
    else {
        // write values in a JSON file
        std::fstream ofs(filenames[1],std::ios_base::out);
        if(ofs.good()) {
            fmt::println(Logger::out("I/O"),"Saving file {}...",filenames[1]); Logger::out("I/O").flush();
            ofs << std::setw(4) << output_JSON << std::endl;
        }
        else {
            fmt::println(Logger::err("I/O"),"Cannot write into {}",filenames[1]); Logger::err("I/O").flush();
        }
    }

    return 0;
}