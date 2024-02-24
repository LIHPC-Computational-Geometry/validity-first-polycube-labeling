// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <geogram/basic/command_line.h>			// for CmdLine::parse()
#include <geogram/basic/command_line_args.h>    // for CmdLine::import_arg_group()
#include <geogram/basic/logger.h>               // for GEO::Logger::*
#include <geogram/mesh/mesh.h>                  // for Mesh
#include <geogram/mesh/mesh_io.h>               // for mesh_load()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <vector>

#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges
#include "glTF.h"
#include "labeling.h"

using GEO::index_t; // to use the FOR() macro of Geogram

int main(int argc, char **argv)
{
    GEO::initialize();
    GEO::CmdLine::import_arg_group("sys"); // declares sys:compression_level, needed by mesh_save() for .geogram files

    std::vector<std::string> filenames;
    if(!CmdLine::parse(argc,argv,filenames,"<input_mesh>"))	{
		return 1;
	}

    GEO::Mesh M;
    if(!filenames.empty()) { // if <input_mesh> is provided
        mesh_load(filenames[0],M);
    }
    else {
        filenames.push_back("box"); // will be use as output filename
        // create a simple mesh (a box)
        M.vertices.create_vertices(8);
        M.vertices.point(0) = GEO::vec3(-0.5, -0.5,  0.5);
        M.vertices.point(1) = GEO::vec3( 0.5, -0.5,  0.5);
        M.vertices.point(2) = GEO::vec3(-0.5,  0.5,  0.5);
        M.vertices.point(3) = GEO::vec3( 0.5,  0.5,  0.5);
        M.vertices.point(4) = GEO::vec3( 0.5, -0.5, -0.5);
        M.vertices.point(5) = GEO::vec3(-0.5, -0.5, -0.5);
        M.vertices.point(6) = GEO::vec3( 0.5,  0.5, -0.5);
        M.vertices.point(7) = GEO::vec3(-0.5,  0.5, -0.5);
        M.facets.create_triangles(12);
        // triangle 0
        M.facets.set_vertex(0,0,0);
        M.facets.set_vertex(0,1,1);
        M.facets.set_vertex(0,2,2);
        // triangle 1
        M.facets.set_vertex(1,0,3);
        M.facets.set_vertex(1,1,2);
        M.facets.set_vertex(1,2,1);
        // triangle 2
        M.facets.set_vertex(2,0,1);
        M.facets.set_vertex(2,1,0);
        M.facets.set_vertex(2,2,4);
        // triangle 3
        M.facets.set_vertex(3,0,5);
        M.facets.set_vertex(3,1,4);
        M.facets.set_vertex(3,2,0);
        // triangle 4
        M.facets.set_vertex(4,0,3);
        M.facets.set_vertex(4,1,1);
        M.facets.set_vertex(4,2,6);
        // triangle 5
        M.facets.set_vertex(5,0,4);
        M.facets.set_vertex(5,1,6);
        M.facets.set_vertex(5,2,1);
        // triangle 6
        M.facets.set_vertex(6,0,2);
        M.facets.set_vertex(6,1,3);
        M.facets.set_vertex(6,2,7);
        // triangle 7
        M.facets.set_vertex(7,0,6);
        M.facets.set_vertex(7,1,7);
        M.facets.set_vertex(7,2,3);
        // triangle 8
        M.facets.set_vertex(8,0,0);
        M.facets.set_vertex(8,1,2);
        M.facets.set_vertex(8,2,5);
        // triangle 9
        M.facets.set_vertex(9,0,7);
        M.facets.set_vertex(9,1,5);
        M.facets.set_vertex(9,2,2);
        // triangle 10
        M.facets.set_vertex(10,0,5);
        M.facets.set_vertex(10,1,7);
        M.facets.set_vertex(10,2,4);
        // triangle 11
        M.facets.set_vertex(11,0,6);
        M.facets.set_vertex(11,1,4);
        M.facets.set_vertex(11,2,7);
        // compute facets adjacency
        M.facets.connect();
    }

    // ensure facet normals are outward
	if(facet_normals_are_inward(M)) {
		flip_facet_normals(M);
		fmt::println(GEO::Logger::warn("normals dir."),"Facet normals of the input mesh were inward");
		fmt::println(GEO::Logger::warn("normals dir."),"You should flip them and update the file for consistency.");
		GEO::Logger::warn("normals dir.").flush();
	}

    // compute facet normals
	std::vector<vec3> normals(M.facets.nb());
	FOR(f,M.facets.nb()) {
		normals[f] = normalize(Geom::mesh_facet_normal(M,f));
	}

    std::vector<std::vector<AdjacentFacetOfVertex>> per_vertex_adj_facets;
    compute_adjacent_facets_of_vertices(M,per_vertex_adj_facets);

    // compute naive labeling
    naive_labeling(M,normals,"label");

    // export to glTF
    write_glTF__labeled_triangle_mesh(filenames[0],M,"label",per_vertex_adj_facets);
    return 0;
}