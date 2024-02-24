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

    std::vector<std::string> filenames;
    if(!CmdLine::parse(argc,argv,filenames,"input_mesh"))	{
		return 1;
	}

    GEO::Mesh M;
    mesh_load(filenames[0],M);

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