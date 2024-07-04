// Flip the normals of a triangle mesh
// that is, outward -> inward or inward -> outward (assuming all have the same direction)
// If no output mesh filename is provided, determine the normals direction

#include <geogram/mesh/mesh.h>  // for GEO::Mesh
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/numeric.h>              // for max_float64()

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>

#include <vector>
#include <string>
#include <set>

#include "geometry_halfedges.h"     // for CustomMeshHalfedges
#include "geometry.h"               // for facet_normals_are_inward() & flip_facet_normals()

using namespace GEO;

int main(int argc, char** argv) {

    std::vector<std::string> filenames;
	GEO::initialize();
	if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh <output_mesh>"
		))
	{
        fmt::println(Logger::err("I/O"),"Usage should be\n{} input_mesh <output_mesh>",argv[0]); Logger::err("I/O").flush();
		return 1;
	}

    Mesh M;
    geo_assert(FileSystem::is_file(filenames[0]));

    if(!GEO::mesh_load(filenames[0],M)) {
        fmt::println("Unable to open {}",filenames[0]);
        return 1;
    }
    geo_assert(M.cells.nb() == 0); // must be a surface mesh
    geo_assert(M.facets.nb() != 0); // must have facets
    geo_assert(M.facets.are_simplices()); // must be a triangle mesh

    if(filenames.size()==1) {
        // no output mesh -> determine normals direction
        if(facet_normals_are_inward(M)) {
            fmt::println(Logger::out("normals dir."),"The facet normals are inward"); Logger::out("normals dir.").flush();
        }
        else {
            fmt::println(Logger::out("normals dir."),"The facet normals are outward"); Logger::out("normals dir.").flush();
        }
        return 0;
    }

    flip_facet_normals(M);

    mesh_save(M,filenames[1]);

    return 0;
}