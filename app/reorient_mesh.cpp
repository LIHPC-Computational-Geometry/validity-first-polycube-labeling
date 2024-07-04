/* Reorient a mesh according to the principal axes of its point cloud.
 * Not as good as "Alignment of 3D models" by Chaouch & Verroust-Blondet https://inria.hal.science/hal-00804653
 * But in our case, misalignment is good because there might be less surfaces where 2 signed axes are equaly close to the normals
 */

#include <geogram/mesh/mesh.h>              // for Mesh
#include <geogram/mesh/mesh_io.h>           // for mesh_load(), mesh_save()
#include <geogram/basic/file_system.h>      // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <array>

#include "stats.h" // for BasicStats
#include "geometry.h"

using namespace GEO;

int main(int argc, char** argv) {
    
    GEO::initialize();
    std::vector<std::string> filenames;
    if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh output_mesh"
		))
	{
		return 1;
	}

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

    if(input_mesh.vertices.nb() == 0) {
        fmt::println(Logger::err("reorient"),"The input mesh has no vertices"); Logger::err("reorient");
        return 1;
    }

    rotate_mesh_according_to_principal_axes(input_mesh);

    mesh_save(input_mesh,filenames[1]);

    return 0;
}