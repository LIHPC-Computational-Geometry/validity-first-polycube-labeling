/* Reorient a mesh according to the principal axes of its point cloud.
 * Not as good as "Alignment of 3D models" by Chaouch & Verroust-Blondet https://inria.hal.science/hal-00804653
 * But in our case, misalignment is good because there might be less surfaces where 2 signed axes are very close to the normals
 */

#include <geogram/mesh/mesh.h>              // for Mesh
#include <geogram/mesh/mesh_io.h>           // for mesh_load(), mesh_save()
#include <geogram/basic/file_system.h>      // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/points/principal_axes.h>  // for PrincipalAxes3d
#include <geogram/basic/vecg.h>             // for vec3, length()
#include <geogram/basic/matrix.h>           // for mat3

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <array>

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

    ////////////////////////////////
    // Compute principal axes
    ////////////////////////////////
    
    if(input_mesh.vertices.nb() == 0) {
        fmt::println(Logger::err("reorient"),"The input mesh has no vertices"); Logger::err("reorient");
        return 1;
    }
    PrincipalAxes3d principal_axes;
    principal_axes.begin();
    FOR(v,input_mesh.vertices.nb()) {
        principal_axes.add_point(input_mesh.vertices.point(v));
    }
    principal_axes.end();

    ////////////////////////////////
    // Compute vector toward (1,1,1) in the principal axes
    ////////////////////////////////

    vec3 rotated_vector = normalize(principal_axes.axis(0)) +
                          normalize(principal_axes.axis(1)) + 
                          normalize(principal_axes.axis(2));
    
    ////////////////////////////////
    // Compute rotation matrix that align (1,1,1) (regarding init axes) to rotated_vector
    ////////////////////////////////

    // thanks Jur van den Berg https://math.stackexchange.com/a/476311
    geo_assert(length(rotated_vector - vec3(-1.0,-1.0,-1.0)) > 10e-3); // formula not applicable if the 2 vectors are opposite
    vec3 v = cross(vec3(1.0,1.0,1.0),rotated_vector);
    mat3 skew_symmetric_cross_product;
    skew_symmetric_cross_product(0,0) = 0;
    skew_symmetric_cross_product(0,1) = -v[2];
    skew_symmetric_cross_product(0,2) = v[1];
    skew_symmetric_cross_product(1,0) = v[2];
    skew_symmetric_cross_product(1,1) = 0;
    skew_symmetric_cross_product(1,2) = -v[0];
    skew_symmetric_cross_product(2,0) = -v[1];
    skew_symmetric_cross_product(2,1) = v[0];
    skew_symmetric_cross_product(2,2) = 0;
    double c = dot(vec3(1.0,1.0,1.0),rotated_vector);
    mat3 R; // initialized with identity
    R += skew_symmetric_cross_product + skew_symmetric_cross_product * skew_symmetric_cross_product * (1/(1+c));

    ////////////////////////////////
    // Rotate point cloud
    ////////////////////////////////

    FOR(v,input_mesh.vertices.nb()) {
        input_mesh.vertices.point(v) = mult(R,input_mesh.vertices.point(v));
    }

    ////////////////////////////////
    // Save output mesh
    ////////////////////////////////

    mesh_save(input_mesh,filenames[1]);

    return 0;
}