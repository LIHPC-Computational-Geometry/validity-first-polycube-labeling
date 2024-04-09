// Computes distortion metrics of the polycube deformation
// Inputs: triangle mesh + polycube mesh (generated with https://github.com/fprotais/fastbndpolycube)
// Output: stretch, area distortion, angle distortion and isometric distortion - as JSON written on stdout
// Based on
// - https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/app/measurement.cpp
// - https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp
// There:
// - F, V1 the facets and vertices of the input mesh
// - F2, V2 the facets and vertices of the polycube mesh
// - A1 the per facet area of the input mesh, and A_m the sum
// - A2 the per facet area of the polycube mesh, and A_d the sum
// - N the per facet normal of the input mesh
// - N_def the per facet normal of the polycube mesh
// - vecA1 the per facet local transformation of the input mesh
// - vecA2 the per facet local transformation of the polycube mesh

#include <geogram/mesh/mesh.h>              // for Mesh
#include <geogram/mesh/mesh_io.h>           // for mesh_load()
#include <geogram/mesh/mesh_geometry.h>     // for mesh_facet_area(), mesh_facet_normal()
#include <geogram/basic/file_system.h>      // for FileSystem::is_file()
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>             // for length()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <Eigen/SVD>

#include <vector>

#include "geometry.h" // for facet_normals_are_inward(), flip_facet_normals(), center_mesh(), compute_jacobians()

using namespace GEO;

// Geogram matrix to Eigen matrix
// Only for square matrices
template <index_t DIM, class T>
void fill_Eigen_from_Geogram_matrix(const GEO::Matrix<DIM,T>& in, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& out) {
    out.resize(DIM,DIM);
    FOR(r,DIM) { // for each row
        FOR(c,DIM) { // for each column
            out(r,c) = in(r,c);
        }
    }
}

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
		"input_mesh polycube_mesh"
		))
	{
		return 1;
	}

    // Check if the input mesh and the polycube mesh exist

    if(!FileSystem::is_file(filenames[0])) {
        fmt::println(Logger::err("I/O"),"{} does not exist",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }

    if(!FileSystem::is_file(filenames[1])) {
        fmt::println(Logger::err("I/O"),"{} does not exist",filenames[1]); Logger::err("I/O").flush();
        return 1;
    }

    // load them

    mesh_io_initialize();
    Mesh input_mesh, polycube_mesh;
    if(!mesh_load(filenames[0],input_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to load the mesh from {}",filenames[0]); Logger::err("I/O").flush();
        return 1;
    }
    if(!mesh_load(filenames[1],polycube_mesh)) {
        fmt::println(Logger::err("I/O"),"Unable to load the mesh from {}",filenames[1]); Logger::err("I/O").flush();
        return 1;
    }

    // ensure the normals are outward

    if(facet_normals_are_inward(input_mesh)) {
        flip_facet_normals(input_mesh);
        fmt::println(Logger::warn("I/O"),"Facet normals of {} were inward and have been flipped", filenames[0]); Logger::warn("I/O").flush();
    }
    if(facet_normals_are_inward(polycube_mesh)) {
        flip_facet_normals(polycube_mesh);
        fmt::println(Logger::warn("I/O"),"Facet normals of {} were inward and have been flipped", filenames[1]); Logger::warn("I/O").flush();
    }

    // check consistency between the meshes

    geo_assert(input_mesh.vertices.nb() == polycube_mesh.vertices.nb());
    geo_assert(input_mesh.facets.nb() == polycube_mesh.facets.nb());
    geo_assert(input_mesh.facets.are_simplices());
    geo_assert(polycube_mesh.facets.are_simplices());
    geo_assert(input_mesh.cells.nb() == 0);
    geo_assert(polycube_mesh.cells.nb() == 0);
    FOR(f,input_mesh.facets.nb()) { // for each facet index
        FOR(lv,3) { // for each local vertex
            geo_assert(input_mesh.facets.vertex(f,lv) == polycube_mesh.facets.vertex(f,lv));
        }
    }

    // center and normalize the meshes

    center_mesh(input_mesh,true);
    center_mesh(polycube_mesh,true);

    // compute facets normal

	std::vector<vec3> input_mesh_per_facet_normals(input_mesh.facets.nb());
	FOR(f,input_mesh.facets.nb()) {
		input_mesh_per_facet_normals[f] = normalize(Geom::mesh_facet_normal(input_mesh,f));
	}

    std::vector<vec3> polycube_mesh_per_facet_normals(polycube_mesh.facets.nb());
    FOR(f,polycube_mesh.facets.nb()) {
		polycube_mesh_per_facet_normals[f] = normalize(Geom::mesh_facet_normal(polycube_mesh,f));
	}

    // compute facets area

    std::vector<double> input_mesh_per_facet_area(input_mesh.facets.nb());
    double input_mesh_total_area = 0.0;
    FOR(f,input_mesh.facets.nb()) {
        input_mesh_per_facet_area[f] = mesh_facet_area(input_mesh,f);
        input_mesh_total_area += input_mesh_per_facet_area[f];
    }

    std::vector<double> polycube_mesh_per_facet_area(polycube_mesh.facets.nb());
    double polycube_mesh_total_area = 0.0;
    FOR(f,polycube_mesh.facets.nb()) {
        polycube_mesh_per_facet_area[f] = mesh_facet_area(polycube_mesh,f);
        polycube_mesh_total_area += polycube_mesh_per_facet_area[f];
    }

    // compute jacobians

    std::vector<mat2> jacobians(input_mesh.facets.nb());
    compute_jacobians(input_mesh,polycube_mesh,input_mesh_per_facet_normals,polycube_mesh_per_facet_normals,jacobians);

    // compute per facet singular values
    std::vector<std::pair<double, double>> per_facet_singular_values(input_mesh.facets.nb(),{0.0,0.0});
    Eigen::MatrixXd current_matrix;
    FOR(f,input_mesh.facets.nb()){
        fill_Eigen_from_Geogram_matrix(jacobians[f],current_matrix);
        // Note: We could use an Eigen::Matrix2d
        // but Evocube configure the solver with Eigen::ComputeThinU | Eigen::ComputeThinV (which doesn't seem necessary)
        // and it's incompatible: "thin U and V are only available when your matrix has a dynamic number of columns"
        // So I use an Eigen::MatrixXd...
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(current_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
        per_facet_singular_values[f].first = svd.singularValues()(0);
        per_facet_singular_values[f].second = svd.singularValues()(1);
    }

    // TODO compute Stretch
    // TODO compute Area distortion
    // TODO compute Angle distortion
    // TODO compute Isometric distortion

    return 0;
}