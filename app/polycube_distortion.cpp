// Computes distortion metrics of the polycube deformation
// Inputs: triangle mesh + polycube mesh (generated with https://github.com/fprotais/fastbndpolycube)
// Output: stretch, area distortion, angle distortion and isometric distortion
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
//
// Usage:
//
//   ./bin/polycube_distortion triangle_mesh.obj polycube_mesh.obj
//     -> write to stdout
//
//   ./bin/polycube_distortion triangle_mesh.obj polycube_mesh.obj output.json
//     -> write to output.json

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

#include <Eigen/SVD> // for Eigen::JacobiSVD<>

#include <vector>
#include <fstream>
#include <iomanip> // for std::setw()

#include <nlohmann/json.hpp>

#ifndef NDEBUG
    #include <dbg.h>
#endif

#include "geometry_distortion.h" // for compute_jacobians(), compute_stretch(), compute_area_distortion(), compute_angle_distortion(), compute_isometric_distortion()
#include "containers_macros.h" // for VECTOR_SUM()
#include "containers_Geogram.h" // to print a GEO::mat2 with {fmt}

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
    
    GEO::initialize();

    std::vector<std::string> filenames;
    if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh polycube_mesh <output_JSON>"
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

    #ifndef NDEBUG
        fmt::println("first 5 triangles");
        FOR(f,5) {
            index_t v0 = input_mesh.facets.vertex(f,0);
            index_t v1 = input_mesh.facets.vertex(f,1);
            index_t v2 = input_mesh.facets.vertex(f,2);
            fmt::println("v0 = {} = ({}, {}, {})",v0,input_mesh.vertices.point(v0).x,input_mesh.vertices.point(v0).y,input_mesh.vertices.point(v0).z);
            fmt::println("v1 = {} = ({}, {}, {})",v1,input_mesh.vertices.point(v1).x,input_mesh.vertices.point(v1).y,input_mesh.vertices.point(v1).z);
            fmt::println("v2 = {} = ({}, {}, {})",v2,input_mesh.vertices.point(v2).x,input_mesh.vertices.point(v2).y,input_mesh.vertices.point(v2).z);
            fmt::println("");
        }
    #endif

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
        input_mesh_per_facet_area[f] = Geom::mesh_facet_area(input_mesh,f);
        input_mesh_total_area += input_mesh_per_facet_area[f];
    }
    #ifndef NDEBUG
        dbg(input_mesh_total_area);
        geo_assert(input_mesh_total_area == VECTOR_SUM(input_mesh_per_facet_area));
    #endif

    std::vector<double> polycube_mesh_per_facet_area(polycube_mesh.facets.nb());
    double polycube_mesh_total_area = 0.0;
    FOR(f,polycube_mesh.facets.nb()) {
        polycube_mesh_per_facet_area[f] = Geom::mesh_facet_area(polycube_mesh,f);
        polycube_mesh_total_area += polycube_mesh_per_facet_area[f];
    }
    #ifndef NDEBUG
        dbg(polycube_mesh_total_area);
        geo_assert(polycube_mesh_total_area == VECTOR_SUM(polycube_mesh_per_facet_area));
    #endif

    // scale the polycube mesh so that both meshes have the same total surface area

    FOR(v,polycube_mesh.vertices.nb()) {
        polycube_mesh.vertices.point(v) *= std::sqrt(input_mesh_total_area / polycube_mesh_total_area);
    }
    #ifndef NDEBUG
        dbg("polycube mesh rescaled");
    #endif

    // recompute facets area of the polycube mesh

    polycube_mesh_total_area = 0.0;
    FOR(f,polycube_mesh.facets.nb()) {
        polycube_mesh_per_facet_area[f] = Geom::mesh_facet_area(polycube_mesh,f);
        polycube_mesh_total_area += polycube_mesh_per_facet_area[f];
    }
    #ifndef NDEBUG
        dbg(polycube_mesh_total_area);
        geo_assert(polycube_mesh_total_area == VECTOR_SUM(polycube_mesh_per_facet_area));
    #endif

    // compute jacobians

    std::vector<mat2> jacobians(input_mesh.facets.nb());
    compute_jacobians(input_mesh,polycube_mesh,input_mesh_per_facet_normals,polycube_mesh_per_facet_normals,jacobians);
    #ifndef NDEBUG
        fmt::println("first 5 Jacobians");
        for(std::vector<mat2>::size_type f_id = 0; f_id < 5; f_id++) {
            fmt::println("{}",jacobians[f_id]);
        }
        fmt::println("");
    #endif

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
    #ifndef NDEBUG
        fmt::println("first 5 singular values");
        for(std::vector<std::pair<double, double>>::size_type f_id = 0; f_id < 5; f_id++) {
            fmt::println("{}, {}", per_facet_singular_values[f_id].first, per_facet_singular_values[f_id].second);
        }
        fmt::println("");
    #endif

    double stretch = compute_stretch(input_mesh_per_facet_area,input_mesh_total_area,polycube_mesh_total_area,per_facet_singular_values);
    double area_distortion = compute_area_distortion(input_mesh_per_facet_area,input_mesh_total_area,per_facet_singular_values);
    double angle_distortion = compute_angle_distortion(input_mesh_per_facet_area,input_mesh_total_area,per_facet_singular_values);
    double isometric_distortion = compute_isometric_distortion(input_mesh_per_facet_area,input_mesh_total_area,per_facet_singular_values);

    // write to stdout or file, according to the number of arguments

    if (filenames.size() <= 2) {
        // no output JSON filename provided
        fmt::println(Logger::out("distortions"), "             stretch = {}", stretch);
        fmt::println(Logger::out("distortions"), "     area_distortion = {}", area_distortion);
        fmt::println(Logger::out("distortions"), "    angle_distortion = {}", angle_distortion);
        fmt::println(Logger::out("distortions"), "isometric_distortion = {}", isometric_distortion);
        Logger::out("distortions").flush();
    }
    else {
        // write values in a JSON file
        nlohmann::json output_JSON;
        output_JSON["stretch"] = stretch;
        output_JSON["area_distortion"] = area_distortion;
        output_JSON["angle_distortion"] = angle_distortion;
        output_JSON["isometric_distortion"] = isometric_distortion;

        std::fstream ofs(filenames[2],std::ios_base::out);
        if(ofs.good()) {
            fmt::println(Logger::out("I/O"),"Saving file {}...",filenames[2]); Logger::out("I/O").flush();
            ofs << std::setw(4) << output_JSON << std::endl;
        }
        else {
            fmt::println(Logger::err("I/O"),"Cannot write into {}",filenames[2]); Logger::err("I/O").flush();
        }
    }

    return 0;
}