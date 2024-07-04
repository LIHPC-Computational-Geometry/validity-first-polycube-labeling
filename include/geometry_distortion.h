#pragma once

#include <geogram/mesh/mesh.h> // for GEO::Mesh

#include <vector>

using namespace GEO;

void compute_per_facet_local_transfo(const Mesh& M, const std::vector<vec3>& facets_normal, std::vector<mat2>& out);

void compute_jacobians(
    const Mesh& M1,
    const Mesh& M2,
    const std::vector<vec3>& M1_normals,
    const std::vector<vec3>& M2_normals,
    std::vector<mat2>& out
);

double compute_stretch(
    const std::vector<double>& input_mesh_per_facet_area,
    double input_mesh_total_area,
    double polycube_mesh_total_area,
    const std::vector<std::pair<double, double>>& per_facet_singular_values
);

double compute_area_distortion(
    const std::vector<double>& input_mesh_per_facet_area,
    double input_mesh_total_area,
    const std::vector<std::pair<double, double>>& per_facet_singular_values
);

double compute_angle_distortion(
    const std::vector<double>& input_mesh_per_facet_area,
    double input_mesh_total_area,
    const std::vector<std::pair<double, double>>& per_facet_singular_values
);

double compute_isometric_distortion(
    const std::vector<double>& input_mesh_per_facet_area,
    double input_mesh_total_area,
    const std::vector<std::pair<double, double>>& per_facet_singular_values
);