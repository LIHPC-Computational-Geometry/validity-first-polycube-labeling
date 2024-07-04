#include "geometry_distortion.h"

void compute_per_facet_local_transfo(const Mesh& M, const std::vector<vec3>& facets_normal, std::vector<mat2>& out) {
    // Based on Evocube > src/distortion.cpp > localTransfo()
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp#L13
    // There:
    // - F(f_id, 1) is the vertex index at the local vertex 1 -> v1
    // - F(f_id, 0) is the vertex index at the local vertex 0 -> v0
    // - local_axis1 is the normalized 3D vector going from v0 to v1
    // - local_axis2 is the cross product between the facet normal and local_axis1
    // - local_coords1[0] is the norm of local_axis1 (before normalization), and local_coords1[1] is 0
    // - local_coords2[0] is dot product of the 3D vector going from v0 to v2 with local_axis1
    // - local_coords2[1] is dot product of the 3D vector going from v0 to v2 with local_axis2
    // - A is the resulting 2x2 matrix for the current facet -> out[f]
    // - column 0 of A is local_coords1
    // - column 1 of A is local_coords2
    out.resize(M.facets.nb());
    index_t v0 = index_t(-1);
    index_t v1 = index_t(-1);
    index_t v2 = index_t(-1);
    vec3 local_axis1;
    vec3 local_axis2;
    vec2 local_coords1;
    vec2 local_coords2;
    mat2 A;
    FOR(f,M.facets.nb()) {
        v0 = M.facets.vertex(f,0);
        v1 = M.facets.vertex(f,1);
        v2 = M.facets.vertex(f,2);
        local_axis1 = M.vertices.point(v1) - M.vertices.point(v0);
        local_coords1 = vec2(local_axis1.length(),0.0);
        local_axis1 = normalize(local_axis1);
        local_axis2 = cross(facets_normal[f],local_axis1);
        local_coords2[0] = dot(M.vertices.point(v2) - M.vertices.point(v0), local_axis1);
        local_coords2[1] = dot(M.vertices.point(v2) - M.vertices.point(v0), local_axis2);
        // fill column 0
        out[f](0,0) = local_coords1[0];
        out[f](1,0) = local_coords1[1];
        // fill column 1
        out[f](0,1) = local_coords2[0];
        out[f](1,1) = local_coords2[1];
    }
}

void compute_jacobians(const Mesh& M1, const Mesh& M2, const std::vector<vec3>& M1_normals, const std::vector<vec3>& M2_normals, std::vector<mat2>& out) {
    // Based on Evocube > src/distortion.cpp > computeJacobians()
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp#L38
    
    // check consistency between the meshes

    geo_assert(M1.vertices.nb() == M2.vertices.nb());
    geo_assert(M1.facets.nb() == M2.facets.nb());
    geo_assert(M1.facets.are_simplices());
    geo_assert(M2.facets.are_simplices());
    geo_assert(M1.cells.nb() == 0);
    geo_assert(M2.cells.nb() == 0);
    FOR(f,M1.facets.nb()) { // for each facet index
        FOR(lv,3) { // for each local vertex
            geo_assert(M1.facets.vertex(f,lv) == M2.facets.vertex(f,lv));
        }
    }

    // compute local transformation

    std::vector<mat2> M1_per_facet_local_transfo;
    compute_per_facet_local_transfo(M1,M1_normals,M1_per_facet_local_transfo);

    std::vector<mat2> M2_per_facet_local_transfo;
    compute_per_facet_local_transfo(M2,M2_normals,M2_per_facet_local_transfo);

    out.resize(M1.facets.nb());
    FOR(f,M1.facets.nb()) {
        out[f] = M2_per_facet_local_transfo[f] * M1_per_facet_local_transfo[f].inverse();
    }
}

double compute_stretch(const std::vector<double>& input_mesh_per_facet_area, double input_mesh_total_area, double polycube_mesh_total_area, const std::vector<std::pair<double, double>>& per_facet_singular_values) {
    // based on Evocube > src/distortion.cpp > computeStretch()
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp#L97
    // "Spherical Parametrization and Remeshing", Praun & Hoppe, 2003
    // Section 3.3 : Review of planar-domain stretch metric ?
    // Itself based on "Texture Mapping Progressive Meshes", Sander & Snyder & Gortler & Hoppe, 2001
    // Section 3 : Texture stretch metric ?
    geo_assert(input_mesh_per_facet_area.size() == per_facet_singular_values.size());
    double sum = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    FOR(f,per_facet_singular_values.size()) {
        // get singular values of this facet
        s1 = per_facet_singular_values[f].first;
        s2 = per_facet_singular_values[f].second;
        // update sum
        sum += input_mesh_per_facet_area[f] * (s1 * s1 + s2 * s2) / 2.0;
    }
    sum /= input_mesh_total_area;
    return (input_mesh_total_area / polycube_mesh_total_area) * (1.0 / std::pow(sum, 2));
}

double compute_area_distortion(
    const std::vector<double>& input_mesh_per_facet_area, // eq. to A
    double input_mesh_total_area, // eq. to A.sum()
    const std::vector<std::pair<double, double>>& per_facet_singular_values // eq. to per_tri_singular_values
) {
    // based on Evocube > src/distortion.cpp > computeAreaDisto()
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp#L113
    // "PolyCube-Maps", Tarini & Hormann & Cignoni & Montani
    // Section 7, Table 1
    // integration & normalisation of s1*s2+1/(s1*s2)
    // Itself referring to
    // - "An adaptable surface parametrization method", Degener & Meseth & Klein, 2003
    // - "Surface parametrization: a tutorial and survey", Floater & Hormann, 2004
    geo_assert(input_mesh_per_facet_area.size() == per_facet_singular_values.size());
    double sum = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    FOR(f,per_facet_singular_values.size()) {
        // get singular values of this facet
        s1 = per_facet_singular_values[f].first;
        s2 = per_facet_singular_values[f].second;
        // update sum
        sum += input_mesh_per_facet_area[f] * 0.5 * (s1 * s2 + 1.0 / (s1 * s2));
    }
    return sum / input_mesh_total_area;
}

double compute_angle_distortion(
    const std::vector<double>& input_mesh_per_facet_area,
    double input_mesh_total_area,
    const std::vector<std::pair<double, double>>& per_facet_singular_values
) {
    // based on Evocube > src/distortion.cpp > computeAngleDisto()
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp#L124
    // "PolyCube-Maps", Tarini & Hormann & Cignoni & Montani
    // Section 7, Table 1
    // integration & normalisation of (s1/s2)+(s2/s1)
    // Itself referring to
    // - "An adaptable surface parametrization method", Degener & Meseth & Klein, 2003
    // - "Surface parametrization: a tutorial and survey", Floater & Hormann, 2004
    geo_assert(input_mesh_per_facet_area.size() == per_facet_singular_values.size());
    double sum = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    FOR(f,per_facet_singular_values.size()) {
        // get singular values of this facet
        s1 = per_facet_singular_values[f].first;
        s2 = per_facet_singular_values[f].second;
        // update sum
        sum += input_mesh_per_facet_area[f] * 0.5 * (s1 / s2 + s2 / s1);
    }
    return sum / input_mesh_total_area;
}

double compute_isometric_distortion(
    const std::vector<double>& input_mesh_per_facet_area,
    double input_mesh_total_area,
    const std::vector<std::pair<double, double>>& per_facet_singular_values
) {
    // based on Evocube > src/distortion.cpp > computeIsometricDisto()
    // https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/distortion.cpp#L137
    // "Computing Surface PolyCube-Maps by Constrained Voxelization", Yang & Fu & Liu
    geo_assert(input_mesh_per_facet_area.size() == per_facet_singular_values.size());
    double sum = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    FOR(f,per_facet_singular_values.size()) {
        // get singular values of this facet
        s1 = per_facet_singular_values[f].first;
        s2 = per_facet_singular_values[f].second;
        // update sum
        sum += std::max(s1, 1.0 / s2) * input_mesh_per_facet_area[f];
    }
    return sum / input_mesh_total_area;
}