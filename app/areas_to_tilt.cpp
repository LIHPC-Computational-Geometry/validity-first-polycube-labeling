/* An app to identify areas on a 3D triangle mesh where, among {+X,-X,+Y,-Y,+Z,-Z}, several directions are very close to the normal.
 * A suitable solution would be to slightly tilt each area, taking into account dependencies (global consistency of the rotations) to avoid overlapping polycube parts after the tilt.
 * A basic solution is tweaked_naive_labeling() in src/labeling.cpp (apply the same rotation for all of those areas to tilt)
 */

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <DisjointSet.hpp>

#include <vector>
#include <array>
#include <map>
#include <utility> // for std::pair<>
#include <algorithm> // for std::random_shuffle(), std::sort

#include "gui_base.h" // for SimpleMeshApplicationExt, colormaps indices
#include "geometry.h" // for is_a_facet_to_tilt()
#include "labeling_generators.h" // for DEFAULT_SENSITIVITY

using namespace GEO;

class AreasToTiltApp : public SimpleMeshApplicationExt {
public:

    AreasToTiltApp() : SimpleMeshApplicationExt("areas_to_tilt"), mesh_ext_(mesh_) {
        // change default values
        show_surface_ = true;
        show_volume_ = false;
        show_colored_cells_ = false;
        // init own variables
        sensitivity_ = DEFAULT_SENSITIVITY;
        group_by_area_ = false;
        nb_areas = 0;
    }

protected:

    void draw_object_properties() override {
        ImGui::InputDouble("Sensitivity", &sensitivity_, 0.0, 0.0, "%.15f");
        if(ImGui::Checkbox("Group by area",&group_by_area_)) {
            update_drawing_settings();
        }
        ImGui::Text("nb_areas=%lu",nb_areas);
        if(ImGui::Button("Apply")) {
            compute_areas_to_tilt();
        }
    }

    bool load(const std::string& filename) override {

        if(!SimpleMeshApplicationExt::load(filename)) {
            lighting_ = true;
            show_attributes_ = false;
            return false;
        }

        // expect a surface triangle mesh
        geo_assert(mesh_.facets.nb() != 0);
        geo_assert(mesh_.facets.are_simplices());
        geo_assert(mesh_.cells.nb() == 0);

        mesh_.vertices.set_double_precision();

        mesh_ext_.facet_normals.recompute();
        compute_areas_to_tilt();

        return true;
    }

    void compute_areas_to_tilt() {

        std::set<index_t> facets_to_tilt;
        size_t nb_facet_on_areas_to_tilt = get_facets_to_tilt(mesh_ext_,facets_to_tilt,sensitivity_);
        Attribute<bool> on_area_to_tilt(mesh_.facets.attributes(),"on_area_to_tilt");
        FOR(f,mesh_.facets.nb()) { // for each facet
            on_area_to_tilt[f] = facets_to_tilt.contains(f);
        }

        fmt::println(Logger::out("areas to tilt"),"{} facet(s) over {} are on an area to tilt",nb_facet_on_areas_to_tilt,mesh_.facets.nb()); Logger::out("areas to tilt").flush();

        Attribute<index_t> area_index(mesh_.facets.attributes(),"area_index");
        nb_areas = group_facets<Attribute>(mesh_,facets_to_tilt,area_index);

        std::vector<index_t> color_randomizer(nb_areas);
        FOR(a,nb_areas) {
            color_randomizer[a] = a;
        }
        std::random_shuffle(color_randomizer.begin(),color_randomizer.end());
        FOR(f,mesh_.facets.nb()) {
            area_index[f] = color_randomizer[area_index[f]];
        }

        fmt::println(Logger::out("areas to tilt"),"facets grouped into {} area(s)",nb_areas); Logger::out("areas to tilt").flush();
    }

    void GL_initialize() override {
        SimpleMeshApplicationExt::GL_initialize();
        update_drawing_settings();
    }

    void update_drawing_settings() {
        if(group_by_area_) {
            lighting_ = false;
            show_attributes_ = true;
            current_colormap_index_ =COLORMAP_CEI_60757;
            attribute_min_ = 0.0f;
            attribute_max_ = (float) nb_areas+1;
            attribute_ = "facets.area_index";
            attribute_name_ = "area_index";
            attribute_subelements_ = MESH_FACETS;
        }
        else {
            lighting_ = false;
            show_attributes_ = true;
            current_colormap_index_ = COLORMAP_BLUE_RED;
            attribute_min_ = 0.0f;
            attribute_max_ = 1.0f;
            attribute_ = "facets.on_area_to_tilt";
            attribute_name_ = "on_area_to_tilt";
            attribute_subelements_ = MESH_FACETS;
        }
    }

    MeshExt mesh_ext_;
    double sensitivity_;
    bool group_by_area_;
    std::size_t nb_areas;

};

int main(int argc, char** argv) {
    AreasToTiltApp app;
	app.start(argc,argv);
    return 0;
}