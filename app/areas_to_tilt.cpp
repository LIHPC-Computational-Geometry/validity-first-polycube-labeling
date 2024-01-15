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

#include "SimpleMeshApplicationExt.h"   // for colormaps indices
#include "about_window.h"               // for draw_about_window()

using namespace GEO;

class AreasToTiltApp : public SimpleMeshApplicationExt {
public:

    AreasToTiltApp() : SimpleMeshApplicationExt("areas_to_tilt") {
        // change default values
        show_surface_ = true;
        show_volume_ = false;
        show_colored_cells_ = false;
        // init own variables
        sensitivity_ = 0.3f; // 0.02f
        group_by_area_ = false;
        nb_areas = 0;
    }

protected:

    void draw_object_properties() override {
        ImGui::SliderFloat("Sensitivity",&sensitivity_,0.0f,1.0f,"%0.4f");
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

        // compute facet normals
        normals_.resize(mesh_.facets.nb());
        FOR(f,mesh_.facets.nb()) { // for each facet
            normals_[f] = normalize(Geom::mesh_facet_normal(mesh_,f));
        }

        compute_areas_to_tilt();

        update_drawing_settings();

        return true;
    }

    void compute_areas_to_tilt() {

        Attribute<bool> on_area_to_tilt(mesh_.facets.attributes(),"on_area_to_tilt");
        unsigned int nb_facet_on_areas_to_tilt = 0;
        FOR(f,mesh_.facets.nb()) { // for each facet
            const vec3& normal = normals_[f];
            // based on labeling.cpp nearest_label()
            std::array<std::pair<double,index_t>,6> per_label_weights = {
                std::make_pair(normal.x < 0.0 ? 0.0 : normal.x,     0), // +X
                std::make_pair(normal.x < 0.0 ? -normal.x : 0.0,    1), // -X
                std::make_pair(normal.y < 0.0 ? 0.0 : normal.y,     2), // +Y
                std::make_pair(normal.y < 0.0 ? -normal.y : 0.0,    3), // -Y
                std::make_pair(normal.z < 0.0 ? 0.0 : normal.z,     4), // +Z
                std::make_pair(normal.z < 0.0 ? -normal.z : 0.0,    5)  // -Z
            };
            std::sort(per_label_weights.begin(),per_label_weights.end());
            if(std::abs(per_label_weights[5].first-per_label_weights[4].first) < (double) sensitivity_) {
                // the 2 labels with the most weight are too close
                nb_facet_on_areas_to_tilt++;
                on_area_to_tilt[f] = true;
            }
            else {
                on_area_to_tilt[f] = false;
            }
        }

        fmt::println(Logger::out("areas to tilt"),"{} facet(s) over {} are on an area to tilt",nb_facet_on_areas_to_tilt,mesh_.facets.nb()); Logger::out("areas to tilt").flush();

        DisjointSet<index_t> ds(mesh_.facets.nb());
        index_t neighbor_facet = index_t(-1);
        FOR(current_facet,mesh_.facets.nb()) {
            FOR(le,3) { // for each of the 3 local edges
                neighbor_facet = mesh_.facets.adjacent(current_facet,le);
                if(on_area_to_tilt[current_facet] == on_area_to_tilt[neighbor_facet]) {
                    ds.mergeSets(current_facet,neighbor_facet);
                }
            }
        }

        Attribute<index_t> area_index(mesh_.facets.attributes(),"area_index");
        nb_areas = ds.getSetsMap(area_index.data());

        std::vector<index_t> color_randomizer(nb_areas);
        FOR(a,nb_areas)
            color_randomizer[a] = a;
        std::random_shuffle(color_randomizer.begin(),color_randomizer.end());
        FOR(f,mesh_.facets.nb())
            area_index[f] = color_randomizer[area_index[f]];

        fmt::println(Logger::out("areas to tilt"),"facets grouped into {} area(s)",nb_areas); Logger::out("areas to tilt").flush();
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

    std::vector<vec3> normals_; // per facet normals
    float sensitivity_;
    bool group_by_area_;
    std::size_t nb_areas;

};

int main(int argc, char** argv) {
    AreasToTiltApp app;
	app.start(argc,argv);
    return 0;
}