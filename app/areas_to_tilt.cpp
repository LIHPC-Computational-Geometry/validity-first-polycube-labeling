#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/command_line.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <array>
#include <utility> // for std::pair<>
#include <algorithm>

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
    }

private:

    void draw_about() override {
        draw_about_window(name(),Environment::instance()->get_value("version"));
    }

    void draw_object_properties() override {
        
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

        Attribute<bool> on_area_to_tilt(mesh_.facets.attributes(),"on_area_to_tilt");
        unsigned int nb_facet_on_areas_to_tilt = 0;
        FOR(f,mesh_.facets.nb()) { // for each facet
            vec3 normal = normalize(Geom::mesh_facet_normal(mesh_,f));
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
            if(std::abs(per_label_weights[5].first-per_label_weights[4].first) < 0.001) {
                // the 2 labels with the most weight are too close
                nb_facet_on_areas_to_tilt++;
                on_area_to_tilt[f] = true;
            }
            else {
                on_area_to_tilt[f] = false;
            }
        }

        fmt::println("{} facet(s) over {} are on an area to tilt",nb_facet_on_areas_to_tilt,mesh_.facets.nb());

        lighting_ = false;
        show_attributes_ = true;
        current_colormap_texture_ = colormaps_[COLORMAP_BLUE_RED].texture;
        attribute_min_ = 0.0f;
        attribute_max_ = 1.0f;
        attribute_ = "facets.on_area_to_tilt";
        attribute_name_ = "on_area_to_tilt";
        attribute_subelements_ = MESH_FACETS;

        return true;
    }
};

int main(int argc, char** argv) {
    AreasToTiltApp app;
	app.start(argc,argv);
    return 0;
}