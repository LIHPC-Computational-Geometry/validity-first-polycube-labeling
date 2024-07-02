#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "geometry_hexahedra.h"         // for compute_scaled_jacobian()
#include "SimpleMeshApplicationExt.h"   // for colormaps indices
#include "basic_stats.h"                // for BasicStats

using namespace GEO;

class HexMeshViewerApp : public SimpleMeshApplicationExt {
public:

    HexMeshViewerApp() : SimpleMeshApplicationExt("hex_mesh_viewer") {
        // change default values
        show_surface_ = false;
        show_volume_ = true;
        show_colored_cells_ = false;

        // init own variables
        show_SJ_ = false;
        drawing_settings_dirty_ = true;
    }

private:

    void draw_scene() override {
        if(drawing_settings_dirty_) {
            if(show_SJ_) {
                lighting_ = false;
                show_attributes_ = true;
                current_colormap_index_ = COLORMAP_PARULA;
                // reversed colormap:
                attribute_min_ = 1.0f; // SJ=1 -> blue
                attribute_max_ = 0.0f; // SJ=0 -> yellow
                attribute_ = "cells.SJ";
                attribute_name_ = "SJ";
                attribute_subelements_ = MESH_CELLS;
            }
            else {
                lighting_ = true;
                show_attributes_ = false;
            }
            drawing_settings_dirty_ = false;
        }
        SimpleMeshApplicationExt::draw_scene();
    }

    void draw_object_properties() override {
        ImGui::Checkbox("##MeshOnOff", &show_mesh_);
        ImGui::SameLine();
        ImGui::ColorEdit3WithPalette("Mesh", mesh_color_.data());
        ImGui::SliderFloat("Mesh size", &mesh_width_, 0.1f, 2.0f, "%.1f");
        if(ImGui::Checkbox("Cell color by SJ",&show_SJ_)) {
            drawing_settings_dirty_ = true;
        }
    }

    bool load(const std::string& filename) override {

        if(!SimpleMeshApplicationExt::load(filename)) { // if load failed
            show_SJ_ = false;
            drawing_settings_dirty_ = true;
            return false;
        }

        geo_assert(mesh_.cells.nb()!=0); // assert is a volume mesh

        mesh_.vertices.set_double_precision();
        
        compute_scaled_jacobian(mesh_,SJ_stats_);
		fmt::println(Logger::out("I/O"),"minSJ of this hex mesh is {}",SJ_stats_.min()); Logger::out("I/O").flush();
		
        show_SJ_ = true;
        drawing_settings_dirty_ = true;

        return true;
    }

private:

    bool show_SJ_;
    bool drawing_settings_dirty_;
    BasicStats SJ_stats_;
};

int main(int argc, char** argv) {
    HexMeshViewerApp app;
	app.start(argc,argv);
    return 0;
}