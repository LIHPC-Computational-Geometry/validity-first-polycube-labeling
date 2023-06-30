#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/file_system.h> // for is_file(), extension() in load()
#include <geogram/basic/command_line.h> // for get_arg_bool() in load()
#include <geogram/mesh/mesh_io.h> // for MeshIOFlags, mesh_load() in load()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "hex_mesh.h"                   // for compute_scaled_jacobian()
#include "SimpleMeshApplicationExt.h"   // for colormaps indices
#include "about_window.h"               // for draw_about_window()

using namespace GEO;

class HexMeshViewerApp : private SimpleMeshApplication {
public:

    HexMeshViewerApp() : SimpleMeshApplication("hex_mesh_viewer") {
        // change default values
        show_surface_ = false;
        show_volume_ = true;
        show_colored_cells_ = false;

        // init own variables
        show_SJ_ = false;
        drawing_settings_dirty_ = true;
    }

    void start(int argc, char** argv) override {
        SimpleMeshApplication::start(argc,argv);
    }

private:

    void draw_about() override {
        draw_about_window(name(),Environment::instance()->get_value("version"));
    }

    void draw_scene() override {
        if(drawing_settings_dirty_) {
            if(show_SJ_) {
                lighting_ = false;
                show_attributes_ = true;
                current_colormap_texture_ = colormaps_[COLORMAP_PARULA].texture;
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
        SimpleMeshApplication::draw_scene();
    }

    void draw_object_properties() override {
        if(ImGui::Checkbox("Cell color by SJ",&show_SJ_)) {
            drawing_settings_dirty_ = true;
        }
    }

    bool load(const std::string& filename) override {

        //// based on SimpleMeshApplication::load() ////////////////////////////////////////////////////////////////

        if(!FileSystem::is_file(filename)) {
            Logger::out("I/O") << "is not a file" << std::endl;
        }
        mesh_gfx_.set_mesh(nullptr);
        mesh_.clear(true,false); // keep_attributes=true is very important. else runtime error when re-importing files
        if(GEO::CmdLine::get_arg_bool("single_precision")) {
            mesh_.vertices.set_single_precision();
        }
        MeshIOFlags flags;
        if(!mesh_load(filename, mesh_, flags)) {
            show_SJ_ = false;               // added
            drawing_settings_dirty_ = true; // added
            return false;
        }
        mesh_gfx_.set_animate(false);            
        mesh_.vertices.set_dimension(3);
        double xyzmin[3];
        double xyzmax[3];
        get_bbox(mesh_, xyzmin, xyzmax, false);
        set_region_of_interest(
            xyzmin[0], xyzmin[1], xyzmin[2],
            xyzmax[0], xyzmax[1], xyzmax[2]
        );
        mesh_.vertices.set_double_precision(); // added
        mesh_gfx_.set_mesh(&mesh_);
	    current_file_ = filename;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        double minSJ = compute_scaled_jacobian(mesh_);
		fmt::println(Logger::out("I/O"),"minSJ of this hex mesh is {}",minSJ); Logger::out("I/O").flush();
		show_SJ_ = true;
        drawing_settings_dirty_ = true;

        return true;
    }

private:

    bool show_SJ_;
    bool drawing_settings_dirty_;
};

int main(int argc, char** argv) {
    HexMeshViewerApp app;
	app.start(argc,argv);
    return 0;
}