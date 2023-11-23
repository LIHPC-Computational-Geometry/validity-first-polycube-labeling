#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram_gfx/ImGui_ext/icon_font.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "LabelingViewerApp.h"
#include "CustomMeshHalfedges.h"

#define TEXT_GREEN  ImVec4(0.0f,0.8f,0.0f,1.0f)
#define TEXT_RED    ImVec4(0.8f,0.0f,0.0f,1.0f)
#define TEXT_GREY   ImVec4(0.6f,0.6f,0.6f,1.0f)

using namespace GEO;

class HalfedgeMovementsApp : public LabelingViewerApp {
public:

    HalfedgeMovementsApp() : LabelingViewerApp("halfedge_movements"), mesh_he(mesh_), custom_mesh_he(mesh_) {
        // change default values
        show_surface_ = true;
        surface_color_ = vec4f(0.9f,0.9f,0.9f,1.0f);
        show_volume_ = false;
        show_colored_cells_ = false;
        // init own variables
        halfedge.facet = NO_FACET;
        halfedge.corner = NO_CORNER;
        origin[0] = 0.0;
        origin[1] = 0.0;
        origin[2] = 0.0;
        extremity[0] = 0.0;
        extremity[1] = 0.0;
        extremity[2] = 0.0;
        rgba[0] = 1.0f;
        rgba[1] = 0.0f;
        rgba[2] = 0.0f;
        rgba[3] = 1.0f;
        ignore_borders_around_vertices = false;
    }

protected:

    void draw_scene() override {
        LabelingViewerApp::draw_scene();

        if(mesh_he.halfedge_is_valid(halfedge)) {
            glupSetColor4fv(GLUP_FRONT_COLOR, rgba); // use the color of this group of points
            // draw halfedge base vertex
            glupSetPointSize(10.0);
            glupBegin(GLUP_POINTS);
            glupPrivateVertex3dv(origin); // give a pointer to the coordinates to GLUP
            glupEnd();
            // draw halfedge vector
            glupBegin(GLUP_LINES);
            glupPrivateVertex3dv(origin);
            glupPrivateVertex3dv(extremity);
            glupEnd();
        }
    }

    void draw_object_properties() override {

        // display halfedge value

        ImGui::Text("halfedge.facet = {%d}",halfedge.facet);
        ImGui::Text("halfedge.corner = {%d}",halfedge.corner);

        if(mesh_he.halfedge_is_valid(halfedge)) {
            ImGui::Text("%s origin vertex is {%d}",icon_UTF8("star-of-life").c_str(),custom_mesh_he.halfedge_vertex_index_from(mesh_,halfedge));
            ImGui::Text("%s extremity vertex is {%d}",icon_UTF8("star-of-life").c_str(),custom_mesh_he.halfedge_vertex_index_to(mesh_,halfedge));
            ImGui::Text("%s left facet is {%d}",icon_UTF8("star-of-life").c_str(),custom_mesh_he.halfedge_facet_left(mesh_,halfedge));
            ImGui::Text("%s right facet is {%d}",icon_UTF8("star-of-life").c_str(),custom_mesh_he.halfedge_facet_right(mesh_,halfedge));
        }
        else {
            ImGui::Text("%s origin vertex is undefined",icon_UTF8("star-of-life").c_str());
            ImGui::Text("%s extremity vertex is undefined",icon_UTF8("star-of-life").c_str());
            ImGui::Text("%s left facet is undefined",icon_UTF8("star-of-life").c_str());
            ImGui::Text("%s right facet is undefined",icon_UTF8("star-of-life").c_str());
        }

        // display validity & on border

        if(mesh_he.halfedge_is_valid(halfedge)) {
            ImGui::TextColored(TEXT_GREEN,"Halfedge is valid");
            if(mesh_he.halfedge_is_border(halfedge)) { // if called on invalid halfedge -> bad assertion
                ImGui::TextColored(TEXT_GREEN,"Halfedge is on border");
            }
            else {
                ImGui::TextColored(TEXT_RED,"Halfedge is not on border");
            }
        }
        else {
            ImGui::TextColored(TEXT_RED,"Halfedge is invalid");
            ImGui::TextColored(TEXT_RED,"Halfedge is not on border");
        }

        // move around facets

        ImGui::SeparatorText("Around facets");
        if(ImGui::Button("Move to prev around facet")) {
            mesh_he.move_to_prev_around_facet(halfedge);
            update_points_coordinates();
        }
        if(ImGui::Button("Move to next around facet")) {
            mesh_he.move_to_next_around_facet(halfedge);
            update_points_coordinates();
        }

        // move around vertices

        ImGui::SeparatorText("Around vertices");
        if(ImGui::Button("Move to prev around vertex")) {
            mesh_he.move_to_prev_around_vertex(halfedge);
            update_points_coordinates();
        }
        if(ImGui::Button("Move to next around vertex")) {
            mesh_he.move_to_next_around_vertex(halfedge);
            update_points_coordinates();
        }
        ImGui::TextUnformatted(icon_UTF8("star-of-life").c_str());
        ImGui::SameLine();
        ImGui::Checkbox("ignore borders",&ignore_borders_around_vertices);
        if(ImGui::Button(icon_UTF8("star-of-life") + " Custom clockwise around vertex")) {
            custom_mesh_he.move_clockwise_around_vertex(halfedge,ignore_borders_around_vertices);
            update_points_coordinates();
        }
        if(ImGui::Button(icon_UTF8("star-of-life") + "Custom counterclockwise around vertex")) {
            custom_mesh_he.move_counterclockwise_around_vertex(halfedge,ignore_borders_around_vertices);
            update_points_coordinates();
        }

        // move around borders

        ImGui::BeginDisabled(!mesh_he.halfedge_is_valid(halfedge) || !mesh_he.halfedge_is_border(halfedge));
        ImGui::SeparatorText("Around borders");
        if(ImGui::Button("Move to prev around border")) {
            mesh_he.move_to_prev_around_border(halfedge);
            update_points_coordinates();
        }
        if(ImGui::Button("Move to next around border")) {
            mesh_he.move_to_next_around_border(halfedge);
            update_points_coordinates();
        }
        ImGui::EndDisabled();

        // flip halfedge

        ImGui::SeparatorText("Flip");
        if(ImGui::Button("Move to opposite")) {
            mesh_he.move_to_opposite(halfedge);
            update_points_coordinates();
        }

        // legend

        ImGui::Separator();
        ImGui::TextColored(TEXT_GREY,"%s = not in vanilla Geogram",icon_UTF8("star-of-life").c_str());
    }

    bool load(const std::string& filename) override {
        if(LabelingViewerApp::load(filename)) {
            mesh_.vertices.set_double_precision();
            if(!mesh_he.halfedge_is_valid(halfedge)) {
                halfedge.facet = 0;
                halfedge.corner = mesh_.facets.corner(0,0); // first corner of facet 0
            }
            update_points_coordinates();
            if(state_ == labeling) {
                labeling_visu_mode_transition(VIEW_RAW_LABELING);
                mesh_he.set_use_facet_region(std::string(LABELING_ATTRIBUTE_NAME));
                custom_mesh_he.set_use_facet_region(std::string(LABELING_ATTRIBUTE_NAME));
            }
            return true;
        }
        halfedge.facet = NO_FACET;
        halfedge.corner = NO_CORNER;
        return false;
    }

    void update_points_coordinates() {
        if(mesh_he.halfedge_is_valid(halfedge)) {
            origin[0] = Geom::halfedge_vertex_from(mesh_,halfedge)[0];
            origin[1] = Geom::halfedge_vertex_from(mesh_,halfedge)[1];
            origin[2] = Geom::halfedge_vertex_from(mesh_,halfedge)[2];
            extremity[0] = Geom::halfedge_vertex_to(mesh_,halfedge)[0];
            extremity[1] = Geom::halfedge_vertex_to(mesh_,halfedge)[1];
            extremity[2] = Geom::halfedge_vertex_to(mesh_,halfedge)[2];
        }
    }

    MeshHalfedges mesh_he;
    CustomMeshHalfedges custom_mesh_he;
    MeshHalfedges::Halfedge halfedge;
    double origin[3];
    double extremity[3];
    float rgba[4]; // halfedge color, red-green-blue-alpha
    bool ignore_borders_around_vertices;
};

int main(int argc, char** argv) {
    HalfedgeMovementsApp app;
	app.start(argc,argv);
    return 0;
}