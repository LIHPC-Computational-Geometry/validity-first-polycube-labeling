#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram_gfx/imgui_ext/icon_font.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "gui_labeling.h"
#include "geometry_halfedges.h"
#include "geometry.h"               // for facet_normals_are_inward() & flip_facet_normals()

#define TEXT_GREEN  ImVec4(0.0f,0.8f,0.0f,1.0f)
#define TEXT_RED    ImVec4(0.8f,0.0f,0.0f,1.0f)
#define TEXT_GREY   ImVec4(0.6f,0.6f,0.6f,1.0f)

using namespace GEO;

class HalfedgeMovementsApp : public LabelingViewerApp {
public:

    HalfedgeMovementsApp()
    : LabelingViewerApp("halfedge_movements",false), // do not auto flip normals (in order to know the orientation in the file)
    mesh_he(mesh_),
    mesh_he_ext(mesh_) {
        // change default values
        show_surface_ = true;
        surface_color_ = vec4f(0.9f,0.9f,0.9f,1.0f);
        show_volume_ = false;
        show_colored_cells_ = false;
        // init own variables
        halfedge.facet = NO_FACET;
        halfedge.corner = NO_CORNER;
        origin_vertex = index_t(-1);
        origin[0] = 0.0;
        origin[1] = 0.0;
        origin[2] = 0.0;
        extremity_vertex = index_t(-1);
        extremity[0] = 0.0;
        extremity[1] = 0.0;
        extremity[2] = 0.0;
        left_facet = index_t(-1);
        right_facet = index_t(-1);
        rgba[0] = 1.0f;
        rgba[1] = 0.0f;
        rgba[2] = 0.0f;
        rgba[3] = 1.0f;
        ignore_borders_around_vertices = false;
        facet_normals_inward = true; // updated in (load)
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
            // draw arrow head
            if(left_facet != NO_FACET) {
                facet_center = Geom::mesh_facet_center(mesh_,left_facet);
                glupBegin(GLUP_LINES);
                glupPrivateVertex3dv(extremity);
                glupPrivateVertex3dv(facet_center.data());
                glupEnd();
            }
            if(right_facet != NO_FACET) {
                facet_center = Geom::mesh_facet_center(mesh_,right_facet);
                glupBegin(GLUP_LINES);
                glupPrivateVertex3dv(extremity);
                glupPrivateVertex3dv(facet_center.data());
                glupEnd();
            }
        }
    }

    void draw_object_properties() override {

        // facet normal direction & button to flip them

        ImGui::SeparatorText("Direction of facet normals");
        if(facet_normals_inward) {
            ImGui::TextColored(TEXT_RED,"Inward");
        }
        else {
            ImGui::TextColored(TEXT_GREEN,"Outward");
        }
        ImGui::SameLine();
        if(ImGui::Button("Flip normals")) {
            flip_facet_normals(mesh_);
            facet_normals_inward = facet_normals_are_inward(mesh_);
        }

        // halfedge color modification

        ImGui::SeparatorText("Halfedge and its neighborhood");
        ImGui::ColorEdit4WithPalette("halfedge color",rgba);

        // display halfedge value

        ImGui::Text("halfedge.facet = {%d}",halfedge.facet);
        ImGui::Text("halfedge.corner = {%d}",halfedge.corner);

        if(mesh_he.halfedge_is_valid(halfedge)) {
            ImGui::Text("%s origin vertex is {%d}",icon_UTF8("star-of-life").c_str(),origin_vertex);
            ImGui::Text("%s extremity vertex is {%d}",icon_UTF8("star-of-life").c_str(),extremity_vertex);
            ImGui::Text("%s left facet is {%d}",icon_UTF8("star-of-life").c_str(),left_facet);
            ImGui::Text("%s right facet is {%d}",icon_UTF8("star-of-life").c_str(),right_facet);
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
            update_geometry_data();
        }
        if(ImGui::Button("Move to next around facet")) {
            mesh_he.move_to_next_around_facet(halfedge);
            update_geometry_data();
        }

        // move around vertices

        ImGui::SeparatorText("Around vertices");
        if(ImGui::Button("Move to prev around vertex")) {
            mesh_he.move_to_prev_around_vertex(halfedge);
            update_geometry_data();
        }
        if(ImGui::Button("Move to next around vertex")) {
            mesh_he.move_to_next_around_vertex(halfedge);
            update_geometry_data();
        }
        ImGui::TextUnformatted(icon_UTF8("star-of-life").c_str());
        ImGui::SameLine();
        ImGui::Checkbox("ignore borders",&ignore_borders_around_vertices);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Allow \"move (counter)clockwise around vertex\" functions to cross borders");
        if(ImGui::Button(icon_UTF8("star-of-life") + "Move clockwise around vertex")) {
            mesh_he_ext.move_clockwise_around_vertex(halfedge,ignore_borders_around_vertices);
            update_geometry_data();
        }
        if(ImGui::Button(icon_UTF8("star-of-life") + "Move counterclockwise around vertex")) {
            mesh_he_ext.move_counterclockwise_around_vertex(halfedge,ignore_borders_around_vertices);
            update_geometry_data();
        }

        // move around borders

        ImGui::BeginDisabled(!mesh_he.halfedge_is_valid(halfedge) || !mesh_he.halfedge_is_border(halfedge));
        ImGui::SeparatorText("Around borders");
        if(ImGui::Button("Move to prev around border")) {
            mesh_he.move_to_prev_around_border(halfedge);
            update_geometry_data();
        }
        if(ImGui::Button("Move to next around border")) {
            mesh_he.move_to_next_around_border(halfedge);
            update_geometry_data();
        }
        ImGui::EndDisabled();

        // flip halfedge

        ImGui::SeparatorText("Flip");
        if(ImGui::Button("Move to opposite")) {
            mesh_he.move_to_opposite(halfedge);
            update_geometry_data();
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
            update_geometry_data();
            if(state_ == labeling) {
                labeling_visu_mode_transition(VIEW_RAW_LABELING);
                mesh_he.set_use_facet_region(std::string(LABELING_ATTRIBUTE_NAME));
                mesh_he_ext.set_use_facet_region(std::string(LABELING_ATTRIBUTE_NAME));
            }
            facet_normals_inward = facet_normals_are_inward(mesh_);
            return true;
        }
        halfedge.facet = NO_FACET;
        halfedge.corner = NO_CORNER;
        return false;
    }

    void update_geometry_data() {
        if(mesh_he.halfedge_is_valid(halfedge)) {
            origin_vertex = Geom::halfedge_vertex_index_from(mesh_,halfedge);
            origin[0] = Geom::halfedge_vertex_from(mesh_,halfedge)[0];
            origin[1] = Geom::halfedge_vertex_from(mesh_,halfedge)[1];
            origin[2] = Geom::halfedge_vertex_from(mesh_,halfedge)[2];
            extremity_vertex = Geom::halfedge_vertex_index_to(mesh_,halfedge);
            extremity[0] = Geom::halfedge_vertex_to(mesh_,halfedge)[0];
            extremity[1] = Geom::halfedge_vertex_to(mesh_,halfedge)[1];
            extremity[2] = Geom::halfedge_vertex_to(mesh_,halfedge)[2];
            left_facet = Geom::halfedge_facet_left(mesh_,halfedge);
            right_facet = Geom::halfedge_facet_right(mesh_,halfedge);
        }
    }

    MeshHalfedges mesh_he;
    MeshHalfedgesExt mesh_he_ext;
    MeshHalfedges::Halfedge halfedge;
    index_t origin_vertex;
    double origin[3];
    index_t extremity_vertex;
    double extremity[3];
    index_t left_facet;
    index_t right_facet;
    float rgba[4]; // halfedge color, red-green-blue-alpha
    bool ignore_borders_around_vertices;
    vec3 facet_center;
    bool facet_normals_inward;
};

int main(int argc, char** argv) {
    HalfedgeMovementsApp app;
	app.start(argc,argv);
    return 0;
}