#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/vecg.h>
#include <geogram/mesh/mesh_halfedges.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "SimpleMeshApplicationExt.h"   // for colormaps indices

using namespace GEO;

class HalfedgeMovementsApp : public SimpleMeshApplicationExt {
public:

    HalfedgeMovementsApp() : SimpleMeshApplicationExt("halfedge_movements"), mesh_he(mesh_) {
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
        rgba[0] = 1.0f;
        rgba[1] = 0.0f;
        rgba[2] = 0.0f;
        rgba[3] = 1.0f;
    }

protected:

    void draw_scene() override {
        SimpleMeshApplication::draw_scene();

        glupSetColor4fv(GLUP_FRONT_COLOR, rgba); // use the color of this group of points
        // draw halfedge base vertex
        glupSetPointSize(10.0);
        glupBegin(GLUP_POINTS);
        glupPrivateVertex3dv(origin); // give a pointer to the coordinates to GLUP
        glupEnd();
    }

    void draw_object_properties() override {
        ImGui::SeparatorText("Around facets");
        if(ImGui::Button("Move to prev around facet")) {
            mesh_he.move_to_prev_around_facet(halfedge);
            origin[0] = Geom::halfedge_vertex_from(mesh_,halfedge)[0];
            origin[1] = Geom::halfedge_vertex_from(mesh_,halfedge)[1];
            origin[2] = Geom::halfedge_vertex_from(mesh_,halfedge)[2];
        }
        if(ImGui::Button("Move to next around facet")) {
            mesh_he.move_to_next_around_facet(halfedge);
            origin[0] = Geom::halfedge_vertex_from(mesh_,halfedge)[0];
            origin[1] = Geom::halfedge_vertex_from(mesh_,halfedge)[1];
            origin[2] = Geom::halfedge_vertex_from(mesh_,halfedge)[2];
        }
        ImGui::SeparatorText("Around vertices");
        if(ImGui::Button("Move to prev around vertex")) {
            mesh_he.move_to_prev_around_vertex(halfedge);
            origin[0] = Geom::halfedge_vertex_from(mesh_,halfedge)[0];
            origin[1] = Geom::halfedge_vertex_from(mesh_,halfedge)[1];
            origin[2] = Geom::halfedge_vertex_from(mesh_,halfedge)[2];
        }
        if(ImGui::Button("Move to next around vertex")) {
            mesh_he.move_to_next_around_vertex(halfedge);
            origin[0] = Geom::halfedge_vertex_from(mesh_,halfedge)[0];
            origin[1] = Geom::halfedge_vertex_from(mesh_,halfedge)[1];
            origin[2] = Geom::halfedge_vertex_from(mesh_,halfedge)[2];
        }
    }

    bool load(const std::string& filename) override {
        if(SimpleMeshApplication::load(filename)) {
            mesh_.vertices.set_double_precision();
            halfedge.facet = 0;
            halfedge.corner = mesh_.facets.corner(0,0); // first corner of facet 0
            return true;
        }
        return false;
    }

    MeshHalfedges mesh_he;
    MeshHalfedges::Halfedge halfedge;
    double origin[3];
    float rgba[4]; // halfedge color, red-green-blue-alpha
};

int main(int argc, char** argv) {
    HalfedgeMovementsApp app;
	app.start(argc,argv);
    return 0;
}