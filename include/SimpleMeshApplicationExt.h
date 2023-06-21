// Adding functionalities to SimpleMeshApplication:
// - overlay of points and edges, stored in groups (all elements in a group have the same color)
// - modification of an OpenGL texture (to update a colormap from GUI)

#pragma once

#include <geogram_gfx/gui/application.h> // for Application::Gl_initialize()
#include <geogram_gfx/gui/simple_application.h>
#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/vecg.h>

#include <vector>
#include <initializer_list>

#include "containers.h"

// Indices of SimpleApplication::colormaps_
// Same ordering as in ext/geogram/src/lib/geogram_gfx/gui/simple_application.cpp init_colormaps()
#define COLORMAP_FRENCH			0
#define COLORMAP_BLACK_WHITE	1
#define COLORMAP_VIRIDIS		2
#define COLORMAP_RAINBOW		3
#define COLORMAP_CEI_60757		4
#define COLORMAP_INFERNO		5
#define COLORMAP_MAGMA			6
#define COLORMAP_PARULA			7
#define COLORMAP_PLASMA			8
#define COLORMAP_BLUE_RED		9

#define SIMPLE_APPLICATION_NB_COLORMAPS ((COLORMAP_BLUE_RED)+1)

using namespace GEO;

class SimpleMeshApplicationExt : public SimpleMeshApplication {
public:

    struct PointsGroup {
        std::vector<vec3> points;
        const float* color = nullptr; // storing a pointer allows for direct color modification from the GUI
        bool show = false;
        PointsGroup(const float* rgba, bool show) : color(rgba), show(show) {}
        void clear() { points.clear(); color = nullptr; show = false; }
    };

    struct EdgesGroup {
        std::vector<std::pair<vec3,vec3>> edges;
        const float* color = nullptr; // storing a pointer allows for direct color modification from the GUI
        bool show = false;
        EdgesGroup(const float* rgba, bool show) : color(rgba), show(show) {}
        void clear() { edges.clear(); color = nullptr; show = false; }
    };

    SimpleMeshApplicationExt(const std::string &name) : SimpleMeshApplication(name) {}

    ~SimpleMeshApplicationExt() {
        clear_scene_overlay();
    }

    void start(int argc, char** argv) override {
        SimpleMeshApplication::start(argc,argv);
    }

    void clear_scene_overlay() {
        points_groups_.clear();
        edges_groups_.clear();
    }

    std::size_t new_points_group(const float* rgba, bool show) {
        points_groups_.push_back(PointsGroup(rgba,show));
        return index_of_last(points_groups_);
    }

    std::size_t new_edges_group(const float* rgba, bool show) {
        edges_groups_.push_back(EdgesGroup(rgba,show));
        return index_of_last(edges_groups_);
    }

    void add_point_to_group(std::size_t group, double x, double y, double z) {
        points_groups_.at(group).points.push_back({x,y,z});
    }

    void add_edge_to_group(std::size_t group, double x1, double y1, double z1, double x2, double y2, double z2) {
        edges_groups_.at(group).edges.push_back(
            std::make_pair<vec3,vec3>({x1,y1,z1},{x2,y2,z2})
        );
    }

    void set_points_group_color(std::size_t index, const float* new_color) {
        points_groups_.at(index).color = new_color;
    }

    void set_edges_group_color(std::size_t index, const float* new_color) {
        edges_groups_.at(index).color = new_color;
    }

    void set_points_group_visibility(std::size_t index, bool visible) {
        points_groups_.at(index).show = visible;
    }

    void set_edges_group_visibility(std::size_t index, bool visible) {
        edges_groups_.at(index).show = visible;
    }

    void points_groups_show_only(std::initializer_list<std::size_t> indices) {
        // hide all groups
        for(PointsGroup& group : points_groups_) {
            group.show = false;
        }
        // show selected groups
        for(std::size_t group_index : indices) {
            points_groups_.at(group_index).show = true;
        }
    }

    void edges_groups_show_only(std::initializer_list<std::size_t> indices) {
        // hide all groups
        for(EdgesGroup& group : edges_groups_) {
            group.show = false;
        }
        // show selected groups
        for(std::size_t group_index : indices) {
            edges_groups_.at(group_index).show = true;
        }
    }

    void draw_scene() override {
        SimpleMeshApplication::draw_scene();

        // draw points groups
        for(const PointsGroup& group : points_groups_) {
            if(group.show) {
                geo_assert(group.color != nullptr);
                glupSetColor4fv(GLUP_FRONT_COLOR, group.color); // use the color of this group of points
                glupSetPointSize(10.0f);
                for(std::size_t i = 0; i<group.points.size(); ++i) { // for each point in this group
                    glupBegin(GLUP_POINTS);
                    glupPrivateVertex3dv(group.points[i].data()); // give a pointer to the coordinates to GLUP
                    glupEnd();
                }
            }
        }

        // draw edges groups
        for(const auto& group : edges_groups_) {
            if(group.show) {
                geo_assert(group.color != nullptr);
                glupSetColor4fv(GLUP_FRONT_COLOR, group.color); // use the color of this group of edges
                glupBegin(GLUP_LINES);
                for(const auto& edge : group.edges) { // for each edge in this group
                    glupPrivateVertex3dv(edge.first.data()); // give a pointer to the coordinates to GLUP
                    glupPrivateVertex3dv(edge.second.data()); // give a pointer to the coordinates to GLUP
                }
                glupEnd();
            }
        }
    }

    void init_rgba_colormap(const std::string& name, int width, int height, unsigned char * data) {
        // like SimpleApplication::init_colormap() but without XPM data, just an u8 array of RGBA values
        colormaps_.push_back(ColormapInfo());
        colormaps_.rbegin()->name = name;
        glGenTextures(1, &colormaps_.rbegin()->texture); // OpenGL will create a new texture and write its id in texture
        glBindTexture(GL_TEXTURE_2D, colormaps_.rbegin()->texture);
        glTexImage2D(
			GL_TEXTURE_2D, 0,
			GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data
		);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
        );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void update_GL_texture(unsigned int texture_index, int width, int height, unsigned char * data) {
        // similar to init_rgba_colormap()
        glBindTexture(GL_TEXTURE_2D, texture_index);
        glTexImage2D(
            GL_TEXTURE_2D, 0,
            GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data
        );
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
        );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

private:
    vector<PointsGroup> points_groups_;
    vector<EdgesGroup> edges_groups_;
};