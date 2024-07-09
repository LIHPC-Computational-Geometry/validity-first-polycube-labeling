// Adding functionalities to SimpleMeshApplication:
// - overlay of points and edges, stored in groups (all elements in a group have the same color)
// - modification of an OpenGL texture (to update a colormap from GUI)
// - pick() method, to retrieve the index of the element under the cursor

#pragma once

#include <geogram_gfx/gui/application.h> // for Application::Gl_initialize()
#include <geogram_gfx/gui/simple_application.h>
#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/vecg.h>

#include <geogram_gfx/third_party/imgui/imgui.h> // for ImVec4

#include <vector>
#include <array>
#include <initializer_list>

#include "containers_std.h" // for index_of_last()

// Indices of SimpleApplication::colormaps_
// Same ordering as in ext/geogram/src/lib/geogram_gfx/gui/simple_application.cpp init_colormaps()
// XPM files are in ext/geogram/src/lib/geogram_gfx/gui/colormaps
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
        const float* size = nullptr;
        bool* show = nullptr;
        PointsGroup(const float* rgba, const float* size, bool* show) : color(rgba), size(size), show(show) {}
        void clear() { points.clear(); color = nullptr; size = nullptr; show = nullptr; }
    };

    struct EdgesGroup {
        std::vector<std::pair<vec3,vec3>> edges;
        unsigned int colormap_index = 0;
        double texture_coordinate = 0.0;
        int* width = nullptr;
        bool* show = nullptr;
        EdgesGroup(unsigned int colormap_index, double texture_coordinate, int* width, bool* show) : colormap_index(colormap_index), texture_coordinate(texture_coordinate), width(width), show(show) {}
        void clear() { edges.clear(); colormap_index = 0; texture_coordinate = 0.0; width = nullptr; show = nullptr; }
    };

    // both float (for ImGui) and char (for OpenGL textures) representations
    class ColorArray {
    public:
        ColorArray(std::initializer_list<std::array<float,4>> init_values) {
            chars_.resize(init_values.size()*4);
            floats_.resize(init_values.size()*4);
            std::size_t color_index = 0;
            for(auto &color : init_values) {
                floats_[4*color_index+0] = color[0]; // R
                floats_[4*color_index+1] = color[1]; // G
                floats_[4*color_index+2] = color[2]; // B
                floats_[4*color_index+3] = color[3]; // A
                update_chars_of_color(color_index);
                color_index++;
            }
        }

        float* as_floats() {
            return floats_.data();
        }

        float* color_as_floats(std::size_t color_index) {
            return (floats_.data()+(4*color_index));
        }

        ImVec4 color_as_ImVec4(std::size_t color_index) {
            return ImVec4(
                floats_.data()[4*color_index+0],
                floats_.data()[4*color_index+1],
                floats_.data()[4*color_index+2],
                floats_.data()[4*color_index+3]
            );
        }

        unsigned char* as_chars() {
            return chars_.data();
        }

        unsigned char* color_as_chars(std::size_t color_index) {
            return (chars_.data()+(4*color_index));
        }

        void update_chars_of_color(std::size_t color_index) {
            chars_[4*color_index+0] = static_cast<unsigned char>(floats_[4*color_index+0]*255.0f); // R
            chars_[4*color_index+1] = static_cast<unsigned char>(floats_[4*color_index+1]*255.0f); // G
            chars_[4*color_index+2] = static_cast<unsigned char>(floats_[4*color_index+2]*255.0f); // B
            chars_[4*color_index+3] = static_cast<unsigned char>(floats_[4*color_index+3]*255.0f); // A
        }
    private:
        std::vector<float> floats_;
        std::vector<unsigned char> chars_;
    };

    SimpleMeshApplicationExt(const std::string &name) : SimpleMeshApplication(name) {}

    ~SimpleMeshApplicationExt() {
        clear_scene_overlay();
    }

    // Change window size
    // by order of priority:
    // - command line argument "gfx:geometry"
    // - preprocessor macro WINDOW_SIZE defined in CMakeLists.txt
    // - default value "1024x1024", hard-coded in ext/geogram/src/lib/geogram/basic/command_line_args.cpp import_arg_group_gfx()
#ifdef WINDOW_SIZE
    void geogram_initialize(int argc, char** argv) override {
		SimpleMeshApplication::geogram_initialize(argc,argv);
		if(!phone_screen_ && CmdLine::get_arg("gfx:geometry") == "1024x1024") { // if not a phone screen and default value for gfx:geometry 

            CmdLine::set_arg("gfx:geometry", WINDOW_SIZE); // bigger window

		}
	}
#endif

    void clear_scene_overlay() {
        points_groups_.clear();
        edges_groups_.clear();
    }

    std::size_t new_points_group(const float* rgba, const float* size, bool* show) {
        points_groups_.push_back(PointsGroup(rgba,size,show));
        return index_of_last(points_groups_);
    }

    std::size_t new_edges_group(unsigned int colormap_index, double texture_coordinate, int* width, bool* show) {
        edges_groups_.push_back(EdgesGroup(colormap_index,texture_coordinate,width,show));
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

    void set_edges_group_color(std::size_t index, unsigned int new_colormap_index, double new_texture_coordinate) {
        edges_groups_.at(index).colormap_index = new_colormap_index;
        edges_groups_.at(index).texture_coordinate = new_texture_coordinate;
    }

    void draw_scene() override {
        SimpleMeshApplication::draw_scene();

        // draw points groups
        for(const PointsGroup& group : points_groups_) {
            if(*group.show) {
                geo_assert(group.color != nullptr);
                glupSetColor4fv(GLUP_FRONT_COLOR, group.color); // use the color of this group of points
                glupSetPointSize(*group.size);
                for(std::size_t i = 0; i<group.points.size(); ++i) { // for each point in this group
                    glupBegin(GLUP_POINTS);
                    glupPrivateVertex3dv(group.points[i].data()); // give a pointer to the coordinates to GLUP
                    glupEnd();
                }
            }
        }

        // draw edges groups
        for(const auto& group : edges_groups_) {
            if(*group.show) {
                glupEnable(GLUP_TEXTURING);
                glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
                glBindTexture(GL_TEXTURE_2D, (GLuint) colormaps_[group.colormap_index].texture);
                glupTextureType(GLUP_TEXTURE_2D);
                glupTextureMode(GLUP_TEXTURE_REPLACE);
                glupPrivateTexCoord1d(group.texture_coordinate);
                glupSetMeshWidth(*group.width);
                glupBegin(GLUP_LINES);
                for(const auto& edge : group.edges) { // for each edge in this group
                    glupPrivateVertex3dv(edge.first.data()); // give a pointer to the coordinates to GLUP
                    glupPrivateVertex3dv(edge.second.data()); // give a pointer to the coordinates to GLUP
                }
                glupDisable(GLUP_TEXTURING);
                glupEnd();
            }
        }
    }

    void cursor_pos_callback( double x, double y, int source ) override {
		cursor_pos_ = GEO::vec2(x,y); // store locally the current cursor position
		SimpleMeshApplication::cursor_pos_callback(x,y,source);
	}

    index_t pick(MeshElementsFlags what) {

        // see https://github.com/BrunoLevy/geogram/discussions/88
        // see https://github.com/BrunoLevy/geogram/pull/102
        // based on https://github.com/BrunoLevy/GraphiteThree/blob/main/src/lib/OGF/renderer/context/rendering_context.cpp get_picked_point()

        index_t x = index_t(cursor_pos_.x), y = index_t(cursor_pos_.y); // double to integer conversion of current cursor position
        if(x >= get_width() || y >= get_height()) { // if cursor out of the window
            return index_t(-1);
        }
        y = get_height()-1-y; // change Y axis orientation. glReadPixels() wants pixel coordinates from bottom-left corner
        mesh_gfx()->set_picking_mode(what); // instead of rendering colors, mesh_gfx will render indices
        draw_scene(); // rendering
        // read the index of the picked element using glReadPixels()
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glPixelStorei(GL_PACK_ROW_LENGTH, 1);
        glReadPixels(
            GLint(x),GLint(y),1,1,GL_RGBA,GL_UNSIGNED_BYTE,buffer
        );
        mesh_gfx()->set_picking_mode(MESH_NONE); // go back to color rendering mode
        // decode index from pixel color
        return index_t(buffer[0])        |
              (index_t(buffer[1]) << 8)  |
              (index_t(buffer[2]) << 16) |
              (index_t(buffer[3]) << 24);
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

protected:
    vec2 cursor_pos_; // cursor position, in pixels from top-left corner

private:
    vector<PointsGroup> points_groups_;
    vector<EdgesGroup> edges_groups_;
    Memory::byte buffer[4]; // a 4-bytes buffer to read pixels
};