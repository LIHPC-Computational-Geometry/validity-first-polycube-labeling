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

    // A list of points with the same color and size
    struct PointsGroup {
        std::vector<vec3> points;
        const float* color = nullptr; // storing a pointer allows for direct color modification from the GUI
        const float* size = nullptr;
        bool* show = nullptr;
        PointsGroup(const float* rgba, const float* size, bool* show) : color(rgba), size(size), show(show) {}
        void clear() { points.clear(); color = nullptr; size = nullptr; show = nullptr; }
    };

    // A list of edges with the color and width
    struct EdgesGroup {
        std::vector<std::pair<vec3,vec3>> edges;
        unsigned int colormap_index = 0;
        double texture_coordinate = 0.0;
        int* width = nullptr;
        bool* show = nullptr;
        EdgesGroup(unsigned int colormap_index, double texture_coordinate, int* width, bool* show) : colormap_index(colormap_index), texture_coordinate(texture_coordinate), width(width), show(show) {}
        void clear() { edges.clear(); colormap_index = 0; texture_coordinate = 0.0; width = nullptr; show = nullptr; }
    };

    // A set of colors that can be accessed both
    // as floats (for ImGui) and as chars (for OpenGL textures)
    class ColorArray {
    public:
        ColorArray(std::initializer_list<std::array<float,4>> init_values);

        float* as_floats();

        float* color_as_floats(std::size_t color_index);

        ImVec4 color_as_ImVec4(std::size_t color_index);

        unsigned char* as_chars();

        unsigned char* color_as_chars(std::size_t color_index);

        void update_chars_of_color(std::size_t color_index);
    
    private:
        std::vector<float> floats_;
        std::vector<unsigned char> chars_;
    };

    SimpleMeshApplicationExt(const std::string &name);

    ~SimpleMeshApplicationExt();

    void geogram_initialize(int argc, char** argv) override;

    void clear_scene_overlay();

    std::size_t new_points_group(const float* rgba, const float* size, bool* show);

    std::size_t new_edges_group(unsigned int colormap_index, double texture_coordinate, int* width, bool* show);

    void add_point_to_group(std::size_t group, double x, double y, double z);

    void add_edge_to_group(std::size_t group, double x1, double y1, double z1, double x2, double y2, double z2);

    void set_points_group_color(std::size_t index, const float* new_color);

    void set_edges_group_color(std::size_t index, unsigned int new_colormap_index, double new_texture_coordinate);

    void draw_scene() override;

    void cursor_pos_callback( double x, double y, int source ) override;

    index_t pick(MeshElementsFlags what);

    void init_rgba_colormap(const std::string& name, int width, int height, unsigned char * data);

    void update_GL_texture(unsigned int texture_index, int width, int height, unsigned char * data);

protected:

    vec2 cursor_pos_; // cursor position, in pixels from top-left corner

private:
    vector<PointsGroup> points_groups_;
    vector<EdgesGroup> edges_groups_;
    Memory::byte buffer[4]; // a 4-bytes buffer to read pixels
};