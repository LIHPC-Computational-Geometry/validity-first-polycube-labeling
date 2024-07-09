#pragma once

#include <geogram_gfx/gui/application.h>                // for set_style()
#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/file_system.h>                  // for is_file(), extension() in load()
#include <geogram/basic/command_line.h>                 // for get_arg_bool() in load()
#include <geogram/mesh/mesh_io.h>                       // for MeshIOFlags, mesh_load() in load()
#include <geogram/basic/command_line.h>                 // for CmdLine::get_arg() and CmdLine::set_arg()
#include <geogram/basic/string.h>                       // for String::string_ends_with()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>
#include <tuple>    // for std::tie()
#include <set>
#include <utility>  // for std::pair

#include "gui_base.h"					// for SimpleMeshApplicationExt
#include "labeling_graph.h"             // for LabelingGraph
#include "labeling.h"                   // for flip_labeling(), label2vector(), LABEL2STR(), compute_per_facet_fidelity
#include "labeling_io.h"                // for load_labeling(), save_labeling()
#include "labeling_generators.h"        // for naive_labeling()
#include "stats.h"                      // for IncrementalStats
#include "io_dump.h"                    // for dump_all_boundaries()

#define RED_WHITE_BLUE_LABELING_COLORS // red-white-blue instead of red-green-blue. better for (most) color-deficient users

#define LABELING_ATTRIBUTE_NAME "label" // name of the per-facet attribute in which the polycube labeling will be stored

// values for *_labeling_visu_mode_
#define VIEW_TRIANGLE_MESH		0
#define VIEW_RAW_LABELING		1
#define VIEW_LABELING_GRAPH		2
#define VIEW_FIDELITY			3
#define VIEW_INVALID_CHARTS		4
#define VIEW_INVALID_BOUNDARIES	5
#define VIEW_INVALID_CORNERS	6

// new colormaps
#define COLORMAP_LABELING       (SIMPLE_APPLICATION_NB_COLORMAPS)
#define COLORMAP_VALIDITY       ((SIMPLE_APPLICATION_NB_COLORMAPS)+1)

#define IMGUI_SLIDERS_WIDTH	100.0f

using namespace GEO;

class LabelingViewerApp : public SimpleMeshApplicationExt {
public:

    enum State {
        empty, // no file loaded
        triangle_mesh, // a mesh has been loaded
        labeling // a mesh has been loaded and a labeling has been loaded/computed
	};

    LabelingViewerApp(const std::string name = "labeling_viewer", bool auto_flip_normals = true);

	void ImGui_initialize() override;

protected:

	virtual void state_transition(State new_state); // virtual, so it can be overwritten in subclasses

	virtual void labeling_visu_mode_transition(int new_mode); // virtual, so it can be overwritten in subclasses

	void draw_scene() override;

	void draw_gui() override;

	void draw_menu_bar() override;

    void draw_object_properties() override;

    bool load(const std::string& filename) override;

	std::string supported_write_file_extensions() override;

	bool save(const std::string& filename) override;

	void GL_initialize() override;

	void mouse_button_callback(int button, int action, int mods, int source) override;

    void update_static_labeling_graph(bool allow_boundaries_between_opposite_labels);

protected:

	// variables to declare first, because first init in constructor
	ColorArray labeling_colors_; // per label colors
	ColorArray validity_colors_; // colors for valid and invalid labeling graph components

	// Core variables

	bool auto_flip_normals_; // if the app should ensure outward facet normals
	LabelingGraph lg_; // labeling seen as charts, boundaries and corners
	State state_; // should only be modified by state_transition()
	int labeling_visu_mode_; // not a enum, to be used in ImGui. should only be modified by labeling_visu_mode_transition() and ImGui::RadioButton
	// Additional geometry info
	std::vector<vec3> normals_; // facet normals
	std::vector<std::vector<index_t>> adj_facets_; // for each vertex, store adjacent facets. no ordering
	std::set<std::pair<index_t,index_t>> feature_edges_; // significant feature edges (with include/geometry.h > FEATURE_EDGES_MIN_ANGLE as threshold)

	// GUI

	bool show_ImGui_demo_window_; // if the ImGui demo window is currently displayed
	// Related to the triangle mesh
	bool show_normals_; // display facet normals as segments + vertices at the tips
	vec3f normals_color_; // color of the facet normals
	float normals_length_factor_; // length of the facet normals. multiplier of the normalized facet normal
	vec3 facet_center_; // last facet center computed
	vec3 normal_tip_; // last normal tip computed
	bool show_feature_edges_; // display feature edges
	int feature_edges_width_; // width of feature edges
	// Related to labelings
	bool allow_boundaries_between_opposite_labels_; // if boundaries between charts labeled to the same axis, but on > 180° angles, are allowed (when computing the labeling graph)
	bool show_boundaries_; // display labeling graph boundaries
	int boundaries_width_; // width of labeling graph boundaries
	std::size_t X_boundaries_group_index_; // index of the edges group in the overlay storing edges of boundaries mapped to the X axis
	std::size_t Y_boundaries_group_index_; // index of the edges group in the overlay storing edges of boundaries mapped to the Y axis
	std::size_t Z_boundaries_group_index_; // index of the edges group in the overlay storing edges of boundaries mapped to the Z axis
	std::size_t axisless_and_invalid_boundaries_group_index; // index of the edges group in the overlay storing axisless and invalid boundary edges
	std::size_t axisless_but_valid_boundaries_group_index_; // index of the edges group in the overlay storing axisless but valid boundary edges
	bool show_corners_; // display labeling graph corners
	vec3f corners_color_; // color of labeling graph corners
	float corners_size_; // size of labeling graph corners
	bool show_turning_points_; // display labeling graph turning-points
	vec3f turning_points_color_; // color of labeling graph turning-points
	float turning_points_size_; // size of labeling graph turning-points
	std::size_t turning_points_group_index_; // index of the points group in the overlay storing turning-points
	std::size_t valid_corners_group_index_; // index of the points group in the overlay storing valid corners
	std::size_t invalid_corners_group_index_;  // index of the points group in the overlay storing invalid corners
	std::string fidelity_text_label_; // text label describing per-facet fidelity stats
};

// print specialization of LabelingViewerApp::State for {fmt}
template <> struct fmt::formatter<LabelingViewerApp::State>: formatter<string_view> {
	// parse is inherited from formatter<string_view>.
	template <typename FormatContext>
	auto format(LabelingViewerApp::State state, FormatContext& ctx) -> decltype(ctx.out()) const {
		switch (state) {
		case LabelingViewerApp::State::empty			: return formatter<string_view>::format("empty", ctx);
		case LabelingViewerApp::State::triangle_mesh	: return formatter<string_view>::format("triangle_mesh", ctx);
		case LabelingViewerApp::State::labeling			: return formatter<string_view>::format("labeling", ctx);
		}
		return formatter<string_view>::format("?", ctx);
	}
};