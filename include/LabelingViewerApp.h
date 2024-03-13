#pragma once

#include <geogram_gfx/gui/application.h> // for set_style()
#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/file_system.h> // for is_file(), extension() in load()
#include <geogram/basic/command_line.h> // for get_arg_bool() in load()
#include <geogram/mesh/mesh_io.h> // for MeshIOFlags, mesh_load() in load()
#include <geogram/basic/command_line.h> // for CmdLine::get_arg() and CmdLine::set_arg()
#include <geogram/basic/string.h> // for String::string_ends_with()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>
#include <tuple>	// for std::tie()
#include <set>
#include <utility>	// for std::pair

#include "SimpleMeshApplicationExt.h"
#include "LabelingGraph.h"   // for StaticLabelingGraph
#include "labeling.h"		 // for load_labeling(), naive_labeling(), save_labeling()
#include "basic_stats.h"	 // for BasicStats
#include "dump_mesh.h"		 // for dump_all_boundaries()

#define RED_WHITE_BLUE_LABELING_COLORS // better for (most) color-deficient users

#define LABELING_ATTRIBUTE_NAME "label"

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

using namespace GEO;

const float green[4] = {0.0f, 1.0f, 0.0f, 1.0f};
const float dark_blue[4] = {0.0f, 0.0f, 5.0f, 1.0f};

class LabelingViewerApp : public SimpleMeshApplicationExt {
public:

    enum State {
        empty,
        triangle_mesh,
        labeling
	};

    LabelingViewerApp(const std::string name = "labeling_viewer", bool auto_flip_normals = true) : SimpleMeshApplicationExt(name),
						     labeling_colors_({
								{1.0f, 0.0f, 0.0f, 1.0f}, // label 0 -> red
								{0.6f, 0.0f, 0.0f, 1.0f}, // label 1 -> darker red
							#ifdef RED_WHITE_BLUE_LABELING_COLORS
								{1.0f, 1.0f, 1.0f, 1.0f}, // label 2 -> white
								{0.6f, 0.6f, 0.6f, 1.0f}, // label 3 -> grey
							#else
								{0.0f, 1.0f, 0.0f, 1.0f}, // label 2 -> green
								{0.0f, 0.6f, 0.0f, 1.0f}, // label 3 -> darker green
							#endif
								{0.0f, 0.0f, 1.0f, 1.0f}, // label 4 -> blue
								{0.0f, 0.0f, 0.6f, 1.0f}  // label 5 -> darker blue
						   	 }),
							 validity_colors_({
								{1.0f, 0.25f, 0.25f, 1.0f}, // invalid charts/boundaries/corners in light red
								{0.25f, 0.25f, 1.0f, 1.0f}  // valid charts/boundaries/corners in light blue
							 }),
							 auto_flip_normals_(auto_flip_normals)
							 {

		show_ImGui_demo_window_ = false;

		show_vertices_ = false;
        show_surface_ = true;
        show_mesh_ = true;
        show_surface_borders_ = false;
        show_volume_ = false;
		surface_color_ =   vec4f(0.9f, 0.9f, 0.9f, 1.0f); // light grey. default is bleu-ish vec4f(0.5f, 0.5f, 1.0f, 1.0f)

        // init own variables
        allow_boundaries_between_opposite_labels_ = false; // parameter of StaticLabelingGraph::fill_from()

        // corners in black by default
		corners_color_[0] = 0.0f;
		corners_color_[1] = 0.0f;
		corners_color_[2] = 0.0f;
		corners_color_[3] = 1.0f;
		corners_size_ = 10.0f;

		// turning points in yellow by default
		turning_points_color_[0] = 1.0f;
		turning_points_color_[1] = 1.0f;
		turning_points_color_[2] = 0.0f;
		turning_points_color_[3] = 1.0f;
		turning_points_size_ = 10.0f;

		nb_turning_points_ = 0;

		fidelity_text_label_ = "";

		show_normals_ = false;
		normals_length_factor_ = 0.1f;

		show_feature_edges_ = true;

		state_transition(empty);
    }

	void ImGui_initialize() override {
		Application::ImGui_initialize();
		set_style("Light");
		if(GEO::FileSystem::is_file("gui.ini")) {
			// Layout modification, saved with ImGui::SaveIniSettingsToDisk()
			// Larger docked object properties panel
			ImGui::LoadIniSettingsFromDisk("gui.ini");
		}
	}

protected:

	virtual void state_transition(State new_state) {
		switch(new_state) {
			case empty:
				break;
			case triangle_mesh:
				show_mesh_ = true;
				lighting_ = true;
				show_attributes_ = false;
				break;
			case labeling:
				labeling_visu_mode_transition(VIEW_LABELING_GRAPH);
				break;
			default:
				geo_assert_not_reached;
		}
		state_ = new_state;
	}

	virtual void labeling_visu_mode_transition(int new_mode) {
		if(colormaps_.empty()) {
			// GL is not initialized yet
			// the state transition will be triggered later, in GL_initialize()
			return;
		}
		switch(new_mode) {
			case VIEW_TRIANGLE_MESH:
				show_mesh_ = true;
				lighting_ = true;
				show_attributes_ = false;
				points_groups_show_only({}); // show none
				edges_groups_show_only({}); // show none
				break;
			case VIEW_RAW_LABELING:
				show_mesh_ = true;
				lighting_ = false;
				show_attributes_ = true;
				current_colormap_index_ = COLORMAP_LABELING;
				attribute_ = fmt::format("facets.{}",LABELING_ATTRIBUTE_NAME);
				attribute_subelements_ = MESH_FACETS;
				attribute_name_ = LABELING_ATTRIBUTE_NAME;
				attribute_min_ = -0.5f;
				attribute_max_ = 5.5f;
				points_groups_show_only({}); // show none
				edges_groups_show_only({}); // show none
				break;
			case VIEW_LABELING_GRAPH:
				show_mesh_ = false;
				lighting_ = false;
				show_attributes_ = true;
				geo_assert(COLORMAP_LABELING < colormaps_.size());
				current_colormap_index_ = COLORMAP_LABELING;
				attribute_ = fmt::format("facets.{}",LABELING_ATTRIBUTE_NAME);
				attribute_subelements_ = MESH_FACETS;
				attribute_name_ = LABELING_ATTRIBUTE_NAME;
				attribute_min_ = -0.5f;
				attribute_max_ = 5.5f;
				// points in overlay
				set_points_group_color(valid_corners_group_index_,corners_color_);
				set_points_group_color(invalid_corners_group_index_,corners_color_); // no 	distinction between valid and invalid corners in this view
				points_groups_show_only({valid_corners_group_index_, invalid_corners_group_index_, turning_points_group_index_});
				// edges in overlay
				set_edges_group_color(X_boundaries_group_index_,COLORMAP_LABELING,0.084); // axis X -> color of label 0 = +X
				set_edges_group_color(Y_boundaries_group_index_,COLORMAP_LABELING,0.417); // axis Y -> color of label 2 = +Y
				set_edges_group_color(Z_boundaries_group_index_,COLORMAP_LABELING,0.750); // axis Z -> color of label 4 = +Z
				edges_groups_show_only({X_boundaries_group_index_, Y_boundaries_group_index_, Z_boundaries_group_index_});
				break;
			case VIEW_FIDELITY:
				show_mesh_ = false;
				lighting_ = false;
				show_attributes_ = true;
				current_colormap_index_ = COLORMAP_INFERNO;
				attribute_ = "facets.fidelity";
				attribute_subelements_ = MESH_FACETS;
				attribute_name_ = "fidelity";
				attribute_min_ = 0.0f; // the fidelity should not be in [0:0.5] (label too far from the normal), so setting 0.5 as the min of the colormap allows to focus the range of interest [0.5:1], but will display all values in [0:0.5] in black...
				attribute_max_ = 1.0f;
				points_groups_show_only({}); // show none
				edges_groups_show_only({}); // show none
				break;
			case VIEW_INVALID_CHARTS:
				show_mesh_ = false;
				lighting_ = false;
				show_attributes_ = true;
				current_colormap_index_ = COLORMAP_VALIDITY;
				attribute_ = "facets.on_invalid_chart";
				attribute_subelements_ = MESH_FACETS;
				attribute_name_ = "on_invalid_chart";
				attribute_min_ = 1.5;
				attribute_max_ = -0.5;
				// points in overlay
				points_groups_show_only({}); // show none
				// edges in overlay
				set_edges_group_color(X_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				set_edges_group_color(Y_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				set_edges_group_color(Z_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				edges_groups_show_only({X_boundaries_group_index_, Y_boundaries_group_index_, Z_boundaries_group_index_});

				// use mesh_gfx_.draw_surface_borders() ?

				break;
			case VIEW_INVALID_BOUNDARIES:
				show_mesh_ = false;
				lighting_ = false;
				show_attributes_ = false;
				// points in overlay
				points_groups_show_only({}); // show none
				// edges in overlay
				set_edges_group_color(X_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				set_edges_group_color(Y_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				set_edges_group_color(Z_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				edges_groups_show_only({X_boundaries_group_index_, Y_boundaries_group_index_, Z_boundaries_group_index_, invalid_boundaries_group_index_, valid_but_axisless_boundaries_group_index_});
				break;
			case VIEW_INVALID_CORNERS:
				show_mesh_ = false;
				lighting_ = false;
				show_attributes_ = false;
				// points in overlay
				set_points_group_color(valid_corners_group_index_,validity_colors_.color_as_floats(1)); // apply the color of valid LabelingGraph components
				set_points_group_color(invalid_corners_group_index_,validity_colors_.color_as_floats(0)); // apply the color of invalid LabelingGraph components
				points_groups_show_only({valid_corners_group_index_, invalid_corners_group_index_});
				// edges in overlay
				set_edges_group_color(X_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				set_edges_group_color(Y_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				set_edges_group_color(Z_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
				edges_groups_show_only({X_boundaries_group_index_, Y_boundaries_group_index_, Z_boundaries_group_index_});

				// use mesh_gfx_.draw_surface_borders() ?
				
				break;
			default:
				geo_assert_not_reached;
		}
		labeling_visu_mode_ = new_mode;
	}

	void draw_scene() override {
		SimpleMeshApplicationExt::draw_scene();
		if((state_ == triangle_mesh) && show_normals_) {
			glupSetColor4fv(GLUP_FRONT_COLOR, green);
			glupSetPointSize(10.0);
			FOR(f,mesh_.facets.nb()) { // for each 
				facet_center_ = mesh_facet_center(mesh_,f);
				normal_tip_ = facet_center_ + normals_[f] * normals_length_factor_;
                glupBegin(GLUP_LINES);
				glupPrivateVertex3dv(facet_center_.data());
				glupPrivateVertex3dv(normal_tip_.data());
                glupEnd();
				glupBegin(GLUP_POINTS);
				glupPrivateVertex3dv(normal_tip_.data());
				glupEnd();
			}
		}
		if(
			((state_ == triangle_mesh) && show_feature_edges_) ||
			((state_ == labeling) && (labeling_visu_mode_ == VIEW_TRIANGLE_MESH))
		) {
			// dark blue color (coord. 0.0 of parula colormap)
            glupEnable(GLUP_TEXTURING);
			glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
			glBindTexture(GL_TEXTURE_2D, colormaps_[COLORMAP_PARULA].texture);
			glupTextureType(GLUP_TEXTURE_2D);
			glupTextureMode(GLUP_TEXTURE_REPLACE);
			glupPrivateTexCoord1d(0.0);
			// thick lines
			glupSetMeshWidth(5);
			glupBegin(GLUP_LINES);
			for(const std::pair<index_t,index_t>& edge : feature_edges_) { // for each edge in the set of feature edges
				glupPrivateVertex3dv(mesh_.vertices.point_ptr(edge.first)); // draw first vertex
				glupPrivateVertex3dv(mesh_.vertices.point_ptr(edge.second)); // draw second vertex
			}
			glupDisable(GLUP_TEXTURING);
			glupEnd();
		}
	}

	void draw_gui() override {
		SimpleMeshApplicationExt::draw_gui();
		if(show_ImGui_demo_window_)
			ImGui::ShowDemoWindow();
	}

	// add a button on the menu bar to export the labeling graph
	void draw_menu_bar() override {
		SimpleApplication::draw_menu_bar();

		if(ImGui::BeginMainMenuBar()) {
			if(ImGui::BeginMenu("Debug")) {
				if(ImGui::MenuItem("Dump labeling graph as text file")) {
					static_labeling_graph_.dump_to_text_file("StaticLabelingGraph.txt",mesh_);
					fmt::println(Logger::out("I/O"),"Exported to StaticLabelingGraph.txt"); Logger::out("I/O").flush();
				}
				if(ImGui::MenuItem("Dump labeling graph as D3 graph")) {
					static_labeling_graph_.dump_to_D3_graph("StaticLabelingGraph.json");
					fmt::println(Logger::out("I/O"),"Exported to StaticLabelingGraph.json"); Logger::out("I/O").flush();
				}
				if(ImGui::MenuItem("Dump boundaries as mesh")) {
					CustomMeshHalfedges mesh_he(mesh_);
					dump_all_boundaries_with_indices_and_axes("boundaries",mesh_,static_labeling_graph_);
				}
				if (ImGui::MenuItem("Show ImGui demo window", NULL, show_ImGui_demo_window_)) {
					show_ImGui_demo_window_ = !show_ImGui_demo_window_;
				}
				ImGui::EndMenu();
			}
			ImGui::EndMainMenuBar();
		}
	}

    void draw_object_properties() override {
        switch (state_) {

		case triangle_mesh:

			ImGui::Checkbox("Show normals",&show_normals_);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Draw facet normals as segments. The point is the tip of the vector.");
			ImGui::SliderFloat("Normals length factor",&normals_length_factor_,0.0f,50.0f);
			ImGui::Checkbox("Show feature edges",&show_feature_edges_);
			ImGui::SameLine();
			ImGui::Text("(%ld)",feature_edges_.size());
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Draw feature edges in blue on top of the mesh");
			if(ImGui::Button("Rotate mesh according to principal axes")) {
				rotate_mesh_according_to_principal_axes(mesh_);
				mesh_gfx_.set_mesh(&mesh_); // re-link the MeshGfx to the mesh
				// update normals
				FOR(f,mesh_.facets.nb()) {
					normals_[f] = normalize(Geom::mesh_facet_normal(mesh_,f));
				}
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Compute principal axes of the point cloud and rotate the mesh to be aligned with them");

			ImGui::Separator();

			ImGui::Checkbox("Allow boundaries between opposite labels",&allow_boundaries_between_opposite_labels_);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("If on, boundaries between opposite labels (e.g. +X and -X)\ncan be considered valid if they only contain > 180° angles");

			if(ImGui::Button("Compute naive labeling")) {
				naive_labeling(mesh_,normals_,LABELING_ATTRIBUTE_NAME);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				state_transition(labeling);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Associate each facet to the label the closest to its normal");

			if(ImGui::Button("Compute tweaked naive labeling")) {
				tweaked_naive_labeling(mesh_,normals_,LABELING_ATTRIBUTE_NAME);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				state_transition(labeling);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Like the naive labeling, but facets normals close to multiple labels\nare slightly rotated before choosing the closest label,\nto avoid labeling fragmentation on subsurfaces\npoorly aligned with XY, XZ, YZ planes");
			
			if(ImGui::Button("Smart init labeling")) {
				smart_init_labeling(mesh_,normals_,LABELING_ATTRIBUTE_NAME,feature_edges_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				state_transition(labeling);
			}

			break;
		case labeling:

			if(ImGui::Button("Remove labeling")) {
				ImGui::OpenPopup("Remove labeling ?");
			}
			ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f)); // Always center the modal when appearing
			if(ImGui::BeginPopupModal("Remove labeling ?", NULL, ImGuiWindowFlags_AlwaysAutoResize))
			{
				ImGui::TextUnformatted("If not manually saved before, the labeling will be lost.");
				ImGui::TextUnformatted("Are you sure you want to remove the labeling?");
				ImGui::PushStyleColor(ImGuiCol_Button, 			ImVec4(0.9f, 0.4f, 0.4f, 1.0f)); // red
				ImGui::PushStyleColor(ImGuiCol_ButtonHovered,	ImVec4(0.9f, 0.2f, 0.2f, 1.0f)); // darker red
				ImGui::PushStyleColor(ImGuiCol_ButtonActive,	ImVec4(0.9f, 0.0f, 0.0f, 1.0f));
				if(ImGui::Button("Yes, remove the labeling", ImVec2(200, 0))) {
					clear_scene_overlay();
					state_transition(triangle_mesh);
					ImGui::CloseCurrentPopup();
				}
				ImGui::PopStyleColor(3);
				ImGui::SetItemDefaultFocus();
				ImGui::SameLine();
				if(ImGui::Button("No, cancel", ImVec2(120, 0))) {
					ImGui::CloseCurrentPopup();
				}
				ImGui::EndPopup();
			}

			ImGui::Checkbox("Allow boundaries between opposite labels",&allow_boundaries_between_opposite_labels_);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("If on, boundaries between opposite labels (e.g. +X and -X)\ncan be considered valid if they only contain > 180° angles");

			ImGui::BeginDisabled( allow_boundaries_between_opposite_labels_ == static_labeling_graph_.is_allowing_boundaries_between_opposite_labels() ); // allow to recompute only if the UI control value changed
			if(ImGui::Button("Recompute labeling graph")) {
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Update charts, boundaries and corners according to the new value of \"Allow boundaries between opposite labels\"");
			ImGui::EndDisabled();

			ImGui::Separator();

			if(ImGui::RadioButton("View triangle mesh",&labeling_visu_mode_,VIEW_TRIANGLE_MESH))
				labeling_visu_mode_transition(VIEW_TRIANGLE_MESH);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Show the triangle mesh without the labeling");

			if(ImGui::RadioButton("View raw labeling",&labeling_visu_mode_,VIEW_RAW_LABELING))
				labeling_visu_mode_transition(VIEW_RAW_LABELING);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Show the triangle mesh and the labeling");

			ImGui::BeginDisabled( (labeling_visu_mode_!=VIEW_RAW_LABELING) && (labeling_visu_mode_!=VIEW_LABELING_GRAPH) );
			if(ImGui::ColorEdit4WithPalette("Label 0 = +X", labeling_colors_.color_as_floats(0))) {
				labeling_colors_.update_chars_of_color(0);
				update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
			}
			if(ImGui::ColorEdit4WithPalette("Label 1 = -X", labeling_colors_.color_as_floats(1))) {
				labeling_colors_.update_chars_of_color(1);
				update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
			}
			if(ImGui::ColorEdit4WithPalette("Label 2 = +Y", labeling_colors_.color_as_floats(2))) {
				labeling_colors_.update_chars_of_color(2);
				update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
			}
			if(ImGui::ColorEdit4WithPalette("Label 3 = -Y", labeling_colors_.color_as_floats(3))) {
				labeling_colors_.update_chars_of_color(3);
				update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
			}
			if(ImGui::ColorEdit4WithPalette("Label 4 = +Z", labeling_colors_.color_as_floats(4))) {
				labeling_colors_.update_chars_of_color(4);
				update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
			}
			if(ImGui::ColorEdit4WithPalette("Label 5 = -Z", labeling_colors_.color_as_floats(5))) {
				labeling_colors_.update_chars_of_color(5);
				update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
			}
			ImGui::EndDisabled();

			if(ImGui::RadioButton("View labeling graph",&labeling_visu_mode_,VIEW_LABELING_GRAPH))
				labeling_visu_mode_transition(VIEW_LABELING_GRAPH);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Show the charts, the boundaries, the corners and the turning-points computed from the labeling");

			ImGui::Text("%ld charts, %ld boundaries",static_labeling_graph_.nb_charts(),static_labeling_graph_.nb_boundaries());

			ImGui::BeginDisabled(labeling_visu_mode_!=VIEW_LABELING_GRAPH);

			ImGui::ColorEdit4WithPalette("Corners", corners_color_);
			ImGui::SameLine();
			ImGui::Text("(%ld)",static_labeling_graph_.nb_corners());
			ImGui::SliderFloat("Corners size", &corners_size_, 0.0f, 50.0f, "%.1f");

			ImGui::ColorEdit4WithPalette("Turning points", turning_points_color_);
			ImGui::SameLine();
			ImGui::Text("(%ld)",nb_turning_points_);
			ImGui::SliderFloat("Turning points size", &turning_points_size_, 0.0f, 50.0f, "%.1f");

			ImGui::EndDisabled();

			if(ImGui::RadioButton("View fidelity",&labeling_visu_mode_,VIEW_FIDELITY))
				labeling_visu_mode_transition(VIEW_FIDELITY);
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Color the facets according to the angle between the normal and the assigned direction.\nYellow = small angle, black = wide angle.\nClick on a facet to print its fidelity.");
			
			ImGui::TextUnformatted(fidelity_text_label_.c_str());

			if(ImGui::RadioButton("View invalid charts",&labeling_visu_mode_,VIEW_INVALID_CHARTS))
				labeling_visu_mode_transition(VIEW_INVALID_CHARTS);
			ImGui::SameLine();
			ImGui::Text("(%ld)",static_labeling_graph_.nb_invalid_charts());

			if(ImGui::RadioButton("View invalid boundaries",&labeling_visu_mode_,VIEW_INVALID_BOUNDARIES))
				labeling_visu_mode_transition(VIEW_INVALID_BOUNDARIES);
			ImGui::SameLine();
			ImGui::Text("(%ld)",static_labeling_graph_.nb_invalid_boundaries());

			if(ImGui::RadioButton("View invalid corners",&labeling_visu_mode_,VIEW_INVALID_CORNERS))
				labeling_visu_mode_transition(VIEW_INVALID_CORNERS);
			ImGui::SameLine();
			ImGui::Text("(%ld)",static_labeling_graph_.nb_invalid_corners());
			
			ImGui::BeginDisabled( (labeling_visu_mode_!=VIEW_INVALID_CHARTS) && 
								  (labeling_visu_mode_!=VIEW_INVALID_BOUNDARIES) && 
								  (labeling_visu_mode_!=VIEW_INVALID_CORNERS) );
			if(ImGui::ColorEdit4WithPalette("Invalid", validity_colors_.color_as_floats(0))) {
				validity_colors_.update_chars_of_color(0);
				update_GL_texture(COLORMAP_VALIDITY,2,1,validity_colors_.as_chars());
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Color of invalid charts/boundaries/corners");
			if(ImGui::ColorEdit4WithPalette("Valid", validity_colors_.color_as_floats(1))) {
				validity_colors_.update_chars_of_color(1);
				update_GL_texture(COLORMAP_VALIDITY,2,1,validity_colors_.as_chars());
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Color of valid charts/boundaries/corners");
			ImGui::EndDisabled();

			if(static_labeling_graph_.is_valid()) {
				ImGui::TextColored(ImVec4(0.0f,0.5f,0.0f,1.0f),"Valid labeling");
				ImGui::SameLine();
				ImGui::TextDisabled("(?)");
				ImGui::SetItemTooltip("The current labeling is a valid polycube representation");
			}
			else {
				ImGui::TextColored(ImVec4(0.8f,0.0f,0.0f,1.0f),"Invalid labeling");
				ImGui::SameLine();
				ImGui::TextDisabled("(?)");
				ImGui::SetItemTooltip("The current labeling is not valid polycube representation.\nValence or adjacency of some components (charts, boundaries, corners) cannot turn into polycube components.");
			}

			if(nb_turning_points_==0) {
				ImGui::TextColored(ImVec4(0.0f,0.5f,0.0f,1.0f),"All monotone boundaries");
				ImGui::SameLine();
				ImGui::TextDisabled("(?)");
				ImGui::SetItemTooltip("There are no turning-points");
			}
			else {
				ImGui::TextColored(ImVec4(0.8f,0.0f,0.0f,1.0f),"Non-monotone boundaries");
				ImGui::SameLine();
				ImGui::TextDisabled("(?)");
				ImGui::SetItemTooltip("Some boundaries contain turning-points");
			}
			
			break;
		
		default:
			break;
		}
    }

    bool load(const std::string& filename) override {

        //// based on SimpleMeshApplication::load() ////////////////////////////////////////////////////////////////

        if(!FileSystem::is_file(filename)) {
            Logger::out("I/O") << "is not a file" << std::endl;
        }

        // added : load a labeling
        if(FileSystem::extension(filename)=="txt") {

			if(state_ == empty) {
				fmt::println(Logger::err("I/O"),"You need to import a triangle mesh before importing a labeling"); Logger::err("I/O").flush();
				return false;
			}

			if(!load_labeling(filename,mesh_,LABELING_ATTRIBUTE_NAME)) {
				fmt::println("load_labeling() not ok"); fflush(stdout);
				// Should the labeling be removed ?
				// If a labeling was already displayed, it should be restored...
				mesh_.facets.clear(false,true);
				clear_scene_overlay();
				state_transition(triangle_mesh);
				return false;
			}

			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			state_transition(labeling);
			return true;
		}

        mesh_gfx_.set_mesh(nullptr);
		feature_edges_.clear();
		mesh_.clear(false,true);
        MeshIOFlags flags;
        if(!mesh_load(filename, mesh_, flags)) {
			state_transition(empty); // added
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

		if(auto_flip_normals_) {
			// ensure facet normals are outward
			if(facet_normals_are_inward(mesh_)) {
				flip_facet_normals(mesh_);
				fmt::println(Logger::warn("normals dir."),"Facet normals of the input mesh were inward");
				fmt::println(Logger::warn("normals dir."),"You should flip them and update the file for consistency.");
				Logger::warn("normals dir.").flush();
			}
		}

		// compute facet normals
		// TODO use geogram's compute_normals() ? in mesh/mesh_geometry.h
		normals_.resize(mesh_.facets.nb());
        FOR(f,mesh_.facets.nb()) {
            normals_[f] = normalize(Geom::mesh_facet_normal(mesh_,f));
        }

		if(mesh_.edges.nb() > 0) {
			// the loaded mesh contains edges -> expect they are feature edges
			adj_facets_.clear();
			remove_feature_edges_with_low_dihedral_angle(mesh_,adj_facets_);

			// transfert feature edges from mesh_.edges to the feature_edges_ set
			transfer_feature_edges(mesh_,feature_edges_);
		}

		clear_scene_overlay();
		state_transition(triangle_mesh);

        return true;
    }

	std::string supported_write_file_extensions() override {
        return SimpleMeshApplication::supported_write_file_extensions() + ";txt"; // add .txt in supported write file extensions
    }

	bool save(const std::string& filename) override {
		if(String::string_ends_with(filename,".txt")) { // bypass inherited save behavior in case of a .txt file -> save the labeling only
			save_labeling(filename,mesh_,LABELING_ATTRIBUTE_NAME);
			fmt::println(Logger::out("I/O"),"Labeling saved to {}",filename); Logger::out("I/O").flush();
			return true;
		}
		else {
			return SimpleMeshApplication::save(filename);
		}
    }

	void GL_initialize() override {
        SimpleMeshApplicationExt::GL_initialize();
		init_rgba_colormap("labeling",6,1,labeling_colors_.as_chars());
		init_rgba_colormap("validity",2,1,validity_colors_.as_chars());
		if((mesh_.vertices.nb()!=0) && mesh_.facets.attributes().is_defined(LABELING_ATTRIBUTE_NAME)) {
			state_transition(labeling);
		}
    }

	void mouse_button_callback(int button, int action, int mods, int source) override {
		if((action==EVENT_ACTION_DOWN) && (button == 0) ) { // if left click
			if( (state_ == labeling) && ((labeling_visu_mode_ == VIEW_RAW_LABELING) ||
												 (labeling_visu_mode_ == VIEW_FIDELITY)) ) { // if current mode is "view raw labeling" or "view fidelity"
				index_t facet_index = pick(MESH_FACETS);
				if ( (facet_index != index_t(-1)) && (facet_index < mesh_.facets.nb()) ) {
					Attribute<index_t> label(mesh_.facets.attributes(), LABELING_ATTRIBUTE_NAME);
					Attribute<double> per_facet_fidelity(mesh_.facets.attributes(), "fidelity");
					vec3 label_direction = label2vector[label[facet_index]];
					fmt::println(Logger::out("fidelity"),"facet #{} : normal=({:.4f},{:.4f},{:.4f}), label={}=({:.4f},{:.4f},{:.4f}) -> fidelity={:.4f}",
						facet_index,
						normals_[facet_index].x,
						normals_[facet_index].y,
						normals_[facet_index].z,
						LABEL2STR(label[facet_index]),
						label_direction.x,
						label_direction.y,
						label_direction.z,
						per_facet_fidelity[facet_index]
					);
					Logger::out("fidelity").flush();
				}
			}
		}
		SimpleMeshApplication::mouse_button_callback(button,action,mods,source);
	}

    virtual void update_static_labeling_graph(bool allow_boundaries_between_opposite_labels) {

		// compute charts, boundaries and corners of the labeling
		static_labeling_graph_.fill_from(mesh_,LABELING_ATTRIBUTE_NAME,feature_edges_,allow_boundaries_between_opposite_labels);
		nb_turning_points_ = static_labeling_graph_.nb_turning_points();

		clear_scene_overlay();

		valid_corners_group_index_ = new_points_group(corners_color_,&corners_size_,true);
		invalid_corners_group_index_ = new_points_group(corners_color_,&corners_size_,true);
		turning_points_group_index_ = new_points_group(turning_points_color_,&turning_points_size_,true);
		X_boundaries_group_index_ = new_edges_group(COLORMAP_LABELING,0.084,6,false); // axis X -> color of label 0 = +X
		Y_boundaries_group_index_ = new_edges_group(COLORMAP_LABELING,0.417,6,false); // axis Y -> color of label 2 = +Y
		Z_boundaries_group_index_ = new_edges_group(COLORMAP_LABELING,0.750,6,false); // axis Z -> color of label 4 = +Z
		invalid_boundaries_group_index_ = new_edges_group(COLORMAP_VALIDITY,0.0,6,false); // color of invalid LabelingGraph components
		valid_but_axisless_boundaries_group_index_ = new_edges_group(COLORMAP_VALIDITY,1.0,6,false); // color of valid LabelingGraph components

		for(std::size_t i = 0; i < static_labeling_graph_.nb_corners(); ++i) {
			const double* coordinates = mesh_.vertices.point_ptr(
				static_labeling_graph_.corners[i].vertex
			);
			add_point_to_group(static_labeling_graph_.corners[i].is_valid ? valid_corners_group_index_ : invalid_corners_group_index_,coordinates[0], coordinates[1], coordinates[2]);
		}

		for(index_t i : static_labeling_graph_.non_monotone_boundaries) {
			FOR(j,static_labeling_graph_.boundaries[i].turning_points.size()) {
				vec3 coordinates = mesh_vertex(mesh_,static_labeling_graph_.boundaries[i].turning_point_vertex(j,mesh_));
				add_point_to_group(turning_points_group_index_,coordinates.x,coordinates.y,coordinates.z);
			}
		}

		std::size_t group_index;
		for(std::size_t i = 0; i < static_labeling_graph_.nb_boundaries(); ++i) {
			const Boundary& boundary = static_labeling_graph_.boundaries[i];
			for(const auto& be : boundary.halfedges) { // for each boundary edge of this boundary
				const double* coordinates_first_point = halfedge_vertex_from(mesh_,be).data();
				const double* coordinates_second_point = halfedge_vertex_to(mesh_,be).data();
				switch(boundary.axis) {
					case -1: // may or may not be invalid
						if(boundary.is_valid) {
							group_index = valid_but_axisless_boundaries_group_index_;
						}
						else {
							group_index = invalid_boundaries_group_index_;
						}
						break;
					case 0: // always valid
						group_index = X_boundaries_group_index_;
						break;
					case 1: // always valid
						group_index = Y_boundaries_group_index_;
						break;
					case 2: // always valid
						group_index = Z_boundaries_group_index_;
						break;
					default:
						geo_assert_not_reached;
				}
				add_edge_to_group(
					group_index,
					coordinates_first_point[0], coordinates_first_point[1], coordinates_first_point[2], // first point
					coordinates_second_point[0], coordinates_second_point[1], coordinates_second_point[2] // second point
				);
			}
		}

		set_points_group_visibility(valid_corners_group_index_,true);
		set_points_group_visibility(invalid_corners_group_index_,true);
		set_edges_group_visibility(X_boundaries_group_index_,true);
		set_edges_group_visibility(Y_boundaries_group_index_,true);
		set_edges_group_visibility(Z_boundaries_group_index_,true);

		BasicStats stats;
		compute_per_facet_fidelity(mesh_,normals_,LABELING_ATTRIBUTE_NAME,"fidelity",stats);
		fidelity_text_label_ = fmt::format("min={:.4f} | max={:.4f} | avg={:.4f}",stats.min(),stats.max(),stats.avg());
	}

protected:

	bool show_ImGui_demo_window_;
    bool allow_boundaries_between_opposite_labels_;
	ColorArray labeling_colors_;
	float corners_color_[4];
	float turning_points_color_[4];
	std::size_t nb_turning_points_;
	ColorArray validity_colors_;
	StaticLabelingGraph static_labeling_graph_;
	State state_; // should only be modified by state_transition()
	int labeling_visu_mode_; // not a enum, to be used in ImGui. should only be modified by labeling_visu_mode_transition() and ImGui::RadioButton
	std::size_t valid_corners_group_index_;
	std::size_t invalid_corners_group_index_;
	float corners_size_;
	std::size_t turning_points_group_index_;
	float turning_points_size_;
	std::size_t X_boundaries_group_index_;
	std::size_t Y_boundaries_group_index_;
	std::size_t Z_boundaries_group_index_;
	std::size_t invalid_boundaries_group_index_;
	std::size_t valid_but_axisless_boundaries_group_index_;
	std::vector<vec3> normals_; // facet normals
	std::string fidelity_text_label_;
	bool show_normals_; // optional visu overlay when state_ is triangle_mesh
	vec3 facet_center_;
	vec3 normal_tip_;
	float normals_length_factor_;
	bool auto_flip_normals_;
	bool show_feature_edges_;
	std::vector<std::vector<index_t>> adj_facets_; // for each vertex, store adjacent facets. no ordering
	std::set<std::pair<index_t,index_t>> feature_edges_;
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