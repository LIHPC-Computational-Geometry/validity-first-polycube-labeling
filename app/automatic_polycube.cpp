#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/file_system.h> // for is_file(), extension() in load()
#include <geogram/basic/command_line.h> // for get_arg_bool() in load()
#include <geogram/mesh/mesh_io.h> // for MeshIOFlags, mesh_load() in load()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "SimpleMeshApplicationExt.h"
#include "LabelingGraph.h"   // for StaticLabelingGraph
#include "labeling.h"		 // for load_labeling(), naive_labeling()

#define LABELING_ATTRIBUTE_NAME "label"

// values for *_labeling_visu_mode_
#define VIEW_TRIANGLE_MESH		0
#define VIEW_RAW_LABELING		1
#define VIEW_LABELING_GRAPH		2
#define VIEW_INVALID_CHARTS		3
#define VIEW_INVALID_BOUNDARIES	4
#define VIEW_INVALID_CORNERS	5

// new colormaps
#define COLORMAP_LABELING       (SIMPLE_APPLICATION_NB_COLORMAPS)
#define COLORMAP_VALIDITY       ((SIMPLE_APPLICATION_NB_COLORMAPS)+1)

using namespace GEO;

class AutomaticPolycubeApp : private SimpleMeshApplicationExt {
public:

    enum State {
        empty,
        triangle_mesh,
        labeling
	};

    AutomaticPolycubeApp() : SimpleMeshApplicationExt("automatic_polycube") {

		show_vertices_ = false;
        show_surface_ = true;
        show_mesh_ = true;
        show_surface_borders_ = false;
        show_volume_ = false;
		surface_color_ =   vec4f(0.9f, 0.9f, 0.9f, 1.0f); // light grey. default is bleu-ish vec4f(0.5f, 0.5f, 1.0f, 1.0f)

        // init own variables
        allow_boundaries_between_opposite_labels_ = false; // parameter of StaticLabelingGraph::fill_from()

        // default colors
		// label 0 -> red
		labeling_colors_as_char_[4*0+0] = 255;
		labeling_colors_as_char_[4*0+1] = 0;
		labeling_colors_as_char_[4*0+2] = 0;
		labeling_colors_as_char_[4*0+3] = 255;
		// label 1 -> dark red
		labeling_colors_as_char_[4*1+0] = 152;
		labeling_colors_as_char_[4*1+1] = 0;
		labeling_colors_as_char_[4*1+2] = 0;
		labeling_colors_as_char_[4*1+3] = 255;
		// label 2 -> green
		labeling_colors_as_char_[4*2+0] = 0;
		labeling_colors_as_char_[4*2+1] = 255;
		labeling_colors_as_char_[4*2+2] = 0;
		labeling_colors_as_char_[4*2+3] = 255;
		// label 3 -> dark green
		labeling_colors_as_char_[4*3+0] = 0;
		labeling_colors_as_char_[4*3+1] = 152;
		labeling_colors_as_char_[4*3+2] = 0;
		labeling_colors_as_char_[4*3+3] = 255;
		// label 4 -> blue
		labeling_colors_as_char_[4*4+0] = 0;
		labeling_colors_as_char_[4*4+1] = 0;
		labeling_colors_as_char_[4*4+2] = 255;
		labeling_colors_as_char_[4*4+3] = 255;
		// label 5 -> dark blue
		labeling_colors_as_char_[4*5+0] = 0;
		labeling_colors_as_char_[4*5+1] = 0;
		labeling_colors_as_char_[4*5+2] = 152;
		labeling_colors_as_char_[4*5+3] = 255;

        FOR(i,4*6) {
			labeling_colors_as_float_[i] = ((float) labeling_colors_as_char_[i]) / 255.0f;
		}

        // corners in black by default
		corners_color_[0] = 0.0f;
		corners_color_[1] = 0.0f;
		corners_color_[2] = 0.0f;
		corners_color_[3] = 1.0f;

		// turning points in yellow by default
		turning_points_color_[0] = 1.0f;
		turning_points_color_[1] = 1.0f;
		turning_points_color_[2] = 0.0f;
		turning_points_color_[3] = 1.0f;

		nb_turning_points = 0;

		// invalid charts/boundaries/corners in red by default
		validity_colors_as_char_[4*0+0] = 255;
		validity_colors_as_char_[4*0+1] = 64;
		validity_colors_as_char_[4*0+2] = 64;
		validity_colors_as_char_[4*0+3] = 255;

		// valid charts/boundaries/corners in blue by default
		validity_colors_as_char_[4*1+0] = 64;
		validity_colors_as_char_[4*1+1] = 64;
		validity_colors_as_char_[4*1+2] = 255;
		validity_colors_as_char_[4*1+3] = 255;

        FOR(i,2*4) {
			validity_colors_as_float_[i] = ((float) validity_colors_as_char_[i]) / 255.0f;
		}

		previous_labeling_visu_mode_ = VIEW_RAW_LABELING;
		current_labeling_visu_mode_ = VIEW_RAW_LABELING;
    }

    void start(int argc, char** argv) override {
        SimpleMeshApplicationExt::start(argc,argv);
    }
private:

    void draw_scene() override {
        // manage state transitions
		if(previous_state_ != current_state_) {
			switch(current_state_) {
				case empty:
					break;
				case triangle_mesh:
					show_mesh_ = true;
					lighting_ = true;
					show_attributes_ = false;
					break;
				case labeling:
					previous_labeling_visu_mode_ = VIEW_TRIANGLE_MESH;
					current_labeling_visu_mode_ = VIEW_LABELING_GRAPH;
					// will trigger a GUI settings update because !=
					break;
				default:
					geo_assert_not_reached;
			}
			previous_state_ = current_state_;
		}

		// manage labeling visu mode transitions (mode editable when current_state_==labeling)
		if(previous_labeling_visu_mode_ != current_labeling_visu_mode_) {
			switch(current_labeling_visu_mode_) {
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
					current_colormap_texture_ = colormaps_[COLORMAP_LABELING].texture;
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
					current_colormap_texture_ = colormaps_[COLORMAP_LABELING].texture;
					attribute_ = fmt::format("facets.{}",LABELING_ATTRIBUTE_NAME);
					attribute_subelements_ = MESH_FACETS;
					attribute_name_ = LABELING_ATTRIBUTE_NAME;
					attribute_min_ = -0.5f;
					attribute_max_ = 5.5f;
					// points in overlay
					set_points_group_color(valid_corners_group_index,corners_color_);
					set_points_group_color(invalid_corners_group_index,corners_color_); // no 	distinction between valid and invalid corners in this view
					points_groups_show_only({valid_corners_group_index, invalid_corners_group_index, turning_points_group_index});
					// edges in overlay
					set_edges_group_color(X_boundaries_group_index,&labeling_colors_as_float_[4*0]); // axis X -> color of label 0 = +X
					set_edges_group_color(Y_boundaries_group_index,&labeling_colors_as_float_[4*2]); // axis Y -> color of label 2 = +Y
					set_edges_group_color(Z_boundaries_group_index,&labeling_colors_as_float_[4*4]); // axis Z -> color of label 4 = +Z
					edges_groups_show_only({X_boundaries_group_index, Y_boundaries_group_index, Z_boundaries_group_index});
					break;
				case VIEW_INVALID_CHARTS:
					show_mesh_ = false;
					lighting_ = false;
					show_attributes_ = true;
					current_colormap_texture_ = colormaps_[COLORMAP_VALIDITY].texture;
					attribute_ = "facets.on_invalid_chart";
					attribute_subelements_ = MESH_FACETS;
					attribute_name_ = "on_invalid_chart";
					attribute_min_ = 1.5;
					attribute_max_ = -0.5;
					// points in overlay
					points_groups_show_only({}); // show none
					// edges in overlay
					set_edges_group_color(X_boundaries_group_index,&labeling_colors_as_float_[4*0]); // axis X -> color of label 0 = +X
					set_edges_group_color(Y_boundaries_group_index,&labeling_colors_as_float_[4*2]); // axis Y -> color of label 2 = +Y
					set_edges_group_color(Z_boundaries_group_index,&labeling_colors_as_float_[4*4]); // axis Z -> color of label 4 = +Z
					edges_groups_show_only({X_boundaries_group_index, Y_boundaries_group_index, Z_boundaries_group_index});

					// use mesh_gfx_.draw_surface_borders() ?

					break;
				case VIEW_INVALID_BOUNDARIES:
					show_mesh_ = false;
					lighting_ = false;
					show_attributes_ = false;
					// points in overlay
					points_groups_show_only({}); // show none
					// edges in overlay
					set_edges_group_color(X_boundaries_group_index,&validity_colors_as_float_[4*1]); // apply the color of valid LabelingGraph components
					set_edges_group_color(Y_boundaries_group_index,&validity_colors_as_float_[4*1]); // apply the color of valid LabelingGraph components
					set_edges_group_color(Z_boundaries_group_index,&validity_colors_as_float_[4*1]); // apply the color of valid LabelingGraph components
					edges_groups_show_only({X_boundaries_group_index, Y_boundaries_group_index, Z_boundaries_group_index, invalid_boundaries_group_index, valid_but_axisless_boundaries_group_index});
					break;
				case VIEW_INVALID_CORNERS:
					show_mesh_ = false;
					lighting_ = false;
					show_attributes_ = false;
					// points in overlay
					set_points_group_color(valid_corners_group_index,&validity_colors_as_float_[4*1]); // apply the color of valid LabelingGraph components
					set_points_group_color(invalid_corners_group_index,&validity_colors_as_float_[4*0]); // apply the color of invalid LabelingGraph components
					points_groups_show_only({valid_corners_group_index, invalid_corners_group_index});
					// edges in overlay
					set_edges_group_color(X_boundaries_group_index,&labeling_colors_as_float_[4*0]); // axis X -> color of label 0 = +X
					set_edges_group_color(Y_boundaries_group_index,&labeling_colors_as_float_[4*2]); // axis Y -> color of label 2 = +Y
					set_edges_group_color(Z_boundaries_group_index,&labeling_colors_as_float_[4*4]); // axis Z -> color of label 4 = +Z
					edges_groups_show_only({X_boundaries_group_index, Y_boundaries_group_index, Z_boundaries_group_index});

					// use mesh_gfx_.draw_surface_borders() ?
					
					break;
				default:
					geo_assert_not_reached;
			}
			previous_labeling_visu_mode_ = current_labeling_visu_mode_;
		}
        SimpleMeshApplicationExt::draw_scene();
    }

    void draw_object_properties() override {
        switch (current_state_) {

		case triangle_mesh:

			ImGui::Checkbox("Allow boundaries between opposite labels",&allow_boundaries_between_opposite_labels_);

			if(ImGui::Button("Compute naive labeling")) {

				naive_labeling(mesh_,LABELING_ATTRIBUTE_NAME);

				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);

				current_state_ = labeling;
			}
			
			break;
		case labeling:

			ImGui::Checkbox("Allow boundaries between opposite labels",&allow_boundaries_between_opposite_labels_);

			ImGui::BeginDisabled( allow_boundaries_between_opposite_labels_ == static_labeling_graph.is_allowing_boundaries_between_opposite_labels() ); // allow to recompute only if the UI control value changed
			if(ImGui::Button("Recompute labeling graph")) {
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::EndDisabled();

			ImGui::Separator();

			ImGui::RadioButton("View triangle mesh",&current_labeling_visu_mode_,VIEW_TRIANGLE_MESH);

			ImGui::RadioButton("View raw labeling",&current_labeling_visu_mode_,VIEW_RAW_LABELING);

			ImGui::BeginDisabled( (current_labeling_visu_mode_!=VIEW_RAW_LABELING) && (current_labeling_visu_mode_!=VIEW_LABELING_GRAPH) );
			if(ImGui::ColorEdit4WithPalette("Label 0 = +X", &labeling_colors_as_float_[4*0])) {
				FOR(i,4)
					labeling_colors_as_char_[4*0+i] = static_cast<unsigned char>(labeling_colors_as_float_[4*0+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_LABELING].texture,6,1,labeling_colors_as_char_);
			}
			if(ImGui::ColorEdit4WithPalette("Label 1 = -X", &labeling_colors_as_float_[4*1])) {
				FOR(i,4)
					labeling_colors_as_char_[4*1+i] = static_cast<unsigned char>(labeling_colors_as_float_[4*1+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_LABELING].texture,6,1,labeling_colors_as_char_);
			}
			if(ImGui::ColorEdit4WithPalette("Label 2 = +Y", &labeling_colors_as_float_[4*2])) {
				FOR(i,4)
					labeling_colors_as_char_[4*2+i] = static_cast<unsigned char>(labeling_colors_as_float_[4*2+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_LABELING].texture,6,1,labeling_colors_as_char_);
			}
			if(ImGui::ColorEdit4WithPalette("Label 3 = -Y", &labeling_colors_as_float_[4*3])) {
				FOR(i,4)
					labeling_colors_as_char_[4*3+i] = static_cast<unsigned char>(labeling_colors_as_float_[4*3+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_LABELING].texture,6,1,labeling_colors_as_char_);
			}
			if(ImGui::ColorEdit4WithPalette("Label 4 = +Z", &labeling_colors_as_float_[4*4])) {
				FOR(i,4)
					labeling_colors_as_char_[4*4+i] = static_cast<unsigned char>(labeling_colors_as_float_[4*4+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_LABELING].texture,6,1,labeling_colors_as_char_);
			}
			if(ImGui::ColorEdit4WithPalette("Label 5 = -Z", &labeling_colors_as_float_[4*5])) {
				FOR(i,4)
					labeling_colors_as_char_[4*5+i] = static_cast<unsigned char>(labeling_colors_as_float_[4*5+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_LABELING].texture,6,1,labeling_colors_as_char_);
			}
			ImGui::EndDisabled();

			ImGui::RadioButton("View labeling graph",&current_labeling_visu_mode_,VIEW_LABELING_GRAPH);

			ImGui::BeginDisabled(current_labeling_visu_mode_!=VIEW_LABELING_GRAPH);
			ImGui::ColorEdit4WithPalette("Corners", corners_color_);
			ImGui::ColorEdit4WithPalette(fmt::format("Turning points ({})",nb_turning_points).c_str(), turning_points_color_);
			ImGui::EndDisabled();

			ImGui::RadioButton(fmt::format("View invalid charts ({})",static_labeling_graph.nb_invalid_charts()).c_str(),&current_labeling_visu_mode_,VIEW_INVALID_CHARTS);

			ImGui::RadioButton(fmt::format("View invalid boundaries ({})",static_labeling_graph.nb_invalid_boundaries()).c_str(),&current_labeling_visu_mode_,VIEW_INVALID_BOUNDARIES);

			ImGui::RadioButton(fmt::format("View invalid corners ({})",static_labeling_graph.nb_invalid_corners()).c_str(),&current_labeling_visu_mode_,VIEW_INVALID_CORNERS);
			
			ImGui::BeginDisabled( (current_labeling_visu_mode_!=VIEW_INVALID_CHARTS) && 
								  (current_labeling_visu_mode_!=VIEW_INVALID_BOUNDARIES) && 
								  (current_labeling_visu_mode_!=VIEW_INVALID_CORNERS) );
			if(ImGui::ColorEdit4WithPalette("Invalid", &validity_colors_as_float_[4*0])) {
				FOR(i,4)
					validity_colors_as_char_[4*0+i] = static_cast<unsigned char>(validity_colors_as_float_[4*0+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_VALIDITY].texture,2,1,validity_colors_as_char_);
			}
			if(ImGui::ColorEdit4WithPalette("Valid", &validity_colors_as_float_[4*1])) {
				FOR(i,4)
					validity_colors_as_char_[4*1+i] = static_cast<unsigned char>(validity_colors_as_float_[4*1+i]*255.0f);
				update_GL_texture(colormaps_[COLORMAP_VALIDITY].texture,2,1,validity_colors_as_char_);
			}
			ImGui::EndDisabled();

			if(static_labeling_graph.is_valid()) {
				ImGui::TextColored(ImVec4(0.0f,0.5f,0.0f,1.0f),"Valid labeling");
			}
			else {
				ImGui::TextColored(ImVec4(0.8f,0.0f,0.0f,1.0f),"Invalid labeling");
			}

			if(nb_turning_points==0) {
				ImGui::TextColored(ImVec4(0.0f,0.5f,0.0f,1.0f),"All monotone boundaries");
			}
			else {
				ImGui::TextColored(ImVec4(0.8f,0.0f,0.0f,1.0f),"Non-monotone boundaries");
			}

			ImGui::Separator();

			ImGui::Text("Fix labeling");

			if(ImGui::Button("Remove surrounded charts")) {
				unsigned int nb_chart_modified = remove_surrounded_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				fmt::println(Logger::out("fix_labeling"),"{} chart(s) modified with the surrounding label",nb_chart_modified); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid boundaries")) {
				unsigned int nb_chart_modified = fix_invalid_boundaries(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				fmt::println(Logger::out("fix_labeling"),"{} chart(s) added to fix invalid boundaries",nb_chart_modified); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid corners")) {
				unsigned int nb_chart_modified = fix_invalid_corners(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				fmt::println(Logger::out("fix_labeling"),"{} chart(s) added to fix invalid corners",nb_chart_modified); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Remove invalid charts")) {
				remove_invalid_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				fmt::println(Logger::out("fix_labeling"),"invalid charts removed using gco"); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
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

			if(current_state_ == empty) {
				fmt::println(Logger::err("I/O"),"You need to import a triangle mesh before importing a labeling"); Logger::err("I/O").flush();
				return false;
			}

			if(!load_labeling(filename,mesh_,LABELING_ATTRIBUTE_NAME)) {
				fmt::println("load_labeling() not ok"); fflush(stdout);
				// Should the labeling be removed ?
				// If a labeling was already displayed, it should be restored...
				clear_scene_overlay();
				current_state_ = triangle_mesh;
				return false;
			}

			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);

			current_state_ = labeling;

			return true;
		}

        mesh_gfx_.set_mesh(nullptr);
        mesh_.clear(true,false); // keep_attributes=true is very important. else runtime error when re-importing files
        MeshIOFlags flags;
        if(!mesh_load(filename, mesh_, flags)) {
            current_state_ = empty; // added
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

		clear_scene_overlay();
        current_state_ = triangle_mesh;

        return true;
    }

	void GL_initialize() override {
        SimpleMeshApplicationExt::GL_initialize();
		init_rgba_colormap("labeling",6,1,labeling_colors_as_char_);
		init_rgba_colormap("validity",2,1,validity_colors_as_char_);
    }

    void update_static_labeling_graph(bool allow_boundaries_between_opposite_labels) {

		// compute charts, boundaries and corners of the labeling
		static_labeling_graph.fill_from(mesh_,LABELING_ATTRIBUTE_NAME,allow_boundaries_between_opposite_labels);
		nb_turning_points = static_labeling_graph.nb_turning_points();

		fmt::println(Logger::out("I/O"),"There are {} charts, {} corners and {} boundaries in this labeling.",static_labeling_graph.nb_charts(),static_labeling_graph.nb_corners(),static_labeling_graph.nb_boundaries());  Logger::out("I/O").flush();

		// static_labeling_graph.dump_to_file("StaticLabelingGraph.txt");

		clear_scene_overlay();

		valid_corners_group_index = new_points_group(corners_color_,true);
		invalid_corners_group_index = new_points_group(corners_color_,true);
		turning_points_group_index = new_points_group(turning_points_color_,true);
		X_boundaries_group_index = new_edges_group(&labeling_colors_as_float_[4*0],false); // axis X -> color of label 0 = +X
		Y_boundaries_group_index = new_edges_group(&labeling_colors_as_float_[4*2],false); // axis Y -> color of label 2 = +Y
		Z_boundaries_group_index = new_edges_group(&labeling_colors_as_float_[4*4],false); // axis Z -> color of label 4 = +Z
		invalid_boundaries_group_index = new_edges_group(&validity_colors_as_float_[4*0],false); // color of invalid LabelingGraph components
		valid_but_axisless_boundaries_group_index = new_edges_group(&validity_colors_as_float_[4*1],false); // color of valid LabelingGraph components

		for(std::size_t i = 0; i < static_labeling_graph.nb_corners(); ++i) {
			const double* coordinates = mesh_.vertices.point_ptr(
				static_labeling_graph.corners[i].vertex
			);
			add_point_to_group(static_labeling_graph.corners[i].is_valid ? valid_corners_group_index : invalid_corners_group_index,coordinates[0], coordinates[1], coordinates[2]);
		}

		for(index_t i : static_labeling_graph.non_monotone_boundaries) {
			for(index_t j : static_labeling_graph.boundaries[i].turning_points) {
				vec3 coordinates = mesh_vertex(mesh_, mesh_.facet_corners.vertex(static_labeling_graph.boundaries[i].halfedges[j].corner));
				add_point_to_group(turning_points_group_index,coordinates.x,coordinates.y,coordinates.z);
			}
		}

		std::size_t group_index;
		for(std::size_t i = 0; i < static_labeling_graph.nb_boundaries(); ++i) {
			const Boundary& boundary = static_labeling_graph.boundaries[i];
			for(const auto& be : boundary.halfedges) { // for each boundary edge of this boundary
				const double* coordinates_first_point = halfedge_vertex_from(mesh_,be).data();
				const double* coordinates_second_point = halfedge_vertex_to(mesh_,be).data();
				switch(boundary.axis) {
					case -1: // may or may not be invalid
						if(boundary.is_valid) {
							group_index = valid_but_axisless_boundaries_group_index;
						}
						else {
							group_index = invalid_boundaries_group_index;
						}
						break;
					case 0: // always valid
						group_index = X_boundaries_group_index;
						break;
					case 1: // always valid
						group_index = Y_boundaries_group_index;
						break;
					case 2: // always valid
						group_index = Z_boundaries_group_index;
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

		set_points_group_visibility(valid_corners_group_index,true);
		set_points_group_visibility(invalid_corners_group_index,true);
		set_edges_group_visibility(X_boundaries_group_index,true);
		set_edges_group_visibility(Y_boundaries_group_index,true);
		set_edges_group_visibility(Z_boundaries_group_index,true);
	}

private:

    bool allow_boundaries_between_opposite_labels_;
	unsigned char labeling_colors_as_char_[6*4]; // for the colormap (facets colors)
	float labeling_colors_as_float_[6*4]; // for ImGui controls & custom edges (boundaries colors)
	float corners_color_[4];
	float turning_points_color_[4];
	std::size_t nb_turning_points;
	unsigned char validity_colors_as_char_[2*4]; // for facets (invalid/valid charts)
	float validity_colors_as_float_[2*4]; // for edges and points (invalid/valid boundaries, invalid/valid corners)
	StaticLabelingGraph static_labeling_graph;
	State previous_state_, current_state_;
	int previous_labeling_visu_mode_, current_labeling_visu_mode_; // not a enum, to be used in ImGui
	std::size_t valid_corners_group_index;
	std::size_t invalid_corners_group_index;
	std::size_t turning_points_group_index;
	std::size_t X_boundaries_group_index;
	std::size_t Y_boundaries_group_index;
	std::size_t Z_boundaries_group_index;
	std::size_t invalid_boundaries_group_index;
	std::size_t valid_but_axisless_boundaries_group_index;
};

// print specialization of AutomaticPolycubeApp::State for {fmt}
template <> struct fmt::formatter<AutomaticPolycubeApp::State>: formatter<string_view> {
	// parse is inherited from formatter<string_view>.
	template <typename FormatContext>
	auto format(AutomaticPolycubeApp::State state, FormatContext& ctx) -> decltype(ctx.out()) const {
		switch (state) {
		case AutomaticPolycubeApp::State::empty			: return formatter<string_view>::format("empty", ctx);
		case AutomaticPolycubeApp::State::triangle_mesh	: return formatter<string_view>::format("triangle_mesh", ctx);
		case AutomaticPolycubeApp::State::labeling		: return formatter<string_view>::format("labeling", ctx);
		}
		return formatter<string_view>::format("?", ctx);
	}
};

int main(int argc, char** argv) {
    AutomaticPolycubeApp app;
	app.start(argc,argv);
    return 0;
}