/*
 * Tricky CLI arguments format with Geogram :
 * Filenames first, then options without dashes prefixes
 * 
 * With GUI:
 *     ./automatic_polycube surface_mesh.obj gui=true
 * 
 * Without GUI:
 *     ./automatic_polycube surface_mesh.obj output_labeling.txt gui=false
 */

#include <geogram/basic/attributes.h>
#include <geogram/basic/command_line.h>			// for declare_arg(), parse(), get_arg_bool()
#include <geogram/basic/command_line_args.h>	// for import_arg_group()
#include <geogram/mesh/mesh_io.h>				// for mesh_load()
#include <geogram/basic/numeric.h>				// for min_float64()

#include <set>
#include <array>
#include <algorithm>			// for std::min(), std::max(), std::min_element()
#include <tuple>				// for std::tie()
#include <deque>

#include "containers.h"			// for std_dev()
#include "LabelingViewerApp.h"
#include "LabelingGraph.h"

#define VIEW_CHARTS_TO_REFINE	((VIEW_INVALID_CORNERS)+1)

#define STDDEV_BASED_CHART_REFINEMENT

class AutomaticPolycubeApp : public LabelingViewerApp {
public:

    AutomaticPolycubeApp() : LabelingViewerApp("automatic_polycube"),
							 selected_chart_(index_t(-1)),
							 selected_chart_mode_(false)
							 {}

protected:

	void labeling_visu_mode_transition(int new_mode) override {
		// handle new option "View charts to refine"
		if(new_mode == VIEW_CHARTS_TO_REFINE) {
			Attribute<index_t> labeling(mesh_.facets.attributes(),LABELING_ATTRIBUTE_NAME); // access the labeling
			vec3 facet_normal(0.0,0.0,0.0);
#ifdef STDDEV_BASED_CHART_REFINEMENT
			Attribute<double> per_chart_stddev_validity(mesh_.facets.attributes(),"per_chart_stddev_validity"); // create new facet attribute
			std::vector<double> per_chart_fidelity_values;
			double global_max = Numeric::min_float64(); // will store the max std dev of all facets
			FOR(c,static_labeling_graph_.nb_charts()) { // for each chart
				double stddev = 0.0;
				per_chart_fidelity_values.resize(static_labeling_graph_.charts[c].facets.size());
				index_t count_facets = 0;
				for(auto f : static_labeling_graph_.charts[c].facets) { // for each facet in this chart
					per_chart_fidelity_values[count_facets] = (GEO::dot(normals_[f],label2vector[labeling[f]]) - 1.0)/0.2; // compute fidelity
					count_facets++;
				}
				stddev = std_dev(per_chart_fidelity_values.begin(),per_chart_fidelity_values.end());
				// write this value for all facets on the current chart
				for(auto f : static_labeling_graph_.charts[c].facets) {
					per_chart_stddev_validity[f] = stddev;
				}
				global_max = std::max(global_max,stddev);
			}

			attribute_ = "facets.per_chart_stddev_validity";
			attribute_name_ = "per_chart_stddev_validity";
			attribute_min_ = (float) global_max;	// in black
			attribute_max_ = 0;						// in white
#else
			Attribute<double> per_chart_min_validity(mesh_.facets.attributes(),"per_chart_min_validity"); // create new facet attribute
			double global_min = 0.0, min_fidelity = 0.0, fidelity = 0.0;
			FOR(c,static_labeling_graph_.nb_charts()) { // for each chart
				min_fidelity = 0.0; // will store the min fidelity of all facets
				for(auto f : static_labeling_graph_.charts[c].facets) { // for each facet in this chart
					fidelity = (GEO::dot(normals_[f],label2vector[labeling[f]]) - 1.0)/0.2;
					min_fidelity = std::min(fidelity,min_fidelity);
				}
				// write this value for all facets on the current chart
				for(auto f : static_labeling_graph_.charts[c].facets) { // for each facet in this chart
					per_chart_min_validity[f] = min_fidelity;
				}
				global_min = std::min(global_min,min_fidelity);
			}

			attribute_ = "facets.per_chart_min_validity";
			attribute_name_ = "per_chart_min_validity";
			attribute_min_ = (float) global_min;	// in black
			attribute_max_ = 0;						// in white
#endif
			show_mesh_ = false;
			lighting_ = false;
			show_attributes_ = true;
			current_colormap_index_ = COLORMAP_BLACK_WHITE;
			attribute_subelements_ = MESH_FACETS;
			
			// points in overlay
			points_groups_show_only({}); // show none
			// edges in overlay
			set_edges_group_color(X_boundaries_group_index_,labeling_colors_.color_as_floats(0)); // axis X -> color of label 0 = +X
			set_edges_group_color(Y_boundaries_group_index_,labeling_colors_.color_as_floats(2)); // axis Y -> color of label 2 = +Y
			set_edges_group_color(Z_boundaries_group_index_,labeling_colors_.color_as_floats(4)); // axis Z -> color of label 4 = +Z
			edges_groups_show_only({X_boundaries_group_index_, Y_boundaries_group_index_, Z_boundaries_group_index_});
			labeling_visu_mode_ = VIEW_CHARTS_TO_REFINE;
		}
		else {
			LabelingViewerApp::labeling_visu_mode_transition(new_mode);
		}
	}

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(state_ == LabelingViewerApp::State::labeling) {
			if(ImGui::RadioButton("View charts to refine",&labeling_visu_mode_,VIEW_CHARTS_TO_REFINE)) {
				labeling_visu_mode_transition(VIEW_CHARTS_TO_REFINE);
			}

			ImGui::SeparatorText("Fix labeling");

			if(ImGui::Button("Remove surrounded charts")) {
				size_t nb_invalid_charts_processed = remove_surrounded_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid charts processed",nb_invalid_charts_processed); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid boundaries (as much as possible)")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				size_t nb_invalid_boundaries_processed = fix_as_much_invalid_boundaries_as_possible(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,normals_,feature_edges_,adj_facets_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid boundaries processed",nb_invalid_boundaries_processed); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Remove charts around invalid boundaries")) {
				remove_charts_around_invalid_boundaries(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid corners (as much as possible)")) {
				size_t nb_invalid_corners_processed = fix_as_much_invalid_corners_as_possible(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,feature_edges_,static_labeling_graph_.facet2chart,adj_facets_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid corners processed",nb_invalid_corners_processed); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Remove invalid charts")) {
				remove_invalid_charts(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("fix_labeling"),"invalid charts removed using gco"); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Increase chart valence")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				bool chart_processed = increase_chart_valence(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,adj_facets_,feature_edges_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid charts processed",size_t(chart_processed)); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Pick the first invalid chart surrounded by feature edges\nand trace a new adjacent chart to increase the valence of the first one");

			ImGui::PushStyleColor(ImGuiCol_Button, 			ImVec4(0.6f, 0.9f, 0.6f, 1.0f)); // green
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered,	ImVec4(0.3f, 0.8f, 0.3f, 1.0f)); // darker green
			ImGui::PushStyleColor(ImGuiCol_ButtonActive,	ImVec4(0.0f, 0.8f, 0.0f, 1.0f));
			if(ImGui::Button("Auto fix validity")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				auto_fix_validity(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,100,feature_edges_,normals_,adj_facets_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::PopStyleColor(3);

			ImGui::SeparatorText("Monotonicity");

			if(ImGui::Button("Move boundaries near turning points")) {
				size_t nb_turning_points_processed = move_boundaries_near_turning_points(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,feature_edges_);
				fmt::println(Logger::out("monotonicity"),"{} turning-points processed",nb_turning_points_processed); Logger::out("monotonicity").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Straighten boundaries")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				straighten_boundaries(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,adj_facets_,feature_edges_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Re-draw boundaries so that they are straighter");

			if(ImGui::Button("Move corners")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				move_corners(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,feature_edges_,adj_facets_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Merge a turning-points and its closest corners")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				bool a_non_monotone_boundary_was_processed = merge_a_turning_point_and_its_closest_corner(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,feature_edges_,adj_facets_);
				fmt::println(Logger::out("monotonicity"),"{} non-monotone boundaries processed",(size_t) a_non_monotone_boundary_was_processed); Logger::out("monotonicity").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Parse non-monotone boundaries, and if one of them has a turning-point on a feature edge, try to pull the closest corner so that they coincide");

			if(ImGui::Button("Join turning-points pair with new chart")) {
				// compute vertex-to-facet adjacency if not already done
				if(adj_facets_.empty()) {
					compute_adjacent_facets_of_vertices(mesh_,adj_facets_);
				}
				join_turning_points_pair_with_new_chart(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,normals_,feature_edges_,adj_facets_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			ImGui::PushStyleColor(ImGuiCol_Button, 			ImVec4(0.6f, 0.9f, 0.6f, 1.0f)); // green
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered,	ImVec4(0.3f, 0.8f, 0.3f, 1.0f)); // darker green
			ImGui::PushStyleColor(ImGuiCol_ButtonActive,	ImVec4(0.0f, 0.8f, 0.0f, 1.0f));
			if(ImGui::Button("Auto fix monotonicity")) {
				auto_fix_monotonicity(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,adj_facets_,feature_edges_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::PopStyleColor(3);
		}
	}

protected:

	index_t selected_chart_;
	unsigned char selected_chart_mode_;
};

int main(int argc, char** argv) {

	std::vector<std::string> filenames;
	GEO::initialize();
	CmdLine::import_arg_group("standard"); // strangely, this line is required for the Logger to work on Release mode when gui=false...
	CmdLine::declare_arg("gui", true, "Show the graphical user interface");
	if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"<input_surface_mesh> <output_labeling>"
		))
	{
		return 1;
	}
	
	if(CmdLine::get_arg_bool("gui")) { // if GUI mode
		AutomaticPolycubeApp app;
		app.start(argc,argv);
		return 0;
	}
	// else: no GUI, auto-process the input mesh

	if(filenames.size() == 0) { // missing filenames[0], that is <input_surface_mesh> argument
		fmt::println(Logger::err("I/O"),"When gui=false, the input filename must be given as CLI argument"); Logger::err("I/O").flush();
		return 1;
	}

	if(filenames.size() < 2) { // missing filenames[1], that is <output_labeling> argument
		fmt::println(Logger::warn("I/O"),"The output filename was not provided, using default labeling.txt"); Logger::warn("I/O").flush();
		filenames.push_back("labeling.txt"); // default output filename
	}

	bool allow_boundaries_between_opposite_labels = true;
	Mesh M;
	if(!mesh_load(filenames[0],M)) {
		fmt::println(Logger::err("I/O"),"Unable to load mesh from {}",filenames[0]);
		return 1;
	}

	if(M.facets.nb() == 0) {
		fmt::println(Logger::err("I/O"),"Input mesh {} has no facets",filenames[0]);
		return 1;
	}

	if(!M.facets.are_simplices()) {
		fmt::println(Logger::err("I/O"),"Input mesh {} is not a triangle mesh",filenames[0]);
		return 1;
	}

	if(M.cells.nb() != 0) {
		fmt::println(Logger::err("I/O"),"Input mesh {} is a volume mesh (#cells is not zero)",filenames[0]);
		return 1;
	}

	//////////////////////////////////////////////////
	// Geometry preprocessing, see LabelingViewerApp::load()
	//////////////////////////////////////////////////

	// TODO mesh rotation according to the PCA, see rotate_mesh_according_to_principal_axes() in geometry.h

	// ensure facet normals are outward
	if(facet_normals_are_inward(M)) {
		flip_facet_normals(M);
		fmt::println(Logger::warn("normals dir."),"Facet normals of the input mesh were inward");
		fmt::println(Logger::warn("normals dir."),"You should flip them and update the file for consistency.");
		Logger::warn("normals dir.").flush();
	}

	// compute facet normals
	std::vector<vec3> normals(M.facets.nb());
	FOR(f,M.facets.nb()) {
		normals[f] = normalize(Geom::mesh_facet_normal(M,f));
	}

	// remove feature edges on edge with small angle
	std::vector<std::vector<index_t>> adj_facets; // for each vertex, store adjacent facets. no ordering
	remove_feature_edges_with_low_dihedral_angle(M,adj_facets);
	// store them as a set
	std::set<std::pair<index_t,index_t>> feature_edges;
	transfer_feature_edges(M,feature_edges);

	//////////////////////////////////////////////////
	// Compute initial labeling
	//////////////////////////////////////////////////

	tweaked_naive_labeling(M,normals,LABELING_ATTRIBUTE_NAME);

	//////////////////////////////////////////////////
	// Construct charts, boundaries & corners
	//////////////////////////////////////////////////

	StaticLabelingGraph slg;
	slg.fill_from(M,LABELING_ATTRIBUTE_NAME,feature_edges,allow_boundaries_between_opposite_labels);

	//////////////////////////////////////////////////
	// Validity & monotonicity correction
	//////////////////////////////////////////////////

	if(auto_fix_validity(M,normals,LABELING_ATTRIBUTE_NAME,slg,100,feature_edges,normals,adj_facets)) {
		// auto-fix the monotonicity only if the validity was fixed
		auto_fix_monotonicity(M,LABELING_ATTRIBUTE_NAME,slg,adj_facets,feature_edges);
	}

	//////////////////////////////////////////////////
	// Write output labeling
	//////////////////////////////////////////////////

	fmt::println(Logger::out("I/O"),"Writing {}...",filenames[1]); Logger::out("I/O").flush();
	save_labeling(filenames[1],M,LABELING_ATTRIBUTE_NAME);
	fmt::println(Logger::out("I/O"),"Done"); Logger::out("I/O").flush();
	
    return 0;
}