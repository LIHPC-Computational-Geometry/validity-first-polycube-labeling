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
				unsigned int nb_chart_modified = remove_surrounded_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("fix_labeling"),"{} chart(s) modified with the surrounding label",nb_chart_modified); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid boundaries")) {
				unsigned int nb_chart_modified = fix_invalid_boundaries(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("fix_labeling"),"{} chart(s) added to fix invalid boundaries",nb_chart_modified); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid corners")) {
				unsigned int nb_chart_modified = fix_invalid_corners(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("fix_labeling"),"{} chart(s) added to fix invalid corners",nb_chart_modified); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Remove invalid charts")) {
				remove_invalid_charts(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("fix_labeling"),"invalid charts removed using gco"); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Auto fix validity")) {
				auto_fix_validity(mesh_,normals_,static_labeling_graph_,100,feature_edges_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			ImGui::SeparatorText("Monotonicity");

			if(ImGui::Button("Move boundaries near turning points")) {
				unsigned int nb_facets_modified = move_boundaries_near_turning_points(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_);
				fmt::println(Logger::out("monotonicity"),"label of {} facets modified",nb_facets_modified); Logger::out("monotonicity").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Straighten boundaries")) {
				FOR(b,static_labeling_graph_.boundaries.size()) {
					straighten_boundary(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,b);
					update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				}
			}

			if(ImGui::Button("Merge corners and turning points")) {
				FOR(b,static_labeling_graph_.non_monotone_boundaries.size()) {
					pull_closest_corner(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,b);
					update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				}
			}

			if(ImGui::Button("Auto fix monotonicity")) {
				auto_fix_monotonicity(mesh_,static_labeling_graph_,500,feature_edges_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			ImGui::SeparatorText("Refinement");

			ImGui::TextUnformatted("Other prefered label");
			ImGui::SameLine();
			if(ImGui::Button(selected_chart_mode_ ? "Pick a chart" : "Select chart")) {
				selected_chart_mode_ = true;
				selected_chart_= index_t(-1);
			}
			ImGui::SameLine();
			ImGui::BeginDisabled(selected_chart_mode_ || (selected_chart_==index_t(-1)));
			if(selected_chart_==index_t(-1))
				ImGui::TextUnformatted("current: none");
			else
				ImGui::Text("current: %d",selected_chart_);
			ImGui::EndDisabled();

			if(ImGui::Button("Trace contour")) {
				trace_contour(mesh_,normals_,LABELING_ATTRIBUTE_NAME,static_labeling_graph_,feature_edges_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
		}
	}

	void mouse_button_callback(int button, int action, int mods, int source) override {
		if(selected_chart_mode_) {
			selected_chart_mode_ = false;
			index_t facet_index = pick(MESH_FACETS);
			if( (facet_index == index_t(-1)) || (facet_index >= mesh_.facets.nb()) ) {
				selected_chart_ = index_t(-1);
			}
			else {
				// get the chart index of the picked facet
				selected_chart_ = static_labeling_graph_.facet2chart[facet_index];
				// get the label of the selected chart
				Attribute<index_t> label(mesh_.facets.attributes(),LABELING_ATTRIBUTE_NAME); // fetch labeling
				index_t previous_label = label[facet_index]; // we could also fetch it from StaticLabelingGraph::Chart::label, but we need the labeling anyway
				// for each facet of the current chart, associate another prefered label as much as possible
				for(auto f : static_labeling_graph_.charts[selected_chart_].facets) {
					const vec3& normal = normals_[f];
					// based on labeling.cpp nearest_label()
					std::array<std::pair<double,index_t>,6> per_label_weights = {
						std::make_pair(normal.x < 0.0 ? 0.0 : normal.x,     0), // +X
						std::make_pair(normal.x < 0.0 ? -normal.x : 0.0,    1), // -X
						std::make_pair(normal.y < 0.0 ? 0.0 : normal.y,     2), // +Y
						std::make_pair(normal.y < 0.0 ? -normal.y : 0.0,    3), // -Y
						std::make_pair(normal.z < 0.0 ? 0.0 : normal.z,     4), // +Z
						std::make_pair(normal.z < 0.0 ? -normal.z : 0.0,    5)  // -Z
					};
					std::sort(per_label_weights.begin(),per_label_weights.end());
					if(per_label_weights[5].second != previous_label) {
						// the prefered label is different from the previous label -> switch to the prefered label
						label[f] = per_label_weights[5].second;
					}
					else if( (per_label_weights[5].second == previous_label) && 
							 (per_label_weights[4].first == per_label_weights[3].first) &&
							 (per_label_weights[3].first == per_label_weights[2].first) &&
							 (per_label_weights[2].first == per_label_weights[1].first) &&
							 (per_label_weights[1].first == per_label_weights[0].first) ) {
						// we have to keep the previous label, because among other ones, no one is better than the others
						continue;
					}
					else {
						label[f] = per_label_weights[4].second; // take the second prefered label
					}
				}
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
		}
		else {
			LabelingViewerApp::mouse_button_callback(button,action,mods,source);
		}
	}

public:

	// return true if successfully found a valid labeling
	static bool auto_fix_validity(Mesh& mesh, std::vector<vec3>& normals, StaticLabelingGraph& slg, unsigned int max_nb_loop, const std::set<std::pair<index_t,index_t>>& feature_edges) {
		unsigned int nb_loops = 0;
		unsigned int nb_fixed_features = 0;
		std::set<std::array<std::size_t,7>> set_of_labeling_features_combinations_encountered;
		while(!slg.is_valid() && nb_loops <= max_nb_loop) { // until valid labeling OR too much steps
			nb_loops++;

			// as much as possible, remove isolated (surrounded) charts
			do {
				nb_fixed_features = remove_surrounded_charts(mesh,LABELING_ATTRIBUTE_NAME,slg);
				// update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);
			} while(nb_fixed_features != 0);

			if(slg.is_valid())
				return true;

			fix_invalid_boundaries(mesh,LABELING_ATTRIBUTE_NAME,slg);
			// update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);

			if(slg.is_valid())
				return true;

			fix_invalid_corners(mesh,normals,LABELING_ATTRIBUTE_NAME,slg);
			// update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);

			if(slg.is_valid())
				return true;

			set_of_labeling_features_combinations_encountered.clear();
			set_of_labeling_features_combinations_encountered.insert({
				slg.nb_charts(),
				slg.nb_boundaries(),
				slg.nb_corners(),
				slg.nb_invalid_charts(),
				slg.nb_invalid_boundaries(),
				slg.nb_invalid_corners(),
				slg.nb_turning_points()
			});

			while(1) {
				if(remove_invalid_charts(mesh,normals,LABELING_ATTRIBUTE_NAME,slg)) {
					// all invalid charts cannot be removed because of feature edges
					// No need to update the labeling graph because the labeling didn't change
					if((slg.nb_invalid_boundaries() == 0) && (slg.nb_invalid_corners() == 0)) {
						// these charts are the only invalid parts of the labeling graph
						fmt::println(Logger::err("fix_labeling"),"the labeling validity cannot be reached because remaining invalid charts are surrounded by feature edges. TODO new operator to fix them"); Logger::err("fix_labeling").flush();
						return false;
					}
					// else : continue and expect the other operators to fix invalid boundaries and corners
				}
				// update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);

				if(slg.is_valid())
					return true;

				std::array<std::size_t,7> features_combination = {
					slg.nb_charts(),
					slg.nb_boundaries(),
					slg.nb_corners(),
					slg.nb_invalid_charts(),
					slg.nb_invalid_boundaries(),
					slg.nb_invalid_corners(),
					slg.nb_turning_points()
				};

				if(VECTOR_CONTAINS(set_of_labeling_features_combinations_encountered,features_combination)) { // wa can use VECTOR_CONTAINS() on sets because they also have find(), cbegin() and cend()
					// we backtracked
					// There is probably small charts that we can remove to help the fixing routine
					remove_charts_around_invalid_boundaries(mesh,normals,LABELING_ATTRIBUTE_NAME,slg);
					slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);
					break; // go back to the beginning of the loop, with other fix operators
				}
				else {
					set_of_labeling_features_combinations_encountered.insert(features_combination); // store the current combination of number of features
				}
			}
			
		}

		if(!slg.is_valid()) {
			fmt::println(Logger::out("fix_labeling"),"auto fix validity stopped (max nb loops reached), no valid labeling found"); Logger::out("fix_labeling").flush();
			return false;
		}
		
		return true;
	}

	static bool auto_fix_monotonicity(Mesh& mesh, StaticLabelingGraph& slg, unsigned int max_nb_steps, const std::set<std::pair<index_t,index_t>>& feature_edges) {
		unsigned int nb_steps = 0;
		while(!slg.non_monotone_boundaries.empty() && nb_steps <= max_nb_steps) {
			move_boundaries_near_turning_points(mesh,LABELING_ATTRIBUTE_NAME,slg);
			// update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);
			remove_surrounded_charts(mesh,LABELING_ATTRIBUTE_NAME,slg);
			// update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			slg.fill_from(mesh,LABELING_ATTRIBUTE_NAME,slg.is_allowing_boundaries_between_opposite_labels(),feature_edges);
		}

		if(!slg.non_monotone_boundaries.empty()) {
			fmt::println(Logger::out("fix_labeling"),"auto fix monotonicity stopped (max nb steps reached), it remains non-monotone boundaries"); Logger::out("fix_labeling").flush();
			return false;
		}

		return true;
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

	bool allow_boundaries_between_opposite_labels = false;
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
	slg.fill_from(M,LABELING_ATTRIBUTE_NAME,allow_boundaries_between_opposite_labels,feature_edges);

	//////////////////////////////////////////////////
	// Validity & monotonicity correction
	//////////////////////////////////////////////////

	if(AutomaticPolycubeApp::auto_fix_validity(M,normals,slg,100,feature_edges)) {
		// auto-fix the monotonicity only if the validity was fixed
		AutomaticPolycubeApp::auto_fix_monotonicity(M,slg,500,feature_edges);
	}

	//////////////////////////////////////////////////
	// Write output labeling
	//////////////////////////////////////////////////

	fmt::println(Logger::out("I/O"),"Writing {}...",filenames[1]); Logger::out("I/O").flush();
	save_labeling(filenames[1],M,LABELING_ATTRIBUTE_NAME);
	fmt::println(Logger::out("I/O"),"Done"); Logger::out("I/O").flush();
	
    return 0;
}