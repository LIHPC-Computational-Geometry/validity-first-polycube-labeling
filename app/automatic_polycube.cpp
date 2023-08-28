#include <geogram/basic/attributes.h>

#include <set>
#include <array>
#include <algorithm>			// for std::min(), std::max()

#include "containers.h"			// for std_dev()
#include "LabelingViewerApp.h"
#include "LabelingGraph.h"

#define VIEW_CHARTS_TO_REFINE	((VIEW_INVALID_CORNERS)+1)

#define STDDEV_BASED_CHART_REFINEMENT

class AutomaticPolycubeApp : public LabelingViewerApp {
public:

    AutomaticPolycubeApp() : LabelingViewerApp("automatic_polycube") {}

private:

	void draw_scene() override {

		// handle new option "View charts to refine"
		if(
		(current_state_ == LabelingViewerApp::State::labeling) && 
		(previous_labeling_visu_mode_ != current_labeling_visu_mode_) &&
		(current_labeling_visu_mode_ == VIEW_CHARTS_TO_REFINE)) {
			Attribute<index_t> labeling(mesh_.facets.attributes(),LABELING_ATTRIBUTE_NAME); // access the labeling
			vec3 facet_normal(0.0,0.0,0.0);
#ifdef STDDEV_BASED_CHART_REFINEMENT
			Attribute<double> per_chart_stddev_validity(mesh_.facets.attributes(),"per_chart_stddev_validity"); // create new facet attribute
			std::vector<index_t> per_chart_fidelity_values;
			double global_max; // will store the max std dev of all facets
			FOR(c,static_labeling_graph.nb_charts()) { // for each chart
				double stddev = 0.0;
				per_chart_fidelity_values.resize(static_labeling_graph.charts[c].facets.size());
				index_t count_facets = 0;
				for(auto f : static_labeling_graph.charts[c].facets) { // for each facet in this chart
					facet_normal = normalize(Geom::mesh_facet_normal(mesh_,f));
					per_chart_fidelity_values[count_facets] = (GEO::dot(facet_normal,label2vector[labeling[f]]) - 1.0)/0.2; // compute fidelity
					count_facets++;
				}
				stddev = std_dev(per_chart_fidelity_values.begin(),per_chart_fidelity_values.end());
				// write this value for all facets on the current chart
				for(auto f : static_labeling_graph.charts[c].facets) {
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
			FOR(c,static_labeling_graph.nb_charts()) { // for each chart
				min_fidelity = 0.0; // will store the min fidelity of all facets
				for(auto f : static_labeling_graph.charts[c].facets) { // for each facet in this chart
					facet_normal = normalize(Geom::mesh_facet_normal(mesh_,f));
					fidelity = (GEO::dot(facet_normal,label2vector[labeling[f]]) - 1.0)/0.2;
					min_fidelity = std::min(fidelity,min_fidelity);
				}
				// write this value for all facets on the current chart
				for(auto f : static_labeling_graph.charts[c].facets) { // for each facet in this chart
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
			current_colormap_texture_ = colormaps_[COLORMAP_BLACK_WHITE].texture;
			attribute_subelements_ = MESH_FACETS;
			
			// points in overlay
			points_groups_show_only({}); // show none
			// edges in overlay
			set_edges_group_color(X_boundaries_group_index,labeling_colors_.color_as_floats(0)); // axis X -> color of label 0 = +X
			set_edges_group_color(Y_boundaries_group_index,labeling_colors_.color_as_floats(2)); // axis Y -> color of label 2 = +Y
			set_edges_group_color(Z_boundaries_group_index,labeling_colors_.color_as_floats(4)); // axis Z -> color of label 4 = +Z
			edges_groups_show_only({X_boundaries_group_index, Y_boundaries_group_index, Z_boundaries_group_index});

			previous_labeling_visu_mode_ = current_labeling_visu_mode_;
			SimpleMeshApplicationExt::draw_scene();
		}
		else {
			// other visualization options managed in LabelingViewerApp
			LabelingViewerApp::draw_scene();
		}
	}

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(current_state_ == LabelingViewerApp::State::labeling) {
			ImGui::RadioButton("View charts to refine",&current_labeling_visu_mode_,VIEW_CHARTS_TO_REFINE);

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

			if(ImGui::Button("Auto fix validity")) {
				if(auto_fix_validity(100)) {
					fmt::println(Logger::out("fix_labeling"),"auto fix validity found a valid labeling"); Logger::out("fix_labeling").flush();
				}
				// else : auto_fix_validity() already printed a message
			}

			ImGui::Separator();

			ImGui::Text("Monotonicity");

			if(ImGui::Button("Move boundaries near turning points")) {
				unsigned int nb_facets_modified = move_boundaries_near_turning_points(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				fmt::println(Logger::out("monotonicity"),"label of {} facets modified",nb_facets_modified); Logger::out("monotonicity").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Auto fix monotonicity")) {
				if(auto_fix_monotonicity(500)) {
					fmt::println(Logger::out("fix_labeling"),"auto fix monotonicity found a labeling with all-monotone boundaries"); Logger::out("fix_labeling").flush();
				}
				// else : auto_fix_monotonicity() already printed a message
			}
		}
	}

	// return true if successfully found a valid labeling
	bool auto_fix_validity(unsigned int max_nb_loop) {
		unsigned int nb_loops = 0;
		unsigned int nb_fixed_features = 0;
		std::set<std::array<std::size_t,7>> set_of_labeling_features_combinations_encountered;
		while(!static_labeling_graph.is_valid() && nb_loops <= max_nb_loop) { // until valid labeling OR too much steps
			nb_loops++;

			// as much as possible, remove isolated (surrounded) charts
			do {
				nb_fixed_features = remove_surrounded_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			} while(nb_fixed_features != 0);

			if(static_labeling_graph.is_valid())
				return true;

			fix_invalid_boundaries(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);

			if(static_labeling_graph.is_valid())
				return true;

			fix_invalid_corners(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);

			if(static_labeling_graph.is_valid())
				return true;

			set_of_labeling_features_combinations_encountered.clear();
			set_of_labeling_features_combinations_encountered.insert({
				static_labeling_graph.nb_charts(),
				static_labeling_graph.nb_boundaries(),
				static_labeling_graph.nb_corners(),
				static_labeling_graph.nb_invalid_charts(),
				static_labeling_graph.nb_invalid_boundaries(),
				static_labeling_graph.nb_invalid_corners(),
				static_labeling_graph.nb_turning_points()
			});

			while(1) {
				remove_invalid_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);

				if(static_labeling_graph.is_valid())
					return true;

				std::array<std::size_t,7> features_combination = {
					static_labeling_graph.nb_charts(),
					static_labeling_graph.nb_boundaries(),
					static_labeling_graph.nb_corners(),
					static_labeling_graph.nb_invalid_charts(),
					static_labeling_graph.nb_invalid_boundaries(),
					static_labeling_graph.nb_invalid_corners(),
					static_labeling_graph.nb_turning_points()
				};

				if(VECTOR_CONTAINS(set_of_labeling_features_combinations_encountered,features_combination)) { // wa can use VECTOR_CONTAINS() on sets because they also have find(), cbegin() and cend()
					// we backtracked
					break; // go back to the beginning of the loop, with other fix operators
				}
				else {
					set_of_labeling_features_combinations_encountered.insert(features_combination); // store the current combination of number of features
				}
			}
			
		}

		if(!static_labeling_graph.is_valid()) {
			fmt::println(Logger::out("fix_labeling"),"auto fix validity stopped (max nb loops reached), no valid labeling found"); Logger::out("fix_labeling").flush();
			return false;
		}
		
		return true;
	}

	bool auto_fix_monotonicity(unsigned int max_nb_steps) {
		unsigned int nb_steps = 0;
		while(!static_labeling_graph.non_monotone_boundaries.empty() && nb_steps <= max_nb_steps) {
			move_boundaries_near_turning_points(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			remove_surrounded_charts(mesh_,LABELING_ATTRIBUTE_NAME,static_labeling_graph);
			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
		}

		if(!static_labeling_graph.non_monotone_boundaries.empty()) {
			fmt::println(Logger::out("fix_labeling"),"auto fix monotonicity stopped (max nb steps reached), it remains non-monotone boundaries"); Logger::out("fix_labeling").flush();
			return false;
		}

		return true;
	}
};

int main(int argc, char** argv) {
    AutomaticPolycubeApp app;
	app.start(argc,argv);
    return 0;
}