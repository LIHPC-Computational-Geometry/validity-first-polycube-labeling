#include <set>
#include <array>

#include "LabelingViewerApp.h"

class AutomaticPolycubeApp : public LabelingViewerApp {
public:

    AutomaticPolycubeApp() : LabelingViewerApp("automatic_polycube") {}

private:

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(current_state_ == LabelingViewerApp::State::labeling) {
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