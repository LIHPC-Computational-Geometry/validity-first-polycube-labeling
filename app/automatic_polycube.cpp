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
		}
	}
};

int main(int argc, char** argv) {
    AutomaticPolycubeApp app;
	app.start(argc,argv);
    return 0;
}