// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/graphcut_labeling.cpp
// using Geogram instead of libigl, and with a GUI.

// TODO if LabelingViewerApp::current_state_ is LabelingViewerApp::State::triangle_mesh,
// compute graph-cut labeling with default paremeters value (instead of manual, naive labeling)

// TODO do not lauch GUI if both parameters are given by CLI

#include "LabelingViewerApp.h"
#include "labeling.h"

class GraphCutLabelingApp : public LabelingViewerApp {
public:

    GraphCutLabelingApp() : LabelingViewerApp("graphcut_labeling") {
		// PolyCut [1, Section 3.1, page 5] and Evocube [2, section 4, p. 7] use a fidelity/compatness ratio of 3 for the first try
		// [1] Livesu, Vining, Sheffer, Gregson, Scateni, "Polycut: Monotone graph-cuts for polycube base-complex construction", ACM Trans. on Graphics, 2013
		// [2] Dumery, Protais, Mestrallet, Bourcier, Ledoux, "Evocube: a Genetic Labeling Framework for Polycube-Maps", Computer Graphics Forum, 2022
		compactness_coeff = 1;
		fidelity_coeff = 3;
	}

protected:

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(current_state_ == LabelingViewerApp::State::labeling) {
			ImGui::Separator();

			ImGui::Text("Graph-cut parameters");

			ImGui::InputInt("Compatness", &compactness_coeff);
			ImGui::InputInt("Fidelity", &fidelity_coeff);

			if(ImGui::Button("Compute solution")) {
				graphcut_labeling(mesh_,LABELING_ATTRIBUTE_NAME,compactness_coeff,fidelity_coeff);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
		}
	}

	int compactness_coeff;
	int fidelity_coeff;
};

int main(int argc, char** argv) {
    GraphCutLabelingApp app;
	app.start(argc,argv);
    return 0;
}