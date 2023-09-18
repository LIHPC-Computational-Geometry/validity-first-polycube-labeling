// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/graphcut_labeling.cpp
// using Geogram instead of libigl, and with a GUI.

// TODO if LabelingViewerApp::current_state_ is LabelingViewerApp::State::triangle_mesh,
// compute graph-cut labeling with default paremeters value (instead of manual, naive labeling)

// TODO do not lauch GUI if both parameters are given by CLI

#include "LabelingViewerApp.h"
#include "labeling.h"

#include <array>

#include "GraphCutLabeling.h"

class GraphCutLabelingApp : public LabelingViewerApp {
public:

    GraphCutLabelingApp() : LabelingViewerApp("graphcut_labeling") {
		// PolyCut [1, Section 3.1, page 5] and Evocube [2, section 4, p. 7] use a fidelity/compatness ratio of 3 for the first try
		// [1] Livesu, Vining, Sheffer, Gregson, Scateni, "Polycut: Monotone graph-cuts for polycube base-complex construction", ACM Trans. on Graphics, 2013
		// [2] Dumery, Protais, Mestrallet, Bourcier, Ledoux, "Evocube: a Genetic Labeling Framework for Polycube-Maps", Computer Graphics Forum, 2022
		compactness_coeff = 1;
		fidelity_coeff = 3;
		_smooth_cost__fill(smooth_cost); // default value filled with dedicated function
	}

protected:

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(current_state_ == LabelingViewerApp::State::labeling) {
			ImGui::Separator();

			ImGui::Text("Graph-cut parameters");

			ImGui::InputInt("Compactness", &compactness_coeff);
			ImGui::InputInt("Fidelity", &fidelity_coeff);

			ImGui::Text("Smooth cost");
			if (ImGui::BeginTable("Smooth cost table", 7, ImGuiTableFlags_Borders)) { // 7 columns
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::TextUnformatted(" ");
				for (int column = 0; column < 6; column++) {
					ImGui::TableSetColumnIndex(column+1);
					ImGui::TextUnformatted(LABEL2STR(column));
				}
				for (int row = 0; row < 6; row++)
				{
					ImGui::TableNextRow();
					ImGui::TableSetColumnIndex(0);
					ImGui::TextUnformatted(LABEL2STR(row));
					for (int column = 0; column < 6; column++)
					{
						ImGui::TableSetColumnIndex(column+1);
						char label[32];
						sprintf(label, "##smooth:%d,%d", column, row);
						ImGui::SetNextItemWidth(100.0f);
						ImGui::InputInt(label,&smooth_cost[column*6+row]);
					}
				}
				ImGui::EndTable();
			}

			if(ImGui::Button("Compute solution")) {
				graphcut_labeling(mesh_,LABELING_ATTRIBUTE_NAME,compactness_coeff,fidelity_coeff,smooth_cost);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
		}
	}

	int compactness_coeff;
	int fidelity_coeff;
	std::array<int,6*6> smooth_cost;
};

int main(int argc, char** argv) {
    GraphCutLabelingApp app;
	app.start(argc,argv);
    return 0;
}