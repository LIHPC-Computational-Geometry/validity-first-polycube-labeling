// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/graphcut_labeling.cpp
// using Geogram instead of libigl, and with a GUI.

// TODO if LabelingViewerApp::state_ is LabelingViewerApp::State::triangle_mesh,
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
		selected_chart = index_t(-1);
		selected_chart_mode = false;
	}

protected:

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(state_ == LabelingViewerApp::State::labeling) {
			ImGui::SeparatorText("Global edition");

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
						ImGui::InputInt(label,&smooth_cost[(std::array<int, 36>::size_type) (column*6+row)]);
					}
				}
				ImGui::EndTable();
			}

			if(ImGui::Button("Compute solution")) {
				graphcut_labeling(mesh_,LABELING_ATTRIBUTE_NAME,compactness_coeff,fidelity_coeff,smooth_cost);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			ImGui::SeparatorText("Per-chart edition");

			if(ImGui::Button("Select chart")) {
				selected_chart_mode = true;
			}
			ImGui::SameLine();
			ImGui::BeginDisabled(selected_chart_mode || (selected_chart==index_t(-1)));
			if(selected_chart==index_t(-1))
				ImGui::TextUnformatted("current: none");
			else
				ImGui::Text("current: %d",selected_chart);
			ImGui::EndDisabled();
		}
	}

	void mouse_button_callback(int button, int action, int mods, int source) override {
		if(selected_chart_mode) {
			selected_chart_mode = false;
			index_t x = index_t(cursor_pos_.x), y = index_t(cursor_pos_.y); // double to integer conversion of current cursor position
			if(x >= get_width() || y >= get_height()) { // if cursor out of the window
				return;
			}
			y = get_height()-1-y; // change Y axis orientation. glReadPixels() wants pixel coordinates from bottom-left corner
			mesh_gfx()->set_picking_mode(MESH_FACETS); // instead of rendering colors, mesh_gfx will render facet indices
			draw_scene(); // rendering
			// read the index of the picked element using glReadPixels()
			Memory::byte picked_mesh_element_as_pixel[4];
			glPixelStorei(GL_PACK_ALIGNMENT, 1);
			glPixelStorei(GL_PACK_ROW_LENGTH, 1);
			glReadPixels(
				GLint(x),GLint(y),1,1,GL_RGBA,GL_UNSIGNED_BYTE,picked_mesh_element_as_pixel
			);
			mesh_gfx()->set_picking_mode(MESH_NONE); // go back to color rendering mode
			// decode facet index from pixel color
			index_t facet_index = index_t(picked_mesh_element_as_pixel[0])        |
								 (index_t(picked_mesh_element_as_pixel[1]) << 8)  |
								 (index_t(picked_mesh_element_as_pixel[2]) << 16) |
								 (index_t(picked_mesh_element_as_pixel[3]) << 24);
			if( (facet_index == index_t(-1)) || (facet_index >= mesh_.facets.nb()) )
				selected_chart = index_t(-1);
			else
				selected_chart = static_labeling_graph_.facet2chart[facet_index];
		}
		else {
			LabelingViewerApp::mouse_button_callback(button,action,mods,source);
		}
	}

	int compactness_coeff;
	int fidelity_coeff;
	std::array<int,6*6> smooth_cost;
	index_t selected_chart;
	bool selected_chart_mode;
};

int main(int argc, char** argv) {
    GraphCutLabelingApp app;
	app.start(argc,argv);
    return 0;
}