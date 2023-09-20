// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/graphcut_labeling.cpp
// using Geogram instead of libigl, and with a GUI.

// TODO do not lauch GUI if both parameters are given by CLI

#include "LabelingViewerApp.h"
#include "labeling.h"

#include <array>
#include <algorithm> // for std::max_element(), std::min_element()

#include "GraphCutLabeling.h"

class GraphCutLabelingApp : public LabelingViewerApp {
public:

    GraphCutLabelingApp() : LabelingViewerApp("graphcut_labeling") {
		// PolyCut [1, Section 3.1, page 5] and Evocube [2, section 4, p. 7] use a fidelity/compatness ratio of 3 for the first try
		// [1] Livesu, Vining, Sheffer, Gregson, Scateni, "Polycut: Monotone graph-cuts for polycube base-complex construction", ACM Trans. on Graphics, 2013
		// [2] Dumery, Protais, Mestrallet, Bourcier, Ledoux, "Evocube: a Genetic Labeling Framework for Polycube-Maps", Computer Graphics Forum, 2022
		compactness_coeff_ = 1;
		fidelity_coeff_ = 3;
		// fill smooth_cost_. equivalent to what GraphCutLabeling::smooth_cost__set__default() does.
		FOR(label1,6) {
        	FOR(label2,6) {
				// same label = very smooth edge, different label = less smooth
				smooth_cost_[label1+label2*6] = (label1==label2) ? 0 : 1;
			}
		}
		selected_chart_ = index_t(-1);
		selected_chart_mode_ = false;
		// new_data_cost_ auto-initialized to 0
		new_data_cost_upper_bound_ = 1000.0f;
	}

protected:

	void state_transition(State new_state) override {
		if(new_state == triangle_mesh) {
			// if state "triangle mesh but no labeling", auto-compute labeling with graph cut optimisation & default parameters
			// Same code as in labeling.cpp graphcut_labeling(), copied here because we need the GraphCutLabeling object for update_per_chart_avg_data_cost()
			Attribute<index_t> label(mesh_.facets.attributes(), LABELING_ATTRIBUTE_NAME);
			auto gcl = GraphCutLabeling(mesh_);
			gcl.data_cost__set__fidelity_based(fidelity_coeff_);
			gcl.smooth_cost__set__default();
			gcl.neighbors__set__compactness_based(compactness_coeff_);
			gcl.compute_solution(label);
			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			update_per_chart_avg_data_cost(gcl);
			new_state = labeling;
			geo_assert(!SimpleApplication::colormaps_.empty()); // should I wait for GL to be initialized? see LabelingViewerApp::GL_initialize()
		}
		LabelingViewerApp::state_transition(new_state);
	}

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(state_ == LabelingViewerApp::State::labeling) {
			ImGui::SeparatorText("Global edition");

			ImGui::Text("Graph-cut parameters");

			ImGui::InputInt("Compactness", &compactness_coeff_);
			ImGui::InputInt("Fidelity", &fidelity_coeff_);

			ImGui::Text("Smooth cost");
			if (ImGui::BeginTable("Smooth cost table", 7, ImGuiTableFlags_Borders)) { // 7 columns
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::TextUnformatted(" ");
				for (int column = 0; column < 6; column++) {
					ImGui::TableSetColumnIndex(column+1);
					ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) column)); // change the text color according to the current label
					ImGui::TextUnformatted(LABEL2STR(column));
					ImGui::PopStyleColor();
				}
				for (int row = 0; row < 6; row++)
				{
					ImGui::TableNextRow();
					ImGui::TableSetColumnIndex(0);
					ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) row)); // change the text color according to the current label
					ImGui::TextUnformatted(LABEL2STR(row));
					ImGui::PopStyleColor();
					for (int column = 0; column < 6; column++)
					{
						ImGui::TableSetColumnIndex(column+1);
						char label[32];
						sprintf(label, "##smooth:%d,%d", column, row);
						ImGui::SetNextItemWidth(100.0f);
						ImGui::InputInt(label,&smooth_cost_[(std::array<int, 36>::size_type) (column*6+row)]);
					}
				}
				ImGui::EndTable();
			}

			if(ImGui::Button("Compute solution")) {
				auto gcl = GraphCutLabeling(mesh_); // TODO give normals, to not compute them each time
				gcl.data_cost__set__fidelity_based(fidelity_coeff_);
				gcl.smooth_cost__set__custom(smooth_cost_);
				gcl.neighbors__set__compactness_based(compactness_coeff_);
				Attribute<index_t> label(mesh_.facets.attributes(), LABELING_ATTRIBUTE_NAME);
				gcl.compute_solution(label);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				update_per_chart_avg_data_cost(gcl);
			}

			ImGui::SeparatorText("Per-chart edition");

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
			ImGui::TextUnformatted("Average data costs");
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 0)); // change the text color
			ImGui::SliderFloat("+X",&new_data_cost_[0],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 1)); // change the text color
			ImGui::SliderFloat("-X",&new_data_cost_[1],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 2)); // change the text color
			ImGui::SliderFloat("+Y",&new_data_cost_[2],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 3)); // change the text color
			ImGui::SliderFloat("-Y",&new_data_cost_[3],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 4)); // change the text color
			ImGui::SliderFloat("+Z",&new_data_cost_[4],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 5)); // change the text color
			ImGui::SliderFloat("-Z",&new_data_cost_[5],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::EndDisabled();
		}
	}

	void mouse_button_callback(int button, int action, int mods, int source) override {
		if(selected_chart_mode_) {
			selected_chart_mode_ = false;
			index_t facet_index = pick(MESH_FACETS);
			if( (facet_index == index_t(-1)) || (facet_index >= mesh_.facets.nb()) )
				selected_chart_ = index_t(-1);
			else {
				selected_chart_ = static_labeling_graph_.facet2chart[facet_index];
				geo_assert(selected_chart_ < per_chart_avg_data_cost_.size()); // assert per_chart_avg_data_cost_ is up to date
				new_data_cost_ = per_chart_avg_data_cost_[selected_chart_];
				geo_assert(min(new_data_cost_) >= 0.0f);
				new_data_cost_upper_bound_ = 2.0f * max(new_data_cost_);
			}
		}
		else {
			LabelingViewerApp::mouse_button_callback(button,action,mods,source);
		}
	}

	bool load(const std::string& filename) override {

		if(FileSystem::extension(filename)=="txt") { // if loading of a labeling
			fmt::println(Logger::err("I/O"),"Loading of a labeling is disabled in graphcut_labeling app"); Logger::err("I/O").flush();
        	return false;
		}
		return LabelingViewerApp::load(filename);
	}

	void update_per_chart_avg_data_cost(const GraphCutLabeling& gcl) {
		per_chart_avg_data_cost_.resize(static_labeling_graph_.nb_charts()); // each vec6f is initialized to 0 in constructor
		FOR(chart_index,static_labeling_graph_.nb_charts()) { // for each chart
			for(index_t f : static_labeling_graph_.charts[chart_index].facets) { // for each facet of the current chart
				per_chart_avg_data_cost_[chart_index] += (vec6f) gcl.data_cost__get__for_facet(f); // add per-label data cost of current facet
			}
			per_chart_avg_data_cost_[chart_index] /= (float) mesh_.facets.nb(); // divide by the number of facets to get the average
		}
	}

	int compactness_coeff_;
	int fidelity_coeff_;
	std::array<int,6*6> smooth_cost_;
	index_t selected_chart_;
	bool selected_chart_mode_;
	std::vector<vec6f> per_chart_avg_data_cost_;
	vec6f new_data_cost_; // per label
	float new_data_cost_upper_bound_; // max value of GUI sliders
};

int main(int argc, char** argv) {
    GraphCutLabelingApp app;
	app.start(argc,argv);
    return 0;
}