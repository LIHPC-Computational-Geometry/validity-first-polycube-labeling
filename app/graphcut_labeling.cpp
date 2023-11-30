// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/graphcut_labeling.cpp
// using Geogram instead of libigl, and with a GUI.

// TODO do not lauch GUI if both parameters are given by CLI

#include "LabelingViewerApp.h"
#include "labeling.h"

#include <array>
#include <algorithm>	// for std::max_element(), std::min_element()
#include <cmath>		// for std::round()

#include "GraphCutLabeling.h"

template <typename T>
struct StatsComponents {
	T min;
	T max;
	T avg;

	void reset() {
		// when T is a GEO::vecng<>, the constructor is a vector with zeros
		min = T();
		max = T();
		avg = T();
	}
};

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
		selected_chart_data_cost_stats_.reset();
		// new_data_cost_ auto-initialized to 0
		new_data_cost_upper_bound_ = 1000.0f;
	}

protected:

	void state_transition(State new_state) override {
		if(new_state == triangle_mesh) {
			// if state "triangle mesh but no labeling", auto-compute labeling with graph cut optimisation & default parameters
			if(SimpleApplication::colormaps_.empty()) {
				fmt::println(Logger::err("I/O"),"A mesh has been loaded but GL is still not initialized -> cannot auto-compute the labeling from graph-cut optimization (colormaps not defined). Will retry later."); Logger::err("I/O").flush();
				return;
			}
			data_cost_.resize(mesh_.facets.nb()*6);
			GraphCutLabeling::fill_data_cost__fidelity_based(mesh_,normals_,data_cost_,fidelity_coeff_);
			graphcut_labeling(mesh_,normals_,LABELING_ATTRIBUTE_NAME,compactness_coeff_,fidelity_coeff_); // compute graph-cut with init value of parameters
			update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			new_state = labeling;
		}
		LabelingViewerApp::state_transition(new_state);
	}

	void GL_initialize() override {
		LabelingViewerApp::GL_initialize();
		if(mesh_.vertices.nb()) { // if a mesh has been loaded
			state_transition(triangle_mesh);
		}
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
				auto gcl = GraphCutLabeling(mesh_,normals_);
				gcl.data_cost__set__fidelity_based(fidelity_coeff_);
				gcl.smooth_cost__set__custom(smooth_cost_);
				gcl.neighbors__set__compactness_based(compactness_coeff_);
				Attribute<index_t> label(mesh_.facets.attributes(), LABELING_ATTRIBUTE_NAME);
				gcl.compute_solution(label);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				selected_chart_ = index_t(-1); // avoid slider values to be obsolete (chart n°selected_chart_ may not be the same, or be out-of-range)
				new_data_cost_ = vec6f(); // zero values
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
			ImGui::SameLine();
			ImGui::Text("min:%6.2f max:%6.2f avg:%6.2f",(double) selected_chart_data_cost_stats_.min[0],
														(double) selected_chart_data_cost_stats_.max[0],
														(double) selected_chart_data_cost_stats_.avg[0]);
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 1)); // change the text color
			ImGui::SliderFloat("-X",&new_data_cost_[1],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::SameLine();
			ImGui::Text("min:%6.2f max:%6.2f avg:%6.2f",(double) selected_chart_data_cost_stats_.min[1],
														(double) selected_chart_data_cost_stats_.max[1],
														(double) selected_chart_data_cost_stats_.avg[1]);
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 2)); // change the text color
			ImGui::SliderFloat("+Y",&new_data_cost_[2],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::SameLine();
			ImGui::Text("min:%6.2f max:%6.2f avg:%6.2f",(double) selected_chart_data_cost_stats_.min[2],
														(double) selected_chart_data_cost_stats_.max[2],
														(double) selected_chart_data_cost_stats_.avg[2]);
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 3)); // change the text color
			ImGui::SliderFloat("-Y",&new_data_cost_[3],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::SameLine();
			ImGui::Text("min:%6.2f max:%6.2f avg:%6.2f",(double) selected_chart_data_cost_stats_.min[3],
														(double) selected_chart_data_cost_stats_.max[3],
														(double) selected_chart_data_cost_stats_.avg[3]);
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 4)); // change the text color
			ImGui::SliderFloat("+Z",&new_data_cost_[4],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::SameLine();
			ImGui::Text("min:%6.2f max:%6.2f avg:%6.2f",(double) selected_chart_data_cost_stats_.min[4],
														(double) selected_chart_data_cost_stats_.max[4],
														(double) selected_chart_data_cost_stats_.avg[4]);
			ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) 5)); // change the text color
			ImGui::SliderFloat("-Z",&new_data_cost_[5],0.0f,new_data_cost_upper_bound_);
			ImGui::PopStyleColor();
			ImGui::SameLine();
			ImGui::Text("min:%6.2f max:%6.2f avg:%6.2f",(double) selected_chart_data_cost_stats_.min[5],
														(double) selected_chart_data_cost_stats_.max[5],
														(double) selected_chart_data_cost_stats_.avg[5]);
			if(ImGui::Button("Apply")) {
				std::array<float,6> per_label_shift = {
					new_data_cost_[0] - selected_chart_data_cost_stats_.avg[0],
					new_data_cost_[1] - selected_chart_data_cost_stats_.avg[1],
					new_data_cost_[2] - selected_chart_data_cost_stats_.avg[2],
					new_data_cost_[3] - selected_chart_data_cost_stats_.avg[3],
					new_data_cost_[4] - selected_chart_data_cost_stats_.avg[4],
					new_data_cost_[5] - selected_chart_data_cost_stats_.avg[5]
				};
				auto gcl = GraphCutLabeling(mesh_,normals_);
				gcl.smooth_cost__set__custom(smooth_cost_);
				gcl.neighbors__set__compactness_based(compactness_coeff_);
				for(index_t f : static_labeling_graph_.charts[selected_chart_].facets) { // for each facet of the current chart
					FOR(l,6) {
						GraphCutLabeling::shift_data_cost(data_cost_,(GCoptimization::SiteID) f,l,per_label_shift[l]); // modify the local data cost vector
					}
				}
				gcl.data_cost__set__all_at_once(data_cost_); // use the local data cost vector in the optimization
				Attribute<index_t> label(mesh_.facets.attributes(), LABELING_ATTRIBUTE_NAME);
				gcl.compute_solution(label);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
				selected_chart_ = index_t(-1); // there is no way of finding the "new" chart, because it may have disappeared, moved, merged...
				selected_chart_data_cost_stats_.reset();
				new_data_cost_ = vec6f(); // zero values
			}
			ImGui::EndDisabled();
		}
	}

	void mouse_button_callback(int button, int action, int mods, int source) override {
		if(selected_chart_mode_) {
			selected_chart_mode_ = false;
			index_t facet_index = pick(MESH_FACETS);
			if( (facet_index == index_t(-1)) || (facet_index >= mesh_.facets.nb()) ) {
				selected_chart_ = index_t(-1);
				selected_chart_data_cost_stats_.reset();
			}
			else {
				// get the chart index of the picked facet
				selected_chart_ = static_labeling_graph_.facet2chart[facet_index];
				// compute data cost stats for the current chart
				vec6i current_facet_data_cost;
				for(index_t f : static_labeling_graph_.charts[selected_chart_].facets) { // for each facet of the selected chart
					current_facet_data_cost = GraphCutLabeling::per_site_data_cost_as_vector(data_cost_,(GCoptimization::SiteID) f);
					selected_chart_data_cost_stats_.avg += (vec6f) current_facet_data_cost; // add per-label data cost of current facet
					FOR(label,6) {
						selected_chart_data_cost_stats_.min[label] = std::min(selected_chart_data_cost_stats_.min[label], (float) current_facet_data_cost[label]);
						selected_chart_data_cost_stats_.max[label] = std::max(selected_chart_data_cost_stats_.max[label], (float) current_facet_data_cost[label]);
					}
				}
				selected_chart_data_cost_stats_.avg /= (float) mesh_.facets.nb(); // divide by the number of facets to get the average
				new_data_cost_ = selected_chart_data_cost_stats_.avg;
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

	virtual void update_static_labeling_graph(bool allow_boundaries_between_opposite_labels) {
		LabelingViewerApp::update_static_labeling_graph(allow_boundaries_between_opposite_labels);
		// update data cost sliders upper bound
		float global_max = 0.0f; // max of all data costs, for all facets and all labels
		FOR(chart_index,static_labeling_graph_.nb_charts()) { // for each chart
			for(index_t f : static_labeling_graph_.charts[chart_index].facets) { // for each facet of the current chart
				global_max = std::max(global_max,(float) max(GraphCutLabeling::per_site_data_cost_as_vector(data_cost_,(GCoptimization::SiteID) f)));
			}
		}
		new_data_cost_upper_bound_ = global_max*1.1f;
	}

	int compactness_coeff_;
	int fidelity_coeff_;
	std::vector<int> data_cost_;
	std::array<int,6*6> smooth_cost_;
	index_t selected_chart_;
	bool selected_chart_mode_;
	StatsComponents<vec6f> selected_chart_data_cost_stats_; // 3 vec6f : per-label min, per-label max and per-label avg
	vec6f new_data_cost_; // per label
	float new_data_cost_upper_bound_; // max value of GUI sliders
};

int main(int argc, char** argv) {
    GraphCutLabelingApp app;
	app.start(argc,argv);
    return 0;
}