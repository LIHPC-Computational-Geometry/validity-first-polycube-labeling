// Rewriting of https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/apps/graphcut_labeling.cpp
// using Geogram instead of libigl, and with a GUI.

// TODO do not lauch GUI if both parameters are given by CLI

#include <geogram/basic/command_line_args.h> // for import_arg_group()

#include <algorithm>	// for std::max_element(), std::min_element()
#include <cmath>		// for std::round()

#include "labeling_graphcuts.h" // for DEFAULT_COMPACTNESS, DEFAULT_FIDELITY, DEFAULT_SENSITIVITY, DEFAULT_ANGLE_OF_ROTATION
#include "labeling_generators.h"
#include "gui_labeling.h"
#include "labeling.h"
#include "containers_Geogram.h"	// for max() on a Geogram vector
#include "geometry.h" // for rotation_matrix()

#define DEFAULT_OUTPUT_FILENAME "labeling.txt"

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

    GraphCutLabelingApp(int compactness = DEFAULT_COMPACTNESS, int fidelity = DEFAULT_FIDELITY) : LabelingViewerApp("graphcut_labeling") {
		compactness_coeff_ = compactness;
		fidelity_coeff_ = fidelity;
		smooth_cost_.resize(6*6);
		// fill smooth_cost_. equivalent to what GraphCutLabeling::smooth_cost__set__default() does.
		FOR(label1,6) {
        	FOR(label2,6) {
				// same label = very smooth edge, different label = less smooth
				smooth_cost_[label1+label2*6] = (label1==label2) ? 0 : 1;
			}
		}
		sensitivity_ = DEFAULT_SENSITIVITY;
		angle_of_rotation_ = DEFAULT_ANGLE_OF_ROTATION;
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
			GraphCutLabeling::fill_data_cost__fidelity_based(mesh_ext_.facet_normals.as_vector(),data_cost_,fidelity_coeff_,mesh_);
			graphcut_labeling(mesh_ext_,labeling_,compactness_coeff_,fidelity_coeff_); // compute graph-cut with init value of parameters
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

			ImGui::TextUnformatted("Normals pre-processing:");
			ImGui::InputDouble("Sensitivity", &sensitivity_, 0.0, 0.0, "%.15f");
			ImGui::InputDouble("Angle of rotation", &angle_of_rotation_);

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
						ImGui::InputInt(label,&smooth_cost_[(std::vector<int>::size_type) (column*6+row)]);
					}
				}
				ImGui::EndTable();
			}

			if(ImGui::Button("Compute solution")) {
				auto gcl = GraphCutLabeling(mesh_ext_);
				// TODO remove redundancy with src/labeling_generators.cpp graphcut_labeling()
				if(std::fpclassify(sensitivity_) == FP_ZERO) {
					// sensitivity is null, use the fidelity based data cost
					gcl.data_cost__set__fidelity_based(fidelity_coeff_);
				}
				else {
					std::vector<int> custom_data_cost(mesh_ext_.facets.nb()*6);
					vec3 normal;
					mat3 rotation_to_apply = rotation_matrix(angle_of_rotation_);
					FOR(f,mesh_ext_.facets.nb()) {
						normal = mesh_ext_.facet_normals[f];
						if(is_a_facet_to_tilt(normal,sensitivity_)) {
							normal = mult(rotation_to_apply,normal); // rotation of the normal
						}
						FOR(label,6) {
							double dot = (GEO::dot(normal,label2vector[label]) - 1.0)/0.2;
							double cost = 1.0 - std::exp(-(1.0/2.0)*std::pow(dot,2));
							custom_data_cost[f*6+label] = (int) (fidelity_coeff_*100*cost);
						}
					}
					gcl.data_cost__set__all_at_once(custom_data_cost);
				}
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
				auto gcl = GraphCutLabeling(mesh_ext_);
				gcl.smooth_cost__set__custom(smooth_cost_);
				gcl.neighbors__set__compactness_based(compactness_coeff_);
				for(index_t f : lg_.charts[selected_chart_].facets) { // for each facet of the current chart
					FOR(l,6) {
						GraphCutLabeling::shift_data_cost(data_cost_,(GCoptimization::SiteID) f,(GCoptimization::LabelID) l,per_label_shift[l]); // modify the local data cost vector
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
				selected_chart_ = lg_.facet2chart[facet_index];
				// compute data cost stats for the current chart
				vec6i current_facet_data_cost;
				for(index_t f : lg_.charts[selected_chart_].facets) { // for each facet of the selected chart
					current_facet_data_cost = GraphCutLabeling::per_siteID_data_cost_as_vector(data_cost_,(GCoptimization::SiteID) f,6,mesh_.facets.nb());
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
		FOR(chart_index,lg_.nb_charts()) { // for each chart
			for(index_t f : lg_.charts[chart_index].facets) { // for each facet of the current chart
				global_max = std::max(global_max,(float) max(GraphCutLabeling::per_siteID_data_cost_as_vector(data_cost_,(GCoptimization::SiteID) f,6,mesh_.facets.nb())));
			}
		}
		new_data_cost_upper_bound_ = global_max*1.1f;
	}

	int compactness_coeff_;
	int fidelity_coeff_;
	std::vector<int> data_cost_;
	std::vector<int> smooth_cost_;
	double sensitivity_;
	double angle_of_rotation_;
	index_t selected_chart_;
	bool selected_chart_mode_;
	StatsComponents<vec6f> selected_chart_data_cost_stats_; // 3 vec6f : per-label min, per-label max and per-label avg
	vec6f new_data_cost_; // per label
	float new_data_cost_upper_bound_; // max value of GUI sliders
};

int main(int argc, char** argv) {

	GEO::initialize();

	CmdLine::import_arg_group("standard");
	CmdLine::declare_arg(
		"gui",
		true,
		"Show the graphical user interface"
	);
	CmdLine::declare_arg(
        "output",
        "",
        "where to write the output labeling (in case gui=false)"
    ); // not a positional arg, to have a similar CLI with/without GUI
	CmdLine::declare_arg(
		"compactness",
		DEFAULT_COMPACTNESS,
		"the compactness coefficient"
	);
	CmdLine::declare_arg(
		"fidelity",
		DEFAULT_FIDELITY,
		"the fidelity coefficient"
	);

	std::vector<std::string> filenames;
	if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"<input_surface_mesh>"
		))
	{
		return 1;
	}
	
	if(CmdLine::get_arg_bool("gui")) { // if GUI mode
		GraphCutLabelingApp app(
			CmdLine::get_arg_int("compactness"),
			CmdLine::get_arg_int("fidelity")
		);
		app.start(argc,argv);
		return 0;
	}
    // else: no GUI, auto-process the input mesh

	if(filenames.size() < 1) { // missing filenames[0], that is <input_surface_mesh> argument
		fmt::println(Logger::err("I/O"),"When gui=false, the input filename must be given as CLI argument"); Logger::err("I/O").flush();
		return 1;
	}

	std::string output_labeling_path = GEO::CmdLine::get_arg("output");
	if(output_labeling_path.empty()) {
		fmt::println(Logger::warn("I/O"),"The output filename was not provided, using default '{}'",DEFAULT_OUTPUT_FILENAME); Logger::warn("I/O").flush();
		output_labeling_path = DEFAULT_OUTPUT_FILENAME;
	}

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

	MeshExt M_ext(M); // will compute facet normals, feature edges, and vertex -> facets adjacency

	//////////////////////////////////////////////////
	// Labeling optimization with Graph-cuts
	//////////////////////////////////////////////////

	fmt::println(Logger::out("graph-cut"),"Using {} as compactness coeff",CmdLine::get_arg_int("compactness"));
	fmt::println(Logger::out("graph-cut"),"Using {} as fidelity coeff",CmdLine::get_arg_int("fidelity"));

	Attribute<index_t> labeling(M_ext.facets.attributes(),LABELING_ATTRIBUTE_NAME);
	graphcut_labeling(
		M_ext,
		labeling,
		CmdLine::get_arg_int("compactness"),
		CmdLine::get_arg_int("fidelity")
	);

	//////////////////////////////////////////////////
	// Write output labeling
	//////////////////////////////////////////////////

	fmt::println(Logger::out("I/O"),"Writing {}...",output_labeling_path); Logger::out("I/O").flush();
	save_labeling(output_labeling_path,M_ext,labeling);
	fmt::println(Logger::out("I/O"),"Done"); Logger::out("I/O").flush();
}