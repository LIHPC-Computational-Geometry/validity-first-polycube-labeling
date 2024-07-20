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
#include <geogram/basic/command_line.h>         // for declare_arg(), parse(), get_arg_bool()
#include <geogram/basic/command_line_args.h>    // for import_arg_group()
#include <geogram/mesh/mesh_io.h>               // for mesh_load()
#include <geogram/basic/numeric.h>              // for min_float64()

#include <set>
#include <array>
#include <algorithm>            // for std::min(), std::max(), std::min_element()
#include <tuple>                // for std::tie()
#include <deque>

#include "stats.h" // for std_dev()
#include "gui_labeling.h"
#include "labeling_graph.h"
#include "labeling_io.h"
#include "labeling_generators.h"
#include "labeling_operators_on_invalidity.h"
#include "labeling_operators_on_distortion.h"

#define STDDEV_BASED_CHART_REFINEMENT

class AutomaticPolycubeApp : public LabelingViewerApp {
public:

    AutomaticPolycubeApp() : LabelingViewerApp("automatic_polycube"),
							 selected_chart_(index_t(-1)),
							 selected_chart_mode_(false)
							 {}

protected:

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(state_ == LabelingViewerApp::State::labeling) {

			ImGui::SeparatorText("Fix labeling");

			if(ImGui::Button("Remove surrounded charts")) {
				size_t nb_invalid_charts_processed = remove_surrounded_charts(labeling_,lg_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid charts processed",nb_invalid_charts_processed); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Increase chart valence")) {
				bool chart_processed = increase_chart_valence(mesh_ext_,labeling_,lg_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid charts processed",size_t(chart_processed)); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Pick the first invalid chart surrounded by feature edges\nand trace a new adjacent chart to increase the valence of the first one");

			if(ImGui::Button("Fix invalid boundaries (as much as possible)")) {
				size_t nb_invalid_boundaries_processed = (size_t) fix_an_invalid_boundary(mesh_ext_,labeling_,lg_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid boundaries processed",nb_invalid_boundaries_processed); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Fix invalid corners (as much as possible)")) {
				size_t nb_invalid_corners_processed = fix_as_much_invalid_corners_as_possible(mesh_ext_,labeling_,lg_);
				fmt::println(Logger::out("fix_labeling"),"{} invalid corners processed",nb_invalid_corners_processed); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Remove invalid charts")) {
				remove_invalid_charts(mesh_ext_,labeling_,lg_);
				fmt::println(Logger::out("fix_labeling"),"invalid charts removed using gco"); Logger::out("fix_labeling").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Remove charts around invalid boundaries")) {
				remove_charts_around_invalid_boundaries(mesh_ext_,labeling_,lg_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			ImGui::PushStyleColor(ImGuiCol_Button, 			ImVec4(0.6f, 0.9f, 0.6f, 1.0f)); // green
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered,	ImVec4(0.3f, 0.8f, 0.3f, 1.0f)); // darker green
			ImGui::PushStyleColor(ImGuiCol_ButtonActive,	ImVec4(0.0f, 0.8f, 0.0f, 1.0f));
			if(ImGui::Button("Auto fix validity")) {
				auto_fix_validity(mesh_ext_,labeling_,lg_,100);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::PopStyleColor(3);

			ImGui::SeparatorText("Monotonicity");

			if(ImGui::Button("Join turning-points pair with new chart")) {
				join_turning_points_pair_with_new_chart(mesh_,labeling_,lg_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Merge a turning-points and its closest corners")) {
				bool a_non_monotone_boundary_was_processed = merge_a_turning_point_and_its_closest_corner(mesh_,labeling_,lg_);
				fmt::println(Logger::out("monotonicity"),"{} non-monotone boundaries processed",(size_t) a_non_monotone_boundary_was_processed); Logger::out("monotonicity").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Parse non-monotone boundaries, and if one of them has a turning-point on a feature edge, try to pull the closest corner so that they coincide");

			if(ImGui::Button("Move boundaries near turning points")) {
				size_t nb_turning_points_processed = move_boundaries_near_turning_points(mesh_,labeling_,lg_);
				fmt::println(Logger::out("monotonicity"),"{} turning-points processed",nb_turning_points_processed); Logger::out("monotonicity").flush();
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			if(ImGui::Button("Straighten boundaries")) {
				straighten_boundaries(mesh_,labeling_,lg_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::SameLine();
			ImGui::TextDisabled("(?)");
			ImGui::SetItemTooltip("Re-draw boundaries so that they are straighter");

			if(ImGui::Button("Move corners")) {
				move_corners(mesh_,labeling_,lg_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}

			ImGui::PushStyleColor(ImGuiCol_Button, 			ImVec4(0.6f, 0.9f, 0.6f, 1.0f)); // green
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered,	ImVec4(0.3f, 0.8f, 0.3f, 1.0f)); // darker green
			ImGui::PushStyleColor(ImGuiCol_ButtonActive,	ImVec4(0.0f, 0.8f, 0.0f, 1.0f));
			if(ImGui::Button("Auto fix monotonicity")) {
				auto_fix_monotonicity(mesh_,labeling_,lg_);
				update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
			}
			ImGui::PopStyleColor(3);
		}
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

	bool allow_boundaries_between_opposite_labels = true;
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
	// Compute initial labeling
	//////////////////////////////////////////////////

	Attribute<index_t> labeling(M.facets.attributes(),LABELING_ATTRIBUTE_NAME);
	smart_init_labeling(M,labeling);

	//////////////////////////////////////////////////
	// Construct charts, boundaries & corners
	//////////////////////////////////////////////////

	LabelingGraph lg;
	lg.fill_from(M,labeling,allow_boundaries_between_opposite_labels);

	//////////////////////////////////////////////////
	// Validity & monotonicity correction
	//////////////////////////////////////////////////

	if(auto_fix_validity(M,labeling,lg,100)) {
		// auto-fix the monotonicity only if the validity was fixed
		auto_fix_monotonicity(M,labeling,lg);
	}

	//////////////////////////////////////////////////
	// Write output labeling
	//////////////////////////////////////////////////

	fmt::println(Logger::out("I/O"),"Writing {}...",filenames[1]); Logger::out("I/O").flush();
	save_labeling(filenames[1],M,labeling);
	fmt::println(Logger::out("I/O"),"Done"); Logger::out("I/O").flush();
	
    return 0;
}