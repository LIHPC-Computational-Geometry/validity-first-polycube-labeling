#include "gui_labeling.h"

LabelingViewerApp::LabelingViewerApp(const std::string name, bool auto_flip_normals) :
    SimpleMeshApplicationExt(name),
    labeling_colors_({
    {1.0f, 0.0f, 0.0f, 1.0f}, // label 0 -> red
    {0.6f, 0.0f, 0.0f, 1.0f}, // label 1 -> darker red
#ifdef RED_WHITE_BLUE_LABELING_COLORS
    {1.0f, 1.0f, 1.0f, 1.0f}, // label 2 -> white
    {0.6f, 0.6f, 0.6f, 1.0f}, // label 3 -> grey
#else
    {0.0f, 1.0f, 0.0f, 1.0f}, // label 2 -> green
    {0.0f, 0.6f, 0.0f, 1.0f}, // label 3 -> darker green
#endif
    {0.0f, 0.0f, 1.0f, 1.0f}, // label 4 -> blue
    {0.0f, 0.0f, 0.6f, 1.0f}  // label 5 -> darker blue
    }),
    validity_colors_({
    { 1.0f, 0.25f, 0.25f, 1.0f}, // invalid charts/boundaries/corners in light red
    {0.25f, 0.25f,  1.0f, 1.0f}  // valid charts/boundaries/corners in light blue
    }),
    auto_flip_normals_(auto_flip_normals),
    mesh_ext_(mesh_)
{
    // init some inherited variables

    show_vertices_ = false;
    show_surface_ = true;
    show_mesh_ = true;
    show_surface_borders_ = false;
    show_volume_ = false;
    lighting_ = false;
    surface_color_ = vec4f(0.9f, 0.9f, 0.9f, 1.0f); // light grey. default is blue-ish vec4f(0.5f, 0.5f, 1.0f, 1.0f)

    // init own variables

    show_ImGui_demo_window_ = false;
    show_normals_ = false;
    normals_color_ = vec3f(0.0f,1.0f,0.0f); // facet normals in green by default
    normals_length_factor_ = 0.1f;
    // facet_center_ and normal_tip_ init to 0.0,0.0,0.0 by default
    show_feature_edges_ = true;
    feature_edges_width_ = 5;
    allow_boundaries_between_opposite_labels_ = false; // parameter of LabelingGraph::fill_from()
    sensitivity_ = DEFAULT_SENSITIVITY;
	angle_of_rotation_ = DEFAULT_ANGLE_OF_ROTATION;
    fidelity_graphcuts_coeff_ = 3;
	compactness_graphcuts_coeff_ = 1;
    show_boundaries_ = false;
    boundaries_width_ = 6;
    X_boundaries_group_index_ = 0;
    Y_boundaries_group_index_ = 0;
    Z_boundaries_group_index_ = 0;
    axisless_and_invalid_boundaries_group_index = 0;
    axisless_but_valid_boundaries_group_index_ = 0;
    show_corners_ = false;
    corners_color_ = vec3f(0.0f,0.0f,0.0f); // corners in black by default
    corners_size_ = 10.0f;
    show_turning_points_ = false;
    turning_points_color_ = vec3f(1.0f,1.0f,0.0f); // turning-points in yellow by default
    turning_points_size_ = 10.0f;
    turning_points_group_index_ = 0;
    valid_corners_group_index_ = 0;
    invalid_corners_group_index_ = 0;
    fidelity_text_label_ = "";

    state_transition(empty);
}

void LabelingViewerApp::ImGui_initialize() {
    Application::ImGui_initialize();
    set_style("Dark");
    if(GEO::FileSystem::is_file("gui.ini")) {
        // Layout modification, saved with ImGui::SaveIniSettingsToDisk()
        // Larger docked object properties panel
        ImGui::LoadIniSettingsFromDisk("gui.ini");
    }
}

void LabelingViewerApp::state_transition(State new_state) {
    switch(new_state) {
        case empty:
            break;
        case triangle_mesh:
            show_surface_ = true;
            show_mesh_ = true;
            show_normals_ = false;
            show_feature_edges_ = true;
            lighting_ = true;
            show_attributes_ = false;
            show_boundaries_ = false;
            show_corners_ = false;
            show_turning_points_ = false;
            break;
        case labeling:
            mesh_ext_.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);
            labeling_visu_mode_transition(VIEW_LABELING_GRAPH);
            break;
        default:
            geo_assert_not_reached;
    }
    state_ = new_state;
}

void LabelingViewerApp::labeling_visu_mode_transition(int new_mode) {
    if(colormaps_.empty()) {
        // GL is not initialized yet
        // the state transition will be triggered later, in GL_initialize()
        return;
    }
    switch(new_mode) {
        case VIEW_TRIANGLE_MESH:
            show_surface_ = true;
            show_mesh_ = true;
            show_normals_ = false;
            show_feature_edges_ = true;
            lighting_ = true;
            show_attributes_ = false;
            show_boundaries_ = false;
            show_corners_ = false;
            show_turning_points_ = false;
            break;
        case VIEW_RAW_LABELING:
            show_surface_ = true;
            show_mesh_ = true;
            show_normals_ = false;
            show_feature_edges_ = false;
            lighting_ = false;
            show_attributes_ = true;
            current_colormap_index_ = COLORMAP_LABELING;
            attribute_ = fmt::format("facets.{}",LABELING_ATTRIBUTE_NAME);
            attribute_subelements_ = MESH_FACETS;
            attribute_name_ = LABELING_ATTRIBUTE_NAME;
            attribute_min_ = -0.5f;
            attribute_max_ = 5.5f;
            show_boundaries_ = false;
            show_corners_ = false;
            show_turning_points_ = false;
            break;
        case VIEW_LABELING_GRAPH:
            show_surface_ = true;
            show_mesh_ = false;
            show_normals_ = false;
            show_feature_edges_ = false;
            lighting_ = false;
            show_attributes_ = true;
            geo_assert(COLORMAP_LABELING < colormaps_.size());
            current_colormap_index_ = COLORMAP_LABELING;
            attribute_ = fmt::format("facets.{}",LABELING_ATTRIBUTE_NAME);
            attribute_subelements_ = MESH_FACETS;
            attribute_name_ = LABELING_ATTRIBUTE_NAME;
            attribute_min_ = -0.5f;
            attribute_max_ = 5.5f;
            // points in overlay
            set_points_group_color(valid_corners_group_index_,corners_color_.data());
            set_points_group_color(invalid_corners_group_index_,corners_color_.data()); // no distinction between valid and invalid corners in this view
            // edges in overlay
            show_boundaries_ = true;
            set_edges_group_color(X_boundaries_group_index_,COLORMAP_LABELING,0.084); // axis X -> color of label 0 = +X
            set_edges_group_color(Y_boundaries_group_index_,COLORMAP_LABELING,0.417); // axis Y -> color of label 2 = +Y
            set_edges_group_color(Z_boundaries_group_index_,COLORMAP_LABELING,0.750); // axis Z -> color of label 4 = +Z
            set_edges_group_color(axisless_and_invalid_boundaries_group_index,COLORMAP_BLACK_WHITE,0.0); // axisless boundaries in black
            set_edges_group_color(axisless_but_valid_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0); // axisless boundaries in black
            show_corners_ = true;
            show_turning_points_ = true;
            break;
        case VIEW_FIDELITY:
            show_surface_ = true;
            show_mesh_ = false;
            show_normals_ = false;
            show_feature_edges_ = false;
            lighting_ = false;
            show_attributes_ = true;
            current_colormap_index_ = COLORMAP_INFERNO;
            attribute_ = "facets.fidelity";
            attribute_subelements_ = MESH_FACETS;
            attribute_name_ = "fidelity";
            attribute_min_ = 0.0f; // the fidelity should not be in [0:0.5] (label too far from the normal), so setting 0.5 as the min of the colormap allows to focus the range of interest [0.5:1], but will display all values in [0:0.5] in black...
            attribute_max_ = 1.0f;
            show_boundaries_ = false;
            set_edges_group_color(X_boundaries_group_index_,COLORMAP_LABELING,0.084); // axis X -> color of label 0 = +X
            set_edges_group_color(Y_boundaries_group_index_,COLORMAP_LABELING,0.417); // axis Y -> color of label 2 = +Y
            set_edges_group_color(Z_boundaries_group_index_,COLORMAP_LABELING,0.750); // axis Z -> color of label 4 = +Z
            set_edges_group_color(axisless_and_invalid_boundaries_group_index,COLORMAP_BLACK_WHITE,0.0); // axisless boundaries in black
            set_edges_group_color(axisless_but_valid_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0); // axisless boundaries in black
            show_corners_ = false;
            show_turning_points_ = false;
            break;
        case VIEW_INVALID_CHARTS:
            show_surface_ = true;
            show_mesh_ = false;
            show_normals_ = false;
            show_feature_edges_ = false;
            lighting_ = false;
            show_attributes_ = true;
            current_colormap_index_ = COLORMAP_VALIDITY;
            attribute_ = "facets.on_invalid_chart";
            attribute_subelements_ = MESH_FACETS;
            attribute_name_ = "on_invalid_chart";
            attribute_min_ = 1.5;
            attribute_max_ = -0.5;
            show_boundaries_ = true;
            // all boundaries in black
            set_edges_group_color(X_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(Y_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(Z_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(axisless_and_invalid_boundaries_group_index,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(axisless_but_valid_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            show_corners_ = false;
            show_turning_points_ = false;
            break;
        case VIEW_INVALID_BOUNDARIES:
            show_surface_ = true;
            show_mesh_ = false;
            show_normals_ = false;
            show_feature_edges_ = false;
            lighting_ = false;
            show_attributes_ = false;
            show_boundaries_ = true;
            set_edges_group_color(X_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
            set_edges_group_color(Y_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
            set_edges_group_color(Z_boundaries_group_index_,COLORMAP_VALIDITY,1.0); // apply the color of valid LabelingGraph components
            set_edges_group_color(axisless_and_invalid_boundaries_group_index,COLORMAP_VALIDITY,0.0);
            set_edges_group_color(axisless_but_valid_boundaries_group_index_,COLORMAP_VALIDITY,1.0);
            show_corners_ = false;
            show_turning_points_ = false;
            break;
        case VIEW_INVALID_CORNERS:
            show_surface_ = true;
            show_mesh_ = false;
            show_normals_ = false;
            show_feature_edges_ = false;
            lighting_ = false;
            show_attributes_ = false;
            // points in overlay
            set_points_group_color(valid_corners_group_index_,validity_colors_.color_as_floats(1)); // apply the color of valid LabelingGraph components
            set_points_group_color(invalid_corners_group_index_,validity_colors_.color_as_floats(0)); // apply the color of invalid LabelingGraph components
            show_boundaries_ = true;
            // all boundaries in black
            set_edges_group_color(X_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(Y_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(Z_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(axisless_and_invalid_boundaries_group_index,COLORMAP_BLACK_WHITE,0.0);
            set_edges_group_color(axisless_but_valid_boundaries_group_index_,COLORMAP_BLACK_WHITE,0.0);
            show_corners_ = true;
            show_turning_points_ = false;
            break;
        default:
            geo_assert_not_reached;
    }
    labeling_visu_mode_ = new_mode;
}

void LabelingViewerApp::draw_scene() {
    SimpleMeshApplicationExt::draw_scene();
    if((state_ == triangle_mesh) && show_normals_) {
        glupSetColor4fv(GLUP_FRONT_COLOR, normals_color_.data());
        glupSetPointSize(10.0);
        FOR(f,mesh_.facets.nb()) { // for each 
            facet_center_ = mesh_facet_center(mesh_,f);
            normal_tip_ = facet_center_ + mesh_ext_.facet_normals[f] * normals_length_factor_;
            glupBegin(GLUP_LINES);
            glupPrivateVertex3dv(facet_center_.data());
            glupPrivateVertex3dv(normal_tip_.data());
            glupEnd();
            glupBegin(GLUP_POINTS);
            glupPrivateVertex3dv(normal_tip_.data());
            glupEnd();
        }
    }
    if(
        ((state_ == triangle_mesh) && show_feature_edges_) ||
        ((state_ == labeling) && (labeling_visu_mode_ == VIEW_TRIANGLE_MESH))
    ) {
        // dark blue color (coord. 0.0 of parula colormap)
        glupEnable(GLUP_TEXTURING);
        glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
        glBindTexture(GL_TEXTURE_2D, colormaps_[COLORMAP_PARULA].texture);
        glupTextureType(GLUP_TEXTURE_2D);
        glupTextureMode(GLUP_TEXTURE_REPLACE);
        glupPrivateTexCoord1d(0.0);
        glupSetMeshWidth(feature_edges_width_);
        glupBegin(GLUP_LINES);
        for(const std::pair<index_t,index_t>& edge : mesh_ext_.feature_edges) { // for each edge in the set of feature edges
            glupPrivateVertex3dv(mesh_.vertices.point_ptr(edge.first)); // draw first vertex
            glupPrivateVertex3dv(mesh_.vertices.point_ptr(edge.second)); // draw second vertex
        }
        glupDisable(GLUP_TEXTURING);
        glupEnd();
    }
}

void LabelingViewerApp::draw_gui() {
    SimpleMeshApplicationExt::draw_gui();
    if(show_ImGui_demo_window_)
        ImGui::ShowDemoWindow();
}

void LabelingViewerApp::draw_menu_bar() {
    SimpleApplication::draw_menu_bar();

    if(ImGui::BeginMainMenuBar()) {

        if(ImGui::MenuItem("Toggle light/dark mode")) {
            // ignoring 'DarkGray' and 'LightGray' styles
            if(get_style() == "Light") {
                set_style("Dark");
            }
            else {
                set_style("Light");
            }
        }

        if(ImGui::BeginMenu("Debug")) {
            if(ImGui::MenuItem("Dump labeling graph as text file")) {
                lg_.dump_to_text_file("LabelingGraph.txt",mesh_);
                fmt::println(Logger::out("I/O"),"Exported to LabelingGraph.txt"); Logger::out("I/O").flush();
            }
            if(ImGui::MenuItem("Dump labeling graph as D3 graph")) {
                lg_.dump_to_D3_graph("LabelingGraph.json");
                fmt::println(Logger::out("I/O"),"Exported to LabelingGraph.json"); Logger::out("I/O").flush();
            }
            if(ImGui::MenuItem("Dump boundaries as mesh")) {
                MeshHalfedgesExt mesh_he(mesh_);
                dump_all_boundaries_with_indices_and_axes("boundaries",mesh_,lg_);
            }
            if (ImGui::MenuItem("Show ImGui demo window", NULL, show_ImGui_demo_window_)) {
                show_ImGui_demo_window_ = !show_ImGui_demo_window_;
            }
            if(state_ == labeling) {
                if (ImGui::MenuItem("Flip labeling")) {
                    flip_labeling(mesh_,labeling_);
                    update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
                }
            }
            if (ImGui::MenuItem("Random labeling")) {
                random_labeling(mesh_,labeling_);
                mesh_ext_.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);
                update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
            }
            ImGui::EndMenu();
        }

        ImGui::EndMainMenuBar();
    }
}

void LabelingViewerApp::draw_object_properties() {

    if(ImGui::Button(icon_UTF8("sliders-h") + "graphics settings")) {
        ImGui::OpenPopup("Graphics settings");
    }
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f)); // Always center the modal when appearing
    if(ImGui::BeginPopupModal("Graphics settings", NULL, ImGuiWindowFlags_AlwaysAutoResize))
    {
        if (
            (state_ == empty) || 
            (state_ == triangle_mesh) || 
            ((state_ == labeling) && (labeling_visu_mode_ == VIEW_TRIANGLE_MESH))
        ) {
            // Surface
            ImGui::Checkbox("##Show surface",&show_surface_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Surface color", surface_color_.data());
            ImGui::SameLine();
            ImGui::TextUnformatted("Surface");
            // Mesh wireframe
            ImGui::Checkbox("##Show mesh wireframe",&show_mesh_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Mesh wireframe color", mesh_color_.data());
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Mesh wireframe width", &mesh_width_, 0.1f, 2.0f, "%.1f");
            ImGui::SameLine();
            ImGui::TextUnformatted("Mesh wireframe");
            // Facet normals
            ImGui::Checkbox("##Show normals",&show_normals_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Normals color", normals_color_.data());
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Normals length factor",&normals_length_factor_,0.0f,5.0f,"%.1f");
            ImGui::SameLine();
            ImGui::TextUnformatted("Facet normals");
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Draw facet normals as segments. The point is the tip of the vector.");
            // Feature edges
            ImGui::Checkbox("##Show feature edges",&show_feature_edges_);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderInt("##Feature edges width",&feature_edges_width_,1,30);
            ImGui::SameLine();
            ImGui::Text("Feature edges (nb=%d)",mesh_ext_.feature_edges.nb());
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Draw feature edges in blue on top of the mesh");
            // Feature edges color (thick lines color in general) cannot be easily changed, because we have to use a colormap coordinate, not just a float[3]
        }
        else if (
            ((state_ == labeling) && (labeling_visu_mode_ == VIEW_RAW_LABELING)) || 
            ((state_ == labeling) && (labeling_visu_mode_ == VIEW_LABELING_GRAPH))
        ) {
            // Mesh wireframe
            ImGui::Checkbox("##Show mesh wireframe",&show_mesh_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Mesh wireframe color", mesh_color_.data());
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Mesh wireframe width", &mesh_width_, 0.1f, 2.0f, "%.1f");
            ImGui::SameLine();
            ImGui::TextUnformatted("Mesh wireframe");
            // Charts & per label color
            ImGui::Checkbox("Show charts",&show_surface_);
            ImGui::SameLine();
            ImGui::Text("(nb=%ld)",lg_.nb_charts());
            ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(),ImGui::GetFrameHeight()));
            ImGui::SameLine();
            if(ImGui::ColorEdit3WithPalette("Label 0 = +X", labeling_colors_.color_as_floats(0))) {
                labeling_colors_.update_chars_of_color(0);
                update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
            }
            ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(),ImGui::GetFrameHeight()));
            ImGui::SameLine();
            if(ImGui::ColorEdit3WithPalette("Label 1 = -X", labeling_colors_.color_as_floats(1))) {
                labeling_colors_.update_chars_of_color(1);
                update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
            }
            ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(),ImGui::GetFrameHeight()));
            ImGui::SameLine();
            if(ImGui::ColorEdit3WithPalette("Label 2 = +Y", labeling_colors_.color_as_floats(2))) {
                labeling_colors_.update_chars_of_color(2);
                update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
            }
            ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(),ImGui::GetFrameHeight()));
            ImGui::SameLine();
            if(ImGui::ColorEdit3WithPalette("Label 3 = -Y", labeling_colors_.color_as_floats(3))) {
                labeling_colors_.update_chars_of_color(3);
                update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
            }
            ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(),ImGui::GetFrameHeight()));
            ImGui::SameLine();
            if(ImGui::ColorEdit3WithPalette("Label 4 = +Z", labeling_colors_.color_as_floats(4))) {
                labeling_colors_.update_chars_of_color(4);
                update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
            }
            ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(),ImGui::GetFrameHeight()));
            ImGui::SameLine();
            if(ImGui::ColorEdit3WithPalette("Label 5 = -Z", labeling_colors_.color_as_floats(5))) {
                labeling_colors_.update_chars_of_color(5);
                update_GL_texture(COLORMAP_LABELING,6,1,labeling_colors_.as_chars());
            }
            // Boundaries
            ImGui::Checkbox("##Show boundaries",&show_boundaries_);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderInt("##Boundaries width",&boundaries_width_,1,30);
            ImGui::SameLine();
            ImGui::Text("Boundaries (nb=%ld)",lg_.nb_boundaries());
            // Corners
            ImGui::Checkbox("##Show corners",&show_corners_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Corners color", corners_color_.data());
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Corners size", &corners_size_, 1.0f, 50.0f, "%.1f");
            ImGui::SameLine();
            ImGui::Text("Corners (nb=%ld)",lg_.nb_corners());
            // Turning-points
            ImGui::Checkbox("##Show turning-points",&show_turning_points_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Turning-points color", turning_points_color_.data());
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Turning points size", &turning_points_size_, 1.0f, 50.0f, "%.1f");
            ImGui::SameLine();
            ImGui::Text("Turning-points (nb=%ld)",lg_.nb_turning_points());
        }
        else if (
            (state_ == labeling) && (labeling_visu_mode_ == VIEW_FIDELITY)
        ) {
            ImGui::Checkbox("##Show mesh wireframe",&show_mesh_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Mesh wireframe color", mesh_color_.data());
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Mesh wireframe width", &mesh_width_, 0.1f, 2.0f, "%.1f");
            ImGui::SameLine();
            ImGui::TextUnformatted("Mesh wireframe");
            // TODO allow to change colormap?
        }
        else if (
            ((state_ == labeling) && (labeling_visu_mode_ == VIEW_INVALID_CHARTS)) ||
            ((state_ == labeling) && (labeling_visu_mode_ == VIEW_INVALID_BOUNDARIES)) ||
            ((state_ == labeling) && (labeling_visu_mode_ == VIEW_INVALID_CORNERS))
        ) {
            // Surface
            ImGui::Checkbox("##Show surface",&show_surface_);
            ImGui::SameLine();
            ImGui::ColorEdit3WithPalette("##Surface color", surface_color_.data());
            ImGui::SameLine();
            ImGui::TextUnformatted("Surface");
            // Boundaries
            ImGui::Checkbox("##Show boundaries",&show_boundaries_);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderInt("##Boundaries width",&boundaries_width_,1,30);
            ImGui::SameLine();
            ImGui::Text("Boundaries (nb=%ld)",lg_.nb_boundaries());
            // Corners
            ImGui::Checkbox("##Show corners",&show_corners_);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(IMGUI_SLIDERS_WIDTH);
            ImGui::SliderFloat("##Corners size", &corners_size_, 1.0f, 50.0f, "%.1f");
            ImGui::SameLine();
            ImGui::Text("Corners (nb=%ld)",lg_.nb_corners());
            // color for invalid components
            if(ImGui::ColorEdit3WithPalette("Invalid", validity_colors_.color_as_floats(0))) {
                validity_colors_.update_chars_of_color(0);
                update_GL_texture(COLORMAP_VALIDITY,2,1,validity_colors_.as_chars());
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Color of invalid charts/boundaries/corners");
            // color for valid components
            if(ImGui::ColorEdit3WithPalette("Valid", validity_colors_.color_as_floats(1))) {
                validity_colors_.update_chars_of_color(1);
                update_GL_texture(COLORMAP_VALIDITY,2,1,validity_colors_.as_chars());
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Color of valid charts/boundaries/corners");
        }
        

        if(ImGui::Button("Close", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    ImGui::Separator();

    switch (state_) {

    case triangle_mesh:

        ImGui::TextUnformatted("Mesh transformation");

        if(ImGui::Button("Rotate mesh according to principal axes")) {
            rotate_mesh_according_to_principal_axes(mesh_);
            mesh_gfx_.set_mesh(&mesh_); // re-link the MeshGfx to the mesh
            mesh_ext_.facet_normals.recompute();
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Compute principal axes of the point cloud and rotate the mesh to be aligned with them");

        ImGui::Separator();
        ImGui::TextUnformatted("Labeling generation");

        ImGui::Checkbox("Allow boundaries between opposite labels",&allow_boundaries_between_opposite_labels_);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("If on, boundaries between opposite labels (e.g. +X and -X)\ncan be considered valid if they only contain > 180° angles");

        ImGui::TextUnformatted("Normals pre-processing:");
        ImGui::InputDouble("Sensitivity", &sensitivity_, 0.0, 0.0, "%.15f");
        ImGui::InputDouble("Angle of rotation", &angle_of_rotation_);

        if(ImGui::Button("Compute naive labeling")) {
            naive_labeling(mesh_ext_,labeling_,sensitivity_,angle_of_rotation_);
            mesh_ext_.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);
            update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
            state_transition(labeling);
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Associate each facet to the label the closest to its normal.\nIf sensitivity != 0.0, pre-process normals.");

        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputInt("fidelity coeff",&fidelity_graphcuts_coeff_);
        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputInt("compactness coeff",&compactness_graphcuts_coeff_);
        if(ImGui::Button("Compute graph-cuts labeling")) {
            graphcut_labeling(mesh_ext_,labeling_,compactness_graphcuts_coeff_,fidelity_graphcuts_coeff_,sensitivity_,angle_of_rotation_);
            mesh_ext_.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);
            update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
            state_transition(labeling);
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Generate a labeling with a graph-cuts optimization.\nIf sensitivity != 0.0, pre-process normals.");

        break;
    case labeling:

        if(ImGui::Button("Remove labeling")) {
            ImGui::OpenPopup("Remove labeling ?");
        }
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f)); // Always center the modal when appearing
        if(ImGui::BeginPopupModal("Remove labeling ?", NULL, ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::TextUnformatted("If not manually saved before, the labeling will be lost.");
            ImGui::TextUnformatted("Are you sure you want to remove the labeling?");
            ImGui::PushStyleColor(ImGuiCol_Button, 			ImVec4(0.9f, 0.4f, 0.4f, 1.0f)); // red
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered,	ImVec4(0.9f, 0.2f, 0.2f, 1.0f)); // darker red
            ImGui::PushStyleColor(ImGuiCol_ButtonActive,	ImVec4(0.9f, 0.0f, 0.0f, 1.0f));
            if(ImGui::Button("Yes, remove the labeling", ImVec2(200, 0))) {
                clear_scene_overlay();
                state_transition(triangle_mesh);
                ImGui::CloseCurrentPopup();
            }
            ImGui::PopStyleColor(3);
            ImGui::SetItemDefaultFocus();
            ImGui::SameLine();
            if(ImGui::Button("No, cancel", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }

        ImGui::Checkbox("Allow boundaries between opposite labels",&allow_boundaries_between_opposite_labels_);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("If on, boundaries between opposite labels (e.g. +X and -X)\ncan be considered valid if they only contain > 180° angles");

        ImGui::BeginDisabled( allow_boundaries_between_opposite_labels_ == lg_.is_allowing_boundaries_between_opposite_labels() ); // allow to recompute only if the UI control value changed
        if(ImGui::Button("Recompute labeling graph")) {
            update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Update charts, boundaries and corners according to the new value of \"Allow boundaries between opposite labels\"");
        ImGui::EndDisabled();

        ImGui::Separator();

        if(ImGui::RadioButton("View triangle mesh",&labeling_visu_mode_,VIEW_TRIANGLE_MESH))
            labeling_visu_mode_transition(VIEW_TRIANGLE_MESH);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Show the triangle mesh without the labeling");
        ImGui::Text("%d vertices, %d facets, %d feature edges",mesh_.vertices.nb(),mesh_.facets.nb(),mesh_ext_.feature_edges.nb());

        if(ImGui::RadioButton("View raw labeling",&labeling_visu_mode_,VIEW_RAW_LABELING))
            labeling_visu_mode_transition(VIEW_RAW_LABELING);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Show the triangle mesh and the labeling");

        if(ImGui::RadioButton("View labeling graph",&labeling_visu_mode_,VIEW_LABELING_GRAPH))
            labeling_visu_mode_transition(VIEW_LABELING_GRAPH);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Show the charts, the boundaries, the corners and the turning-points computed from the labeling");

        ImGui::Text("%ld charts, %ld boundaries, %ld corners, %ld turning-points",lg_.nb_charts(),lg_.nb_boundaries(),lg_.nb_corners(),lg_.nb_turning_points());

        if(ImGui::RadioButton("View fidelity",&labeling_visu_mode_,VIEW_FIDELITY))
            labeling_visu_mode_transition(VIEW_FIDELITY);
        ImGui::SameLine();
        ImGui::TextDisabled("(?)");
        ImGui::SetItemTooltip("Color the facets according to the angle between the normal and the assigned direction.\nYellow = small angle, black = wide angle.\nClick on a facet to print its fidelity.");
        
        ImGui::TextUnformatted(fidelity_text_label_.c_str());

        if(ImGui::RadioButton("View invalid charts",&labeling_visu_mode_,VIEW_INVALID_CHARTS))
            labeling_visu_mode_transition(VIEW_INVALID_CHARTS);
        ImGui::SameLine();
        ImGui::Text("(nb=%ld)",lg_.nb_invalid_charts());

        if(ImGui::RadioButton("View invalid boundaries",&labeling_visu_mode_,VIEW_INVALID_BOUNDARIES))
            labeling_visu_mode_transition(VIEW_INVALID_BOUNDARIES);
        ImGui::SameLine();
        ImGui::Text("(nb=%ld)",lg_.nb_invalid_boundaries());

        if(ImGui::RadioButton("View invalid corners",&labeling_visu_mode_,VIEW_INVALID_CORNERS))
            labeling_visu_mode_transition(VIEW_INVALID_CORNERS);
        ImGui::SameLine();
        ImGui::Text("(nb=%ld)",lg_.nb_invalid_corners());

        if(lg_.is_valid()) {
            ImGui::TextColored(ImVec4(0.0f,0.5f,0.0f,1.0f),"Valid labeling");
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("The current labeling is a valid polycube representation");
        }
        else {
            ImGui::TextColored(ImVec4(0.8f,0.0f,0.0f,1.0f),"Invalid labeling");
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("The current labeling is not valid polycube representation.\nValence or adjacency of some components (charts, boundaries, corners) cannot turn into polycube components.");
        }

        if(lg_.nb_turning_points()==0) {
            ImGui::TextColored(ImVec4(0.0f,0.5f,0.0f,1.0f),"All monotone boundaries");
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("There are no turning-points");
        }
        else {
            ImGui::TextColored(ImVec4(0.8f,0.0f,0.0f,1.0f),"Non-monotone boundaries");
            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Some boundaries contain turning-points");
        }
        
        break;
    
    default:
        break;
    }
}

bool LabelingViewerApp::load(const std::string& filename) {

    //// based on SimpleMeshApplication::load() ////////////////////////////////////////////////////////////////

    if(!FileSystem::is_file(filename)) {
        Logger::out("I/O") << "is not a file" << std::endl;
    }

    // added : load a labeling
    if(FileSystem::extension(filename)=="txt") {

        if(state_ == empty) {
            fmt::println(Logger::err("I/O"),"You need to import a triangle mesh before importing a labeling"); Logger::err("I/O").flush();
            return false;
        }

        if(!load_labeling(filename,mesh_,labeling_)) {
            fmt::println("load_labeling() not ok"); fflush(stdout);
            // Should the labeling be removed ?
            // If a labeling was already displayed, it should be restored...
            mesh_.facets.clear(false,true);
            clear_scene_overlay();
            state_transition(triangle_mesh);
            return false;
        }

        mesh_ext_.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);
        update_static_labeling_graph(allow_boundaries_between_opposite_labels_);
        state_transition(labeling);
        return true;
    }

    mesh_gfx_.set_mesh(nullptr);
    mesh_ext_.clear();
    if(labeling_.is_bound()) {
        labeling_.unbind();
        // labeling_.unregister_me(mesh_.facets.attributes().find_attribute_store(LABELING_ATTRIBUTE_NAME));
    }
    mesh_.clear(false,true);
    MeshIOFlags flags;
    if(!mesh_load(filename, mesh_, flags)) {
        state_transition(empty); // added
        mesh_ext_.facet_normals.clear(); // added
        mesh_ext_.adj_facet_corners.clear(); // added
        mesh_ext_.feature_edges.clear(); // added
        if(labeling_.is_bound()) {
            labeling_.unbind();
            // labeling_.unregister_me(mesh_.facets.attributes().find_attribute_store(LABELING_ATTRIBUTE_NAME));
        }
        return false;
    }
    mesh_gfx_.set_animate(false);
    mesh_.vertices.set_dimension(3);
    double xyzmin[3];
    double xyzmax[3];
    get_bbox(mesh_, xyzmin, xyzmax, false);
    set_region_of_interest(
        xyzmin[0], xyzmin[1], xyzmin[2],
        xyzmax[0], xyzmax[1], xyzmax[2]
    );
    mesh_.vertices.set_double_precision(); // added
    mesh_gfx_.set_mesh(&mesh_);
    current_file_ = filename;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(auto_flip_normals_) {
        // ensure facet normals are outward
        if(facet_normals_are_inward(mesh_)) {
            flip_facet_normals(mesh_);
            fmt::println(Logger::warn("normals dir."),"Facet normals of the input mesh were inward");
            fmt::println(Logger::warn("normals dir."),"You should flip them and update the file for consistency.");
            Logger::warn("normals dir.").flush();
        }
    }

    if(labeling_.is_bound()) {
        labeling_.unbind();
        // labeling_.unregister_me(mesh_.facets.attributes().find_attribute_store(LABELING_ATTRIBUTE_NAME));
    }
    labeling_.bind(mesh_.facets.attributes(),LABELING_ATTRIBUTE_NAME);
    mesh_ext_.facet_normals.recompute();
    mesh_ext_.adj_facet_corners.recompute();
    mesh_ext_.feature_edges.recompute();

    clear_scene_overlay();
    state_transition(triangle_mesh);

    return true;
}

std::string LabelingViewerApp::supported_write_file_extensions() {
    return SimpleMeshApplication::supported_write_file_extensions() + ";txt"; // add .txt in supported write file extensions
}

bool LabelingViewerApp::save(const std::string& filename) {
    if(String::string_ends_with(filename,".txt")) { // bypass inherited save behavior in case of a .txt file -> save the labeling only
        save_labeling(filename,mesh_,labeling_);
        fmt::println(Logger::out("I/O"),"Labeling saved to {}",filename); Logger::out("I/O").flush();
        return true;
    }
    else {
        return SimpleMeshApplication::save(filename);
    }
}

void LabelingViewerApp::GL_initialize() {
    SimpleMeshApplicationExt::GL_initialize();
    init_rgba_colormap("labeling",6,1,labeling_colors_.as_chars());
    init_rgba_colormap("validity",2,1,validity_colors_.as_chars());
    state_transition(state_); // not all state_transition() code has been executed if GL was not initialized (in particular because missing colormaps)
}

void LabelingViewerApp::mouse_button_callback(int button, int action, int mods, int source) {
    if((action==EVENT_ACTION_DOWN) && (button == 0) ) { // if left click
        if( (state_ == labeling) && ((labeling_visu_mode_ == VIEW_RAW_LABELING) ||
                                     (labeling_visu_mode_ == VIEW_FIDELITY)) ) { // if current mode is "view raw labeling" or "view fidelity"
            index_t facet_index = pick(MESH_FACETS);
            if ( (facet_index != index_t(-1)) && (facet_index < mesh_.facets.nb()) ) {
                Attribute<index_t> label(mesh_.facets.attributes(), LABELING_ATTRIBUTE_NAME);
                Attribute<double> per_facet_fidelity(mesh_.facets.attributes(), "fidelity");
                vec3 label_direction = label2vector[label[facet_index]];
                fmt::println(Logger::out("fidelity"),"facet #{} : normal=({:.4f},{:.4f},{:.4f}), label={}=({:.4f},{:.4f},{:.4f}) -> fidelity={:.4f}",
                    facet_index,
                    mesh_ext_.facet_normals[facet_index].x,
                    mesh_ext_.facet_normals[facet_index].y,
                    mesh_ext_.facet_normals[facet_index].z,
                    LABEL2STR(label[facet_index]),
                    label_direction.x,
                    label_direction.y,
                    label_direction.z,
                    per_facet_fidelity[facet_index]
                );
                Logger::out("fidelity").flush();
            }
        }
    }
    SimpleMeshApplication::mouse_button_callback(button,action,mods,source);
}

void LabelingViewerApp::update_static_labeling_graph(bool allow_boundaries_between_opposite_labels) {

    // compute charts, boundaries and corners of the labeling
    mesh_ext_.halfedges.set_use_facet_region(LABELING_ATTRIBUTE_NAME);
    lg_.fill_from(mesh_ext_,labeling_,allow_boundaries_between_opposite_labels);

    clear_scene_overlay();

    valid_corners_group_index_ = new_points_group(corners_color_.data(),&corners_size_,&show_corners_);
    invalid_corners_group_index_ = new_points_group(corners_color_.data(),&corners_size_,&show_corners_);
    turning_points_group_index_ = new_points_group(turning_points_color_.data(),&turning_points_size_,&show_turning_points_);

    X_boundaries_group_index_ = new_edges_group(COLORMAP_LABELING,0.084,&boundaries_width_,&show_boundaries_); // axis X -> use the color of label 0 = +X
    Y_boundaries_group_index_ = new_edges_group(COLORMAP_LABELING,0.417,&boundaries_width_,&show_boundaries_); // axis Y -> use the color of label 2 = +Y
    Z_boundaries_group_index_ = new_edges_group(COLORMAP_LABELING,0.750,&boundaries_width_,&show_boundaries_); // axis Z -> use the color of label 4 = +Z
    axisless_and_invalid_boundaries_group_index = new_edges_group(COLORMAP_VALIDITY,0.0,&boundaries_width_,&show_boundaries_); // use the color of invalid LabelingGraph components
    axisless_but_valid_boundaries_group_index_ = new_edges_group(COLORMAP_VALIDITY,1.0,&boundaries_width_,&show_boundaries_); // use the color of valid LabelingGraph components

    for(std::size_t i = 0; i < lg_.nb_corners(); ++i) {
        const double* coordinates = mesh_.vertices.point_ptr(
            lg_.corners[i].vertex
        );
        add_point_to_group(lg_.corners[i].is_valid ? valid_corners_group_index_ : invalid_corners_group_index_,coordinates[0], coordinates[1], coordinates[2]);
    }

    for(index_t i : lg_.non_monotone_boundaries) {
        FOR(j,lg_.boundaries[i].turning_points.size()) {
            vec3 coordinates = mesh_vertex(mesh_,lg_.boundaries[i].turning_point_vertex(j,mesh_));
            add_point_to_group(turning_points_group_index_,coordinates.x,coordinates.y,coordinates.z);
        }
    }

    std::size_t group_index;
    for(std::size_t i = 0; i < lg_.nb_boundaries(); ++i) {
        const Boundary& boundary = lg_.boundaries[i];
        for(const auto& be : boundary.halfedges) { // for each boundary edge of this boundary
            const double* coordinates_first_point = halfedge_vertex_from(mesh_,be).data();
            const double* coordinates_second_point = halfedge_vertex_to(mesh_,be).data();
            switch(boundary.axis) {
                case -1: // may or may not be invalid
                    if(boundary.is_valid) {
                        group_index = axisless_but_valid_boundaries_group_index_;
                    }
                    else {
                        group_index = axisless_and_invalid_boundaries_group_index;
                    }
                    break;
                case 0: // always valid
                    group_index = X_boundaries_group_index_;
                    break;
                case 1: // always valid
                    group_index = Y_boundaries_group_index_;
                    break;
                case 2: // always valid
                    group_index = Z_boundaries_group_index_;
                    break;
                default:
                    geo_assert_not_reached;
            }
            add_edge_to_group(
                group_index,
                coordinates_first_point[0], coordinates_first_point[1], coordinates_first_point[2], // first point
                coordinates_second_point[0], coordinates_second_point[1], coordinates_second_point[2] // second point
            );
        }
    }

    IncrementalStats stats;
    compute_per_facet_fidelity(mesh_ext_,labeling_,"fidelity",stats);
    fidelity_text_label_ = fmt::format("min={:.4f} | max={:.4f} | avg={:.4f}",stats.min(),stats.max(),stats.avg());
}