#include <geogram_gfx/gui/simple_application.h>
#include <geogram_gfx/gui/geogram_logo_256.xpm>
#include <geogram/basic/file_system.h>
#include <geogram/basic/string.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/bibliography/bibliography.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/stopwatch.h>

#include <geogram_gfx/full_screen_effects/ambient_occlusion.h>
#include <geogram_gfx/full_screen_effects/unsharp_masking.h>

#ifdef GEOGRAM_WITH_LUA

#  include <geogram_gfx/lua/lua_glup.h>
#  include <geogram_gfx/lua/lua_simple_application.h>
#  include <geogram_gfx/lua/lua_imgui.h>
#  include <geogram/lua/lua_io.h>

extern "C" {
#    include <geogram/third_party/lua/lua.h>    
#    include <geogram/third_party/lua/lauxlib.h>
#    include <geogram/third_party/lua/lualib.h>
}

#endif

#ifdef GEO_OS_EMSCRIPTEN
#  include <emscripten.h>
#endif

#ifdef GEO_OS_ANDROID
#  include <geogram/basic/android_utils.h>
#endif

namespace {
#  include <geogram_gfx/gui/gui_state_v.h>
#  include <geogram_gfx/gui/gui_state_h.h>
#  include <geogram_gfx/gui/gui_state.h>
}

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <fstream>
#include <stdexcept>

#include "LabelingViewerApp.h"
#include "colormaps_array.h" // for colormap_name, colormap_xpm & macros
#include "CustomMeshGfx.h"   // for CustomMeshGfx
#include "LabelingGraph.h"   // for StaticLabelingGraph
#include "labeling.h"		 // for load_labeling(), naive_labeling()

/******************************************************************************/

namespace {

   /**
    * \brief Converts a complete path to a file to a label
    *  displayed in the file browser.
    * \details Strips viewer_path from the input path.
    * \param[in] path the complete path, can be either a directory or
    *  a file
    * \return the label to be displayed in the menu
    */
    std::string path_to_label(
        const std::string& viewer_path, const std::string& path
    ) {
        std::string result = path;
        if(GEO::String::string_starts_with(result, viewer_path)) {
            result = result.substr(
                viewer_path.length(), result.length()-viewer_path.length()
            );
        }
        return result;
    }

    
}

/******************************************************************************/

    LabelingViewerApp::LabelingViewerApp() :
	Application("labeling_viewer"),
	text_editor_(&text_editor_visible_),
	previous_state_(empty),
	current_state_(empty)
    {
	lighting_ = false;
	edit_light_ = false;
	clipping_ = false;
	clip_mode_ = GLUP_CLIP_STRADDLING_CELLS;
	edit_clip_ = false;
	fixed_clip_ = false;
	background_color_ = vec4f(0.0f, 0.0f, 0.0f, 1.0f);	
	effect_ = GLenum(0);
	
	filename_[0] = '\0';
	geogram_logo_texture_ = 0;
	viewer_properties_visible_ = true;
	object_properties_visible_ = true;
	console_visible_ = true;
	use_text_editor_           = false;
	text_editor_visible_       = false;
	menubar_visible_           = true;

        console_ = new Console(&console_visible_);
	console_->hide_command_prompt();
	text_editor_.set_fixed_layout(false);
        status_bar_ = new StatusBar;

	add_key_func("q", [this]() { stop(); }, "quit");
	add_key_func("z", [this]() { zoom_in(); }, "zoom in");
	add_key_func("Z", [this]() { zoom_out(); }, "zoom out");
	add_key_func("H", [this]() { home(); }, "home");		
	add_key_toggle("L",   &lighting_, "light");
	add_key_toggle("l",   &edit_light_, "light edit");
	add_key_toggle("F1",  &clipping_, "clipping");
	add_key_toggle("F2",  &edit_clip_, "clip plane edit");
	add_key_toggle("F3",  &fixed_clip_, "fixed clip plane");
	add_key_toggle("a",   animate_ptr(), "animate");
	add_key_toggle("F6",  &text_editor_visible_, "text editor");
	add_key_toggle("F7",  &viewer_properties_visible_, "viewer properties");
	add_key_toggle("F8",  &object_properties_visible_, "object properties");
	add_key_toggle("F9",  &console_visible_, "console");
	add_key_toggle("F12", &menubar_visible_, "menubar");
	set_region_of_interest(
	    0.0, 0.0, 0.0, 1.0, 1.0, 1.0
	);

	object_translation_ = vec3(0.0, 0.0, 0.0);

	mouse_op_ = MOUSE_NOOP;
	mouse_target_ = MOUSE_NOTARGET;

	three_D_ = true;
	zoom_ = 1.0;
	zoom_down_ = 1.0;

#ifdef GEOGRAM_WITH_LUA	
	lua_error_occured_ = false;
	lua_state_ = luaL_newstate();
	luaL_openlibs(lua_state_);
	init_lua_io(lua_state_);
	init_lua_glup(lua_state_);
	init_lua_simple_application(lua_state_);		
	init_lua_imgui(lua_state_);
#else
	lua_error_occured_ = false;
	lua_state_ = nullptr;
#endif
	geo_cite_with_info(
	    "WEB:ImGUI",
	    "Used to create the GUI of GEOGRAM utilities "
	    "(vorpaview, geobox, geocod)."
	);

	props_pinned_ = false;

    //copied from the constructor of SimpleMeshApplication in ext/geogram/src/lib/geogram_gfx/gui/simple_mesh_application.cpp
    set_default_filename("out.meshb");
	
        anim_speed_ = 1.0f;
        anim_time_ = 0.0f;

        show_vertices_ = false;
        show_vertices_selection_ = true;
        vertices_size_ = 1.0f;
	vertices_color_ = vec4f(0.0f, 1.0f, 0.0f, 1.0f);
	vertices_transparency_ = 0.0f;
	
        show_surface_ = true;
        show_surface_sides_ = false;
        show_mesh_ = false;
	mesh_color_ = vec4f(0.0f, 0.0f, 0.0f, 1.0f);
	mesh_width_ = 0.1f;
	
        show_surface_borders_ = false;
	surface_color_ =   vec4f(0.9f, 0.9f, 0.9f, 1.0f); // light grey. default is bleu-ish vec4f(0.5f, 0.5f, 1.0f, 1.0f)
	surface_color_2_ = vec4f(1.0f, 0.5f, 0.0f, 1.0f); 
	
        show_volume_ = false;
	volume_color_ = vec4f(0.9f, 0.9f, 0.9f, 1.0f);	
        cells_shrink_ = 0.0f;
        show_colored_cells_ = false;
        show_hexes_ = true;
	show_connectors_ = true;
	
        show_attributes_ = false;
        current_colormap_texture_ = 0;
        attribute_min_ = 0.0f;
        attribute_max_ = 0.0f;
        attribute_ = "vertices.point_fp32[0]";
        attribute_name_ = "point_fp32[0]";
        attribute_subelements_ = MESH_VERTICES;

        add_key_toggle("p", &show_vertices_, "vertices");
        add_key_toggle("S", &show_surface_, "surface");
        add_key_toggle("c", &show_surface_sides_, "two-sided");
        add_key_toggle("B", &show_surface_borders_, "borders");
        add_key_toggle("m", &show_mesh_, "mesh");
        add_key_toggle("V", &show_volume_, "volume");
        add_key_toggle("j", &show_hexes_, "hexes");
        add_key_toggle("k", &show_connectors_, "connectors");
        add_key_toggle("C", &show_colored_cells_, "colored cells");

		// added

		// default colors
		// label 0 -> red
		labeling_colors_[0][0] = 1.0f;
		labeling_colors_[0][1] = 0.0f;
		labeling_colors_[0][2] = 0.0f;
		labeling_colors_[0][3] = 1.0f;
		// label 1 -> dark red
		labeling_colors_[1][0] = 0.6f;
		labeling_colors_[1][1] = 0.0f;
		labeling_colors_[1][2] = 0.0f;
		labeling_colors_[1][3] = 1.0f;
		// label 2 -> green
		labeling_colors_[2][0] = 0.0f;
		labeling_colors_[2][1] = 1.0f;
		labeling_colors_[2][2] = 0.0f;
		labeling_colors_[2][3] = 1.0f;
		// label 3 -> dark green
		labeling_colors_[3][0] = 0.0f;
		labeling_colors_[3][1] = 0.6f;
		labeling_colors_[3][2] = 0.0f;
		labeling_colors_[3][3] = 1.0f;
		// label 4 -> blue
		labeling_colors_[4][0] = 0.0f;
		labeling_colors_[4][1] = 0.0f;
		labeling_colors_[4][2] = 1.0f;
		labeling_colors_[4][3] = 1.0f;
		// label 5 -> dark blue
		labeling_colors_[5][0] = 0.0f;
		labeling_colors_[5][1] = 0.0f;
		labeling_colors_[5][2] = 0.6f;
		labeling_colors_[5][3] = 1.0f;

		// give pointers to these colors to mesh_gfx_
		mesh_gfx_.bind_int_attribute_value_to_color(0,labeling_colors_[0]);
		mesh_gfx_.bind_int_attribute_value_to_color(1,labeling_colors_[1]);
		mesh_gfx_.bind_int_attribute_value_to_color(2,labeling_colors_[2]);
		mesh_gfx_.bind_int_attribute_value_to_color(3,labeling_colors_[3]);
		mesh_gfx_.bind_int_attribute_value_to_color(4,labeling_colors_[4]);
		mesh_gfx_.bind_int_attribute_value_to_color(5,labeling_colors_[5]);

		// corners in black by default
		corners_color_[0] = 0.0f;
		corners_color_[1] = 0.0f;
		corners_color_[2] = 0.0f;
		corners_color_[3] = 1.0f;

		// valid boundaries in blue by default
		valid_boundaries_color_[0] = 0.0f;
		valid_boundaries_color_[1] = 0.0f;
		valid_boundaries_color_[2] = 1.0f;
		valid_boundaries_color_[3] = 1.0f;

		// invalid boundaries in red by default
		invalid_boundaries_color_[0] = 1.0f;
		invalid_boundaries_color_[1] = 0.0f;
		invalid_boundaries_color_[2] = 0.0f;
		invalid_boundaries_color_[3] = 1.0f;

		previous_labeling_visu_mode_ = VIEW_RAW_LABELING;
		current_labeling_visu_mode_ = VIEW_RAW_LABELING;
    }

    LabelingViewerApp::~LabelingViewerApp() {
#ifdef GEOGRAM_WITH_LUA
	if(lua_state_ != nullptr) {
	    lua_close(lua_state_);
	    lua_state_ = nullptr;
	}
#endif	
    }
    
    void LabelingViewerApp::home() {
	zoom_ = 1.0;
	object_translation_ = vec3(0.0, 0.0, 0.0);
	object_rotation_.reset();
	light_rotation_.reset();
	clip_rotation_.reset();
	clip_translation_ = vec3(0.0, 0.0, 0.0);
	clipping_ = false;
	edit_clip_ = false;
	fixed_clip_ = false;
	edit_light_ = false;
    }

    void LabelingViewerApp::add_key_func(
	const std::string& key, std::function<void()> cb,
	const char* help
    ) {
	key_funcs_[key] = cb;
	if(help != nullptr) {
	    key_funcs_help_[key] = help;
	}
    }
    
    void LabelingViewerApp::add_key_toggle(
	const std::string& key, bool* p_val,
	const char* help
    ) {
	add_key_func(
	    key, [p_val]() { *p_val = !*p_val; },
	    (help != nullptr) ?
	    (std::string("toggle ") + help).c_str() : nullptr
	);
    }

    void LabelingViewerApp::char_callback(unsigned int c){
	Application::char_callback(c);
	if(text_editor_visible_) {
	    return;
	}
	std::string k=" ";
	k[0] = char(c);
	auto F = key_funcs_.find(k);
	if(F != key_funcs_.end()) {
	    F->second();
	}
    }

    void LabelingViewerApp::key_callback(
	int key, int scancode, int action, int mods
    ) {
	Application::key_callback(key, scancode, action, mods);
	if(action == 0) {
	    auto F = key_funcs_.find(std::string(key_to_string(key)));
	    if(F != key_funcs_.end()) {
		F->second();
	    }
	}
    }
    
    void LabelingViewerApp::set_style(const std::string& style) {
	Application::set_style(style);
	if(String::string_starts_with(style, "Light")) {
	    background_color_ = vec4f(1.0f, 1.0f, 1.0f, 1.0f);
	} else {
	    background_color_ = vec4f(0.0f, 0.0f, 0.0f, 1.0f);	    
	}
    }

    void LabelingViewerApp::set_region_of_interest(
	double xmin, double ymin, double zmin,
	double xmax, double ymax, double zmax
    ) {
	roi_.xyz_min[0] = xmin;
	roi_.xyz_min[1] = ymin;
	roi_.xyz_min[2] = zmin;
	roi_.xyz_max[0] = xmax;
	roi_.xyz_max[1] = ymax;
	roi_.xyz_max[2] = zmax;
	roi_radius_ = sqrt(
	    0.25 * (xmax - xmin) * (xmax - xmin) +
	    0.25 * (ymax - ymin) * (ymax - ymin) +
	    0.25 * (zmax - zmin) * (zmax - zmin)
	);
    }

    void LabelingViewerApp::get_region_of_interest(
	double& xmin, double& ymin, double& zmin,
	double& xmax, double& ymax, double& zmax
    ) const {
	xmin = roi_.xyz_min[0];
	ymin = roi_.xyz_min[1];
	zmin = roi_.xyz_min[2];
	xmax = roi_.xyz_max[0];
	ymax = roi_.xyz_max[1];
	zmax = roi_.xyz_max[2];
    }
    
    void LabelingViewerApp::draw_gui() {
#ifdef GEO_OS_ANDROID
	if(
	    supported_read_file_extensions() != "" ||
	    supported_write_file_extensions() != ""
	) {
	    static bool firsttime = true;
	    static int nb_perms = 2;
	    static const char* perms[] = {
		"READ_EXTERNAL_STORAGE",
		"WRITE_EXTERNAL_STORAGE"
	    };
	    if(firsttime) {
		firsttime = false;
		bool OK = true;
		for(int i=0; i<nb_perms; ++i) {
		    OK = OK && AndroidUtils::has_permission(
			CmdLine::get_android_app(), perms[i]
		    );
		}
		if(!OK) {
		    AndroidUtils::request_permissions(
			CmdLine::get_android_app(), nb_perms, perms
		    );
		}
	    }
	}
#endif	
	draw_menu_bar();
	draw_dock_space();
	draw_viewer_properties_window();
	draw_object_properties_window();
	draw_console();
	draw_command_window();
	if(text_editor_visible_) {
	    text_editor_.draw();
	}
	if(
	    ImGui::FileDialog("##load_dlg", filename_, geo_imgui_string_length)
	) {
	    load(filename_);
	}
	if(
	    ImGui::FileDialog("##save_dlg", filename_, geo_imgui_string_length)
	) {
	    save(filename_);
	}
	if(status_bar_->active()) {
	    float w = float(get_frame_buffer_width());
	    float h = float(get_frame_buffer_height());
	    float STATUS_HEIGHT = status_bar_->get_window_height();
	    if(STATUS_HEIGHT == 0.0f) {
		STATUS_HEIGHT = float(get_font_size());
		if(phone_screen_) {
		    STATUS_HEIGHT *= std::max(w,h)/600.f;
		}
	    }
	    STATUS_HEIGHT *= 1.5f;
	    ImGui::SetNextWindowPos(
		ImVec2(0.0f, h-STATUS_HEIGHT),
		ImGuiCond_Always
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(w,STATUS_HEIGHT-1.0f),
		ImGuiCond_Always
	    );
	    status_bar_->draw();
	}
    }

    void LabelingViewerApp::draw_scene_begin() {

	glClearColor(
	    background_color_.x,
	    background_color_.y,
	    background_color_.z,
	    background_color_.w
	);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glViewport(
	    0, 0,
	    int(double(get_frame_buffer_width())),
	    int(double(get_frame_buffer_height()))
	);

	double zScreen = 5.0; // screen projection plane 
	{
	    glupMatrixMode(GLUP_PROJECTION_MATRIX);
	    glupLoadIdentity();
	    
	    double aspect = double(get_width()) / double(get_height());
	    double zNear = 1.0;   // near clipping plane.
	    double zFar = 10.0;   // far clipping plane.
	    if(three_D_) {
		double camera_aperture = 9.0; // field of view in degrees.
		double view_max_size = zScreen * tan(
		    (camera_aperture * M_PI / 180.0) / 2.0
		);
		double right;
		double top;
		if(aspect < 1.0) {
		    top = view_max_size;
		    right = top * aspect;
		} else {
		    right = view_max_size;
		    top = right / aspect;
		}
		right /= zoom_;
		top   /= zoom_;
		glupFrustum(-right, right, -top, top, zNear, zFar);
	    } else {
		double x = 1.0 / zoom_;
		double y = 1.0 / zoom_;
		if(aspect > 1.0) {
		    x *= aspect;
		} else {
		    y /= aspect;
		}
		glupOrtho(-x, x, -y, y, zNear, zFar);
	    }
	}

	static vec3 light0 = normalize(vec3(1.0, 1.0, 1.0));
	vec3 light = transform_vector(
	    light0, light_rotation_.get_value()
	);
	glupLightVector3f(float(light.x), float(light.y), float(light.z));
    
	glupMatrixMode(GLUP_MODELVIEW_MATRIX);
	glupLoadIdentity();
	glupTranslated(0.0, 0.0, -zScreen);

	// Clipping
	{
	    static double clip_eqn[4];
	    clip_eqn[0] = 0.0;
	    clip_eqn[1] = 0.0;
	    clip_eqn[2] = 0.0;
	    clip_eqn[3] = 0.0;
	    glupClipPlane(clip_eqn);
	    glupDisable(GLUP_CLIPPING);
	    
	    if(clipping_) {
		glupPushMatrix();
		glupTranslated(
		    clip_translation_.x,
		    clip_translation_.y,
		    clip_translation_.z
		);
		glupMultMatrix(clip_rotation_.get_value());
		
		{
		    if(effect_ != 0) {
			glDepthMask(GL_FALSE);
		    }
		    glupSetColor3f(
			GLUP_FRONT_AND_BACK_COLOR,
			1.0f - background_color_.x,
			1.0f - background_color_.y,
			1.0f - background_color_.z
		    );
		    float sq_w = 1.25f / float(zoom_);

		    GLboolean vertex_colors_save =
			glupIsEnabled(GLUP_VERTEX_COLORS);
		    GLboolean texturing_save = glupIsEnabled(GLUP_TEXTURING);
		    glupDisable(GLUP_VERTEX_COLORS);
		    glupDisable(GLUP_TEXTURING);

		    // Draw the cross 
		    glupBegin(GLUP_LINES);
		    glupVertex3f(-sq_w, 0.0f, 0.0f);
		    glupVertex3f(sq_w, 0.0f, 0.0f);
		    glupVertex3f(0.0f, -sq_w, 0.0f);
		    glupVertex3f(0.0f, sq_w, 0.0f);

		    // Draw the square around the cross 
		    for(index_t i=0; i<3; ++i) {
			glupVertex3f(sq_w, -sq_w, 0.0f);
			glupVertex3f(sq_w, sq_w, 0.0f);
			glupVertex3f(sq_w, sq_w, 0.0f);            
			glupVertex3f(-sq_w, sq_w, 0.0f);
			glupVertex3f(-sq_w, sq_w, 0.0f);            
			glupVertex3f(-sq_w, -sq_w, 0.0f);
			glupVertex3f(-sq_w, -sq_w, 0.0f);
			glupVertex3f(sq_w, -sq_w, 0.0f);
			sq_w = sq_w * 1.01f;
		    }
		    glupEnd();
		    if(vertex_colors_save) {
			glupEnable(GLUP_VERTEX_COLORS);
		    }
		    if(texturing_save) {
			glupEnable(GLUP_TEXTURING);
		    }
		}
		clip_eqn[0] = 0.0;
		clip_eqn[1] = 0.0;
		clip_eqn[2] = -1.0;
		clip_eqn[3] = 0.0;
		glupEnable(GLUP_CLIPPING); 
		glupClipPlane(clip_eqn);
		glupClipMode(clip_mode_);
		glupPopMatrix();
		if(effect_ != 0) {
		    glDepthMask(GL_TRUE);
		}
	    } 
	}
	
	glupTranslate(object_translation_);
	glupMultMatrix(object_rotation_.get_value());
	
	glupScaled(
	    1.5 / roi_radius_, 1.5 / roi_radius_, 1.5 / roi_radius_
	);
	glupTranslated(
	    -0.5 * (roi_.xyz_min[0] + roi_.xyz_max[0]),
	    -0.5 * (roi_.xyz_min[1] + roi_.xyz_max[1]),
	    -0.5 * (roi_.xyz_min[2] + roi_.xyz_max[2])
	);

	if(lighting_) {
	    glupEnable(GLUP_LIGHTING);
	} else {
	    glupDisable(GLUP_LIGHTING);	    
	}

	// Save transform for picking
	{
	    glGetIntegerv(GL_VIEWPORT, viewport_);
	    // Note: OpenGL uses column-major order for matrices
	    // (thus what we get is the transpose of each matrix)
	    glupGetMatrixdv(GLUP_MODELVIEW_MATRIX, modelview_transpose_.data());
	    glupGetMatrixdv(GLUP_PROJECTION_MATRIX, project_transpose_.data());
	}
    }

    vec3 LabelingViewerApp::project(const vec3& p) {
	vec3 result;
	glupProject(
	    p.x, p.y, p.z,
	    modelview_transpose_.data(), project_transpose_.data(), viewport_,
	    &result.x, &result.y, &result.z
	);
	return result;
    }

    vec3 LabelingViewerApp::unproject(const vec3& p) {
	vec3 result;
	glupUnProject(
	    p.x, p.y, p.z,
	    modelview_transpose_.data(), project_transpose_.data(), viewport_,
	    &result.x, &result.y, &result.z
	);
	return result;
    }

    vec2 LabelingViewerApp::unproject_2d(const vec2& p) {
	double z = project(vec3(0.0, 0.0, 0.0)).z;
	vec3 result3d = unproject(vec3(p.x, p.y, z));
	return vec2(result3d.x, result3d.y);
    }
    
    void LabelingViewerApp::draw_scene_end() {
    }

    void LabelingViewerApp::draw_scene() {
        // copied from ext/geogram/src/lib/geogram_gfx/gui/simple_mesh_application.cpp
        if(mesh_gfx_.mesh() == nullptr) {
            return;
        }
        
        if(animate()) {
            anim_time_ = float(
                sin(double(anim_speed_) * GEO::SystemStopwatch::now())
            );
            anim_time_ = 0.5f * (anim_time_ + 1.0f);
        }
        
        mesh_gfx_.set_lighting(lighting_);
        mesh_gfx_.set_time(double(anim_time_));

		// manage state transitions
		if(previous_state_ != current_state_) {
			switch(current_state_) {
				case empty:
					break;
				case triangle_mesh:
					show_mesh_ = true;
					lighting_ = true;
					mesh_gfx_.unset_scalar_attribute();
					mesh_gfx_.unset_facets_colors_by_int_attribute();
					break;
				case labeling:
					previous_labeling_visu_mode_ = VIEW_TRIANGLE_MESH;
					current_labeling_visu_mode_ = VIEW_LABELING_GRAPH;
					// will trigger a GUI settings update because !=
					break;
				default:
					geo_assert_not_reached;
			}
			previous_state_ = current_state_;
		}

		// manage labeling visu mode transitions (mode editable when current_state_==labeling)
		if(previous_labeling_visu_mode_ != current_labeling_visu_mode_) {
			switch(current_labeling_visu_mode_) {
				case VIEW_TRIANGLE_MESH:
					show_mesh_ = true;
					lighting_ = true;
					mesh_gfx_.unset_scalar_attribute();
					mesh_gfx_.unset_facets_colors_by_int_attribute();
					mesh_gfx_.set_custom_points_group_visibility(corner_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(X_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(Y_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(Z_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(invalid_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(valid_but_axisless_boundaries_group_index,false);
					break;
				case VIEW_RAW_LABELING:
					show_mesh_ = true;
					lighting_ = false;
					mesh_gfx_.unset_scalar_attribute();
					mesh_gfx_.set_facets_colors_by_int_attribute(LABELING_ATTRIBUTE_NAME);
					mesh_gfx_.set_custom_points_group_visibility(corner_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(X_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(Y_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(Z_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(invalid_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(valid_but_axisless_boundaries_group_index,false);
					break;
				case VIEW_LABELING_GRAPH:
					show_mesh_ = false;
					lighting_ = false;
					mesh_gfx_.unset_scalar_attribute();
					mesh_gfx_.set_facets_colors_by_int_attribute(LABELING_ATTRIBUTE_NAME);
					mesh_gfx_.set_custom_points_group_visibility(corner_group_index,true);
					mesh_gfx_.set_custom_edges_group_color(X_boundaries_group_index,labeling_colors_[0]); // axis X -> color of label 0 = +X
					mesh_gfx_.set_custom_edges_group_color(Y_boundaries_group_index,labeling_colors_[2]); // axis Y -> color of label 2 = +Y
					mesh_gfx_.set_custom_edges_group_color(Z_boundaries_group_index,labeling_colors_[4]); // axis Z -> color of label 4 = +Z
					mesh_gfx_.set_custom_edges_group_visibility(X_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(Y_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(Z_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(invalid_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(valid_but_axisless_boundaries_group_index,false);
					break;
				case VIEW_INVALID_CHARTS:
					show_mesh_ = false;
					lighting_ = false;
					mesh_gfx_.unset_facets_colors_by_int_attribute();
					mesh_gfx_.set_scalar_attribute(MESH_FACETS,"on_invalid_chart",0.0,1.0,TO_GL_TEXTURE_INDEX(COLORMAP_BLUE_RED)); // color according to whether the facet in on an invalid chart or not
					mesh_gfx_.set_custom_points_group_visibility(corner_group_index,false);
					mesh_gfx_.set_custom_edges_group_color(X_boundaries_group_index,labeling_colors_[0]); // axis X -> color of label 0 = +X
					mesh_gfx_.set_custom_edges_group_color(Y_boundaries_group_index,labeling_colors_[2]); // axis Y -> color of label 2 = +Y
					mesh_gfx_.set_custom_edges_group_color(Z_boundaries_group_index,labeling_colors_[4]); // axis Z -> color of label 4 = +Z
					mesh_gfx_.set_custom_edges_group_visibility(X_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(Y_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(Z_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(invalid_boundaries_group_index,false);
					mesh_gfx_.set_custom_edges_group_visibility(valid_but_axisless_boundaries_group_index,false);
					break;
				case VIEW_INVALID_BOUNDARIES:
					show_mesh_ = false;
					lighting_ = false;
					mesh_gfx_.unset_facets_colors_by_int_attribute();
					mesh_gfx_.unset_scalar_attribute();
					mesh_gfx_.set_custom_points_group_visibility(corner_group_index,false);
					mesh_gfx_.set_custom_edges_group_color(X_boundaries_group_index,valid_boundaries_color_);
					mesh_gfx_.set_custom_edges_group_color(Y_boundaries_group_index,valid_boundaries_color_);
					mesh_gfx_.set_custom_edges_group_color(Z_boundaries_group_index,valid_boundaries_color_);
					mesh_gfx_.set_custom_edges_group_visibility(X_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(Y_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(Z_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(invalid_boundaries_group_index,true);
					mesh_gfx_.set_custom_edges_group_visibility(valid_but_axisless_boundaries_group_index,true);
					break;
				default:
					geo_assert_not_reached;
			}
			previous_labeling_visu_mode_ = current_labeling_visu_mode_;
		}

		draw_points();
		draw_surface();
		draw_edges();
		draw_volume();
		mesh_gfx_.draw_custom_points(); // per-group visibility managed in mesh_gfx_
		mesh_gfx_.draw_custom_edges(); // per-group visibility managed in mesh_gfx_
    }

    void LabelingViewerApp::draw_points() {
        if(show_vertices_) {
	    if(vertices_transparency_ != 0.0f) {
		glDepthMask(GL_FALSE);
		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	    }
            mesh_gfx_.set_points_color(
		vertices_color_.x, vertices_color_.y, vertices_color_.z,
		1.0f - vertices_transparency_
	    );
            mesh_gfx_.set_points_size(vertices_size_);
            mesh_gfx_.draw_vertices();

	    if(vertices_transparency_ != 0.0f) {	    
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
	    }
        }

        if(show_vertices_selection_) {
            mesh_gfx_.set_points_color(1.0, 0.0, 0.0);
            mesh_gfx_.set_points_size(2.0f * vertices_size_);
            mesh_gfx_.set_vertices_selection("selection");
            mesh_gfx_.draw_vertices();
            mesh_gfx_.set_vertices_selection("");            
        }
    }

    void LabelingViewerApp::draw_surface() {
	mesh_gfx_.set_mesh_color(0.0, 0.0, 0.0);

	mesh_gfx_.set_surface_color(
	    surface_color_.x, surface_color_.y, surface_color_.z
	);
        if(show_surface_sides_) {
	    mesh_gfx_.set_backface_surface_color(
		surface_color_2_.x, surface_color_2_.y, surface_color_2_.z
	    );
        }
	
        mesh_gfx_.set_show_mesh(show_mesh_);
	mesh_gfx_.set_mesh_color(mesh_color_.x, mesh_color_.y, mesh_color_.z);
	mesh_gfx_.set_mesh_width(index_t(mesh_width_*10.0f));
	
        if(show_surface_) {
	    float specular_backup = glupGetSpecular();
	    glupSetSpecular(0.4f);
            mesh_gfx_.draw_surface();
	    glupSetSpecular(specular_backup);	    
        }
        
        if(show_surface_borders_) {
            mesh_gfx_.draw_surface_borders();
        }
    }

    void LabelingViewerApp::draw_edges() {
        if(show_mesh_) {
            mesh_gfx_.draw_edges();
        }
    }
    
    void LabelingViewerApp::draw_volume() {
        if(show_volume_) {

            if(
                glupIsEnabled(GLUP_CLIPPING) &&
                glupGetClipMode() == GLUP_CLIP_SLICE_CELLS
            ) {
                mesh_gfx_.set_lighting(false);
            }
            
            mesh_gfx_.set_shrink(double(cells_shrink_));
            mesh_gfx_.set_draw_cells(GEO::MESH_HEX, show_hexes_);
            mesh_gfx_.set_draw_cells(GEO::MESH_CONNECTOR, show_connectors_);
	    
            if(show_colored_cells_) {
                mesh_gfx_.set_cells_colors_by_type();
            } else {
                mesh_gfx_.set_cells_color(
		    volume_color_.x, volume_color_.y, volume_color_.z
		);
            }
            mesh_gfx_.draw_volume();

            mesh_gfx_.set_lighting(lighting_);
        }
    }
    
    void LabelingViewerApp::draw_graphics() {
	if(!full_screen_effect_.is_null()) {
	    // Note: on retina, we use window resolution
	    // rather than frame buffer resolution
	    // (we do not need fullres for shadows)
	    full_screen_effect_->pre_render(
		get_width(), get_height()
	    );
	}
	draw_scene_begin();
	draw_scene();
	draw_scene_end();
	if(!full_screen_effect_.is_null()) {
	    full_screen_effect_->post_render();
	}
    }

    void LabelingViewerApp::draw_viewer_properties_window() {
	if(!viewer_properties_visible_) {
	    return;
	}
	if(ImGui::Begin(
	       (icon_UTF8("camera")+" Viewer").c_str(),
	       &viewer_properties_visible_)
	) {
	    draw_viewer_properties();
	}
	ImGui::End();
    }

    void LabelingViewerApp::draw_viewer_properties() {
	if(phone_screen_ && get_height() >= get_width()) {
	    if(props_pinned_) {
		if(ImGui::SimpleButton(
		       icon_UTF8("dot-circle") + "##props_pin"
		)) {
		    props_pinned_ = !props_pinned_;
		}
	    } else {
		if(ImGui::SimpleButton(
		       icon_UTF8("thumbtack") + "##props_pin"
		)) {
		    props_pinned_ = !props_pinned_;
		}
	    }
	    ImGui::SameLine();
	}
	
	if(ImGui::Button(
	       (icon_UTF8("home") + " Home").c_str(),
	       ImVec2(-1.0, 0.0)
	)) {
	    home();
	}
	ImGui::Separator();
	ImGui::Checkbox("Animate", animate_ptr());
	if(three_D_) {
	    ImGui::Checkbox("Lighting", &lighting_);
	    if(lighting_) {
		ImGui::Checkbox("Edit light", &edit_light_);	    
	    }
	    ImGui::Separator();
	    ImGui::Checkbox("Clipping", &clipping_);
	    if(clipping_) {
		ImGui::Combo(
		    "##mode", (int*)&clip_mode_,                
		    "std. GL\0cells\0straddle\0slice\0\0"
		);
		ImGui::Checkbox(
		    "edit clip", &edit_clip_
		);
		ImGui::Checkbox(
		    "fixed clip", &fixed_clip_
		);
	    }
	    ImGui::Separator();
	}
	ImGui::ColorEdit3WithPalette("Background", background_color_.data());
	// Full-screen effects are for now deactivated under Android:
	// SSAO has stripes (and is damned slow)
	// unsharp masking gives a black screen
	if(!phone_screen_) {
	    if(three_D_) {
		if(ImGui::Combo(
		       "sfx",
		       (int*)&effect_,
		       "none\0SSAO\0cartoon\0\0"
		)) {
		    switch(effect_) {
			case 0:
			    full_screen_effect_.reset();
			    break;
			case 1:
			    full_screen_effect_ = new AmbientOcclusionImpl();
			    break;
			case 2:
			    full_screen_effect_ = new UnsharpMaskingImpl();
			    break;
		    }
		}
	    }
	}
    }

    void LabelingViewerApp::draw_object_properties_window() {
	if(!object_properties_visible_) {
	    return;
	}
	if(ImGui::Begin(
	       (icon_UTF8("edit")+" Object").c_str(),
	       &object_properties_visible_)
	) {
	    draw_object_properties();
	}
	ImGui::End();
    }

    void LabelingViewerApp::draw_object_properties() {
		switch (current_state_) {

		case triangle_mesh:
			if(ImGui::Button("Compute naive labeling")) {

				naive_labeling(mesh_,LABELING_ATTRIBUTE_NAME);

				update_static_labeling_graph();

				current_state_ = labeling;
			}
			break;
		case labeling:

			ImGui::RadioButton("View triangle mesh",&current_labeling_visu_mode_,VIEW_TRIANGLE_MESH);

			ImGui::RadioButton("View raw labeling",&current_labeling_visu_mode_,VIEW_RAW_LABELING);

			ImGui::BeginDisabled( (current_labeling_visu_mode_!=VIEW_RAW_LABELING) && (current_labeling_visu_mode_!=VIEW_LABELING_GRAPH) );
			ImGui::ColorEdit4WithPalette("Label 0 = +X", labeling_colors_[0]);
			ImGui::ColorEdit4WithPalette("Label 1 = -X", labeling_colors_[1]);
			ImGui::ColorEdit4WithPalette("Label 2 = +Y", labeling_colors_[2]);
			ImGui::ColorEdit4WithPalette("Label 3 = -Y", labeling_colors_[3]);
			ImGui::ColorEdit4WithPalette("Label 4 = +Z", labeling_colors_[4]);
			ImGui::ColorEdit4WithPalette("Label 5 = -Z", labeling_colors_[5]);
			ImGui::EndDisabled();

			ImGui::RadioButton("View labeling graph",&current_labeling_visu_mode_,VIEW_LABELING_GRAPH);

			ImGui::BeginDisabled(current_labeling_visu_mode_!=VIEW_LABELING_GRAPH);
			ImGui::ColorEdit4WithPalette("Corners", corners_color_);
			ImGui::EndDisabled();

			ImGui::RadioButton(fmt::format("View invalid charts ({})",static_labeling_graph.nb_invalid_charts()).c_str(),&current_labeling_visu_mode_,VIEW_INVALID_CHARTS);

			ImGui::Text("Bleu = valid chart | Red = invalid chart");

			ImGui::RadioButton(fmt::format("View invalid boundaries ({})",static_labeling_graph.nb_invalid_boundaries()).c_str(),&current_labeling_visu_mode_,VIEW_INVALID_BOUNDARIES);

			ImGui::BeginDisabled(current_labeling_visu_mode_!=VIEW_INVALID_BOUNDARIES);
			ImGui::ColorEdit4WithPalette("Valid boundaries", valid_boundaries_color_);
			ImGui::ColorEdit4WithPalette("Invalid boundaries", invalid_boundaries_color_);
			ImGui::EndDisabled();
			
			break;
		
		default:
			break;
		}
    }

    void LabelingViewerApp::draw_command_window() {
	if(Command::current() == nullptr) {
	    return;
	}
	if(!Command::current()->is_visible()) {
	    Command::reset_current();
	    return;
	}
        if(ImGui::Begin(
	    "Command",
            Command::current()->is_visible_ptr()
        )) {
	    Command::current()->draw();
	}
        ImGui::End();
    }

    void LabelingViewerApp::draw_console() {
	console_->draw(&console_visible_);
    }
    
    void LabelingViewerApp::draw_menu_bar() {
	if(!menubar_visible_) {
	    return;
	}

	if(phone_screen_) {
	    if(ImGui::BeginMainMenuBar()) {

		float w = ImGui::GetContentRegionAvail().x;
		if(ImGui::BeginMenu(icon_UTF8("ellipsis-v"))) {
		    draw_application_menus();
		    draw_fileops_menu();
		    draw_about();
		    ImGui::EndMenu();
		}
		w -= ImGui::GetContentRegionAvail().x; // gets btn size

		ImGui::Dummy(ImVec2(0.5f*w, 1.0));
		
		if(supported_read_file_extensions() != "") {
		    if(ImGui::SimpleButton(
			   icon_UTF8("folder-open") + "##menubar_open")
		    ) {
			ImGui::OpenFileDialog(
			    "##load_dlg",
			    supported_read_file_extensions().c_str(),
			    filename_,
			    ImGuiExtFileDialogFlags_Load		
			);
		    }
		}

                if(supported_write_file_extensions() != "") {
		    if(ImGui::SimpleButton(
			   icon_UTF8("save") + "##menubar_save")
		    ) {
			ImGui::OpenFileDialog(
			    "##save_dlg",
			    supported_write_file_extensions().c_str(),
			    filename_,
			    ImGuiExtFileDialogFlags_Save
			);
		    }
		}

		draw_application_icons();
		
		if(use_text_editor_) {
		    ImGui::Dummy(ImVec2(0.5f*w, 1.0));
		    if(ImGui::SimpleButton(icon_UTF8("keyboard"))) {
#ifdef GEO_OS_ANDROID			
			AndroidUtils::show_soft_keyboard(
			    CmdLine::get_android_app()
			);
#endif			
		    }
		}

		ImGui::Dummy(ImVec2(ImGui::GetContentRegionAvail().x - w, 1.0));
		
		if(ImGui::BeginMenu(icon_UTF8("bars"))) {
		    draw_windows_menu();
		    ImGui::EndMenu();
		}

		ImGui::EndMainMenuBar();
	    }
	    return;
	} 
	
        if(ImGui::BeginMainMenuBar()) {
            if(ImGui::BeginMenu("File")) {
                if(supported_read_file_extensions() != "") {
                    draw_load_menu();
                }
#ifndef GEO_OS_EMSCRIPTEN		
		if(current_file_ != "") {
		    if(ImGui::MenuItem(icon_UTF8("save") + " Save")) {
			if(save(current_file_)) {
			    Logger::out("I/O") << "Saved "
					       << current_file_ << std::endl;
			} else {
			    Logger::out("I/O") << "Could not save "
					       << current_file_ << std::endl;
			}
		    }
		}
#endif		
                if(supported_write_file_extensions() != "") {
                    draw_save_menu();
                }
		draw_fileops_menu();
#ifndef GEO_OS_EMSCRIPTEN                        
                ImGui::Separator();
                if(ImGui::MenuItem(icon_UTF8("door-open") + " quit",
				   "[q]", false, true)
		) {
                    this->stop();
                }
#endif
                draw_about();
                // Key shortcuts not really relevant on Android		
		if(!phone_screen_) {
		    draw_help();
		}
                ImGui::EndMenu();
            }
	    if(ImGui::BeginMenu("Windows")) {
		draw_windows_menu();
		ImGui::EndMenu();
	    }
            draw_application_menus();
            
            ImGui::EndMainMenuBar();            
	}
    }

    void LabelingViewerApp::draw_load_menu() {
#ifdef GEO_OS_EMSCRIPTEN
            ImGui::Text("To load a file,");
            ImGui::Text("use the \"Browse\"");
            ImGui::Text("button on the top");
            ImGui::Text("(or \"recent files\"");
            ImGui::Text("below)");
            ImGui::Separator();
            if(ImGui::BeginMenu("Recent files...")) {
                browse(path_);
                ImGui::EndMenu(); 
            }
#else
	    if(ImGui::MenuItem(icon_UTF8("folder-open") + " Load...")) {
		ImGui::OpenFileDialog(
		    "##load_dlg",
		    supported_read_file_extensions().c_str(),
		    filename_,
		    ImGuiExtFileDialogFlags_Load		
		);
	    }
#endif        
    }

    void LabelingViewerApp::draw_save_menu() {
#ifdef GEO_OS_EMSCRIPTEN
        if(ImGui::BeginMenu(icon_UTF8("save") + " Save as...")) {
	    ImGui::MenuItem("Supported extensions:", nullptr, false, false);
            std::vector<std::string> extensions;
            String::split_string(
                supported_write_file_extensions(), ';', extensions
            );
            for(index_t i=0; i<extensions.size(); ++i) {
		ImGui::MenuItem(
		    " ." + extensions[i], nullptr, false, false
		);	    		
	    }
	    ImGui::Separator();
	    static char buff[geo_imgui_string_length];
	    if(current_file_ != "") {
		strcpy(buff, current_file_.c_str());
	    } else if (extensions.size() != 0) {
		strcpy(buff, ("out." + extensions[0]).c_str());		
	    }

	    if(ImGui::InputText(
		   "##MenuFileName",buff,geo_imgui_string_length,
		   ImGuiInputTextFlags_EnterReturnsTrue)
	    ) {
		current_file_ = buff;
		if(String::string_starts_with(current_file_, "/")) {
		    current_file_ = current_file_.substr(
			1,current_file_.length()-1
		    );
		}
		if(save(current_file_)) {
		    std::string command =
			"saveFileFromMemoryFSToDisk(\'" +
			current_file_ +
			"\');" ;
		    emscripten_run_script(command.c_str());
		}
	    }
            ImGui::EndMenu();
        }
#else        
        if(ImGui::MenuItem(icon_UTF8("save") + " Save as...")) {
	    ImGui::OpenFileDialog(
		"##save_dlg",
		supported_write_file_extensions().c_str(),
		filename_,
		ImGuiExtFileDialogFlags_Save
	    );
        }
#endif        
    }

    void LabelingViewerApp::draw_fileops_menu() {
    }

    void LabelingViewerApp::draw_about() {
        ImGui::Separator();
        if(ImGui::BeginMenu(icon_UTF8("info") + " About...")) {
            ImGui::Text("%s : a GEOGRAM application", name().c_str());
	    float sz = float(280.0 * std::min(scaling(), 2.0));
            ImGui::Image(
                convert_to_ImTextureID(geogram_logo_texture_),
                ImVec2(sz, sz)
            );
            ImGui::Text("\n");            
            ImGui::Separator();
            ImGui::Text("\n");
            ImGui::Text("GEOGRAM website: ");
            ImGui::Text("https://github.com/BrunoLevy/geogram");

            ImGui::Text("\n");
            ImGui::Separator();
            ImGui::Text(
                "%s",
                (
                    "GEOGRAM version:" +
                    Environment::instance()->get_value("version")
                ).c_str()
            );
            ImGui::EndMenu();
        }
    }

    void LabelingViewerApp::draw_help() {
	if(ImGui::BeginMenu(icon_UTF8("question") + " help")) {
	    ImGui::Text("Key shortcuts");
	    ImGui::Separator();
	    std::vector<std::string> helps;
	    for(auto it: key_funcs_help_) {
		helps.push_back(it.first + " : " + it.second);
	    }
	    std::sort(helps.begin(), helps.end());
	    for(const std::string& s: helps) {
		ImGui::Text("%s",s.c_str());
	    }
	    ImGui::EndMenu();
	}
    }
    
    void LabelingViewerApp::draw_windows_menu() {
	if(phone_screen_) {
	    ImGui::MenuItem(
		"     " + icon_UTF8("window-restore") + " Windows",
		nullptr,
		false, false
	    );
	}
	if(use_text_editor_) {
	    ImGui::MenuItem(
		icon_UTF8("code") + " text editor",
		phone_screen_ ? nullptr : "[F6]",
		&text_editor_visible_
	    );
	}
	ImGui::MenuItem(
	    icon_UTF8("eye") + " viewer properties",
	    phone_screen_ ? nullptr : "[F7]",
	    &viewer_properties_visible_
	);
	ImGui::MenuItem(
	    icon_UTF8("edit") + " object properties",
	    phone_screen_ ? nullptr : "[F8]",
	    &object_properties_visible_
	);
	ImGui::MenuItem(
	    icon_UTF8("terminal") + " console",
	    phone_screen_ ? nullptr : "[F9]",
	    &console_visible_
	);
	if(!phone_screen_) {
	    ImGui::MenuItem(
		icon_UTF8("bars") + " menubar", "[F12]", &menubar_visible_
	    );
	}
	ImGui::Separator();

	{
	    bool needs_to_close = false;
	    if(phone_screen_) {
		ImGui::MenuItem(
		    "     " + icon_UTF8("font") + " Font size",
		    nullptr,
		    false, false
		);
	    } else {
		needs_to_close =
		    ImGui::BeginMenu(icon_UTF8("font") + " Font size");
	    }

	    if(phone_screen_ || needs_to_close) {
		static index_t font_sizes[] = {10, 12, 14, 16, 18, 22};
		for(index_t i=0; i<sizeof(font_sizes)/sizeof(int); ++i) {
		    bool selected = (get_font_size() == font_sizes[i]);
		    if(ImGui::MenuItem(
			   String::to_string(font_sizes[i]),
			   nullptr,
			   &selected
		    )) {
			set_font_size(font_sizes[i]);
		    }
		}
		if(needs_to_close) {
		    ImGui::EndMenu();
		}
	    }
	}
	if(phone_screen_) {	
	    ImGui::Separator();
	}
	{
	    bool needs_to_close = false;
	    if(phone_screen_) {
		ImGui::MenuItem(
		    "     " + icon_UTF8("palette") + " Style",
		    nullptr,
		    false, false
		);
	    } else {
		needs_to_close = ImGui::BeginMenu(icon_UTF8("cog") + " Style");
	    }
	    if(phone_screen_ || needs_to_close) {
		std::vector<std::string> styles;
		String::split_string(get_styles(), ';', styles);
		for(index_t i=0; i<styles.size(); ++i) {
		    bool selected = (get_style() == styles[i]);
		    if(ImGui::MenuItem(styles[i], nullptr, &selected)) {
			set_style(styles[i]);
		    }
		}
		if(needs_to_close) {
		    ImGui::EndMenu();
		}
	    }
	}
	ImGui::Separator();
	if(ImGui::MenuItem(icon_UTF8("undo") +  " Restore layout")) {
	    set_default_layout();
	}
	if(CmdLine::get_arg_bool("gui:expert")) {
	    if(ImGui::MenuItem("Test android vertical layout")) {
		set_gui_state(default_layout_android_vertical());
	    }
	    if(ImGui::MenuItem("Test android horizontal layout")) {
		set_gui_state(default_layout_android_horizontal());		
	    }
	    if(ImGui::MenuItem("Export gui state to C++")) {
		std::string filename =
		    FileSystem::get_current_working_directory() +
		    "/gui_state.h";
		Logger::out("GUI")
		    << "Exporting current GUI state to C++ file: "
		    << filename
		    << std::endl;
		std::string state = get_gui_state();
		std::ofstream out("gui_state.h");
		out << "// Serialized ImGui windows docking configuration"
		    << std::endl;
		out << "// generated using <geogram_program> gui:expert=true"
		    << std::endl;
		out << "// then Windows->Export gui state to C++"
		    << std::endl;
		out << "const unsigned char gui_state[] = {";
		for(size_t i=0; i<state.length(); ++i) {
		    if((i%10) == 0) {
			out << std::endl;
		    }
		    int x = int((unsigned char)state[i]);
 		    out << x << ",";
 		}
	        out << " 0 };" << std::endl;
	    }
	}
    }
    
    void LabelingViewerApp::set_default_layout() {
	set_gui_state(default_layout());
    }

    const char* LabelingViewerApp::default_layout_android_vertical() const {
	return (const char*)gui_state_v;
    }

    const char* LabelingViewerApp::default_layout_android_horizontal() const {
	return (const char*)gui_state_h;
    }
    
    const char* LabelingViewerApp::default_layout() const {
	const char* result = nullptr;
	if(phone_screen_) {
	    if(get_height() >= get_width()) {
		result = default_layout_android_vertical();
	    } else {
		result = default_layout_android_horizontal();
	    }
	} else {
	    result = (const char*)gui_state;
	}
	return result;
    }

    
    void LabelingViewerApp::resize(
	index_t w, index_t h, index_t fb_w, index_t fb_h
    ) {
	Application::resize(w,h,fb_w,fb_h);
	if(phone_screen_) {
	    set_default_layout();
	}
    }
    
    void LabelingViewerApp::draw_application_menus() {
    }

    void LabelingViewerApp::draw_application_icons() {
	if(phone_screen_) {
	    if(ImGui::SimpleButton(icon_UTF8("sliders-h"))) {
		if(object_properties_visible_) {
		    object_properties_visible_ = false;
		    viewer_properties_visible_ = false;
		} else {
		    object_properties_visible_ = true;
		    viewer_properties_visible_ = true;
		}
	    }
	}
    }
    
    void LabelingViewerApp::post_draw() {
	Command::flush_queue();
    }

    void LabelingViewerApp::mouse_button_callback(
	int button, int action, int mods, int source
    ) {
	geo_argused(mods);
	
	// Hide object and viewer properties if in phone
	// mode and user clicked elsewhere.
	if(phone_screen_ &&
	   !ImGui::GetIO().WantCaptureMouse &&
	   get_height() >= get_width()
	) {
	    if(!props_pinned_) {
		object_properties_visible_ = false;
		viewer_properties_visible_ = false;
	    }
	}

	
	// Swap "buttons" if using fingers (it is more
	// natural to do the rotation with one finger,
	// zoom with two fingers and then translation)
	// Same thing if editing light and event source
	// is stylus.
	if(
	    source == EVENT_SOURCE_FINGER ||
	    (lighting_ && source == EVENT_SOURCE_STYLUS && edit_light_)
	) {
	    if(button == 0) {
		button = 1;
	    } else if(button == 1) {
		button = 0;
	    }
	}
	
	if(action == EVENT_ACTION_DOWN) {
	    mouse_down_xy_ = mouse_xy_;
	    if(button == 1) {
		if(three_D_) {
		    mouse_op_ = MOUSE_ROTATE;
		} else {
		    mouse_op_ = MOUSE_TRANSLATE;
		}
	    } else if(button == 0) {
		mouse_op_ = MOUSE_TRANSLATE;
	    } else if(button == 2) {
		mouse_op_ = MOUSE_ZOOM;
		zoom_down_ = zoom_;
	    }
	    if(clipping_ && edit_clip_) {
		mouse_target_ = MOUSE_CLIP;
	    } else if(lighting_ && edit_light_) {
		mouse_target_ = MOUSE_LIGHT;
	    } else {
		mouse_target_ = MOUSE_OBJECT;
	    }
	    if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
		switch(mouse_target_) {
		    case MOUSE_NOTARGET:
			break;
		    case MOUSE_OBJECT:
			object_rotation_.grab(mouse_xy_);
			if(fixed_clip_) {
			    clip_rotation_.grab(mouse_xy_); 
			}
			break;
		    case MOUSE_LIGHT:
			light_rotation_.grab(mouse_xy_);		    
			break;
		    case MOUSE_CLIP:
			clip_rotation_.grab(mouse_xy_);		    
			break;
		}
	    }
	} else if(action == EVENT_ACTION_UP) {
	    if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
		switch(mouse_target_) {
		    case MOUSE_NOTARGET:
			break;
		    case MOUSE_OBJECT:
			object_rotation_.release(mouse_xy_);
			if(fixed_clip_) {
			    clip_rotation_.release(mouse_xy_); 
			}
			break;
		    case MOUSE_LIGHT:
			light_rotation_.release(mouse_xy_);		    
			break;
		    case MOUSE_CLIP:
			clip_rotation_.release(mouse_xy_);		    
			break;
		}
	    }
	    mouse_op_ = MOUSE_NOOP;
	    mouse_target_ = MOUSE_NOTARGET;
	}
    }

    void LabelingViewerApp::cursor_pos_callback(
	double x, double y, int source
    ) {
	geo_argused(source);
	mouse_xy_ = vec2(
	    double(x) / double(get_width()),
	    double(y) / double(get_height())
	);
	mouse_xy_ *= 2.0;
	mouse_xy_ -= vec2(1.0, 1.0);
	if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
	    switch(mouse_target_) {
		case MOUSE_NOTARGET:
		    break;
		case MOUSE_OBJECT:
		    object_rotation_.drag(mouse_xy_);
		    if(fixed_clip_) {
			clip_rotation_.drag(mouse_xy_);			
		    }
		    break;
		case MOUSE_LIGHT:
		    light_rotation_.drag(mouse_xy_);		    
		    break;
		case MOUSE_CLIP:
		    clip_rotation_.drag(mouse_xy_);		    
		    break;
	    }
	} else if(mouse_op_ == MOUSE_ZOOM && mouse_target_ == MOUSE_OBJECT) {
	    double R = mouse_xy_.y - mouse_down_xy_.y;
	    double fact = (1.0 + R);
	    fact = std::min(fact, 2.0);
	    fact = std::max(fact, 0.1);
	    zoom_ = zoom_down_ * fact;
	} else if( mouse_op_ == MOUSE_TRANSLATE) {
	    double dx = mouse_xy_.x - mouse_down_xy_.x;	    	    
	    double dy = mouse_xy_.y - mouse_down_xy_.y;
	    if(mouse_target_ == MOUSE_OBJECT) {
		object_translation_.x += 2.0 * dx / zoom_;
		object_translation_.y -= 2.0 * dy / zoom_;
		if(fixed_clip_) {
		    clip_translation_.x += 2.0 * dx / zoom_;
		    clip_translation_.y -= 2.0 * dy / zoom_;
		}
	    } else if(mouse_target_ == MOUSE_CLIP) {
		clip_translation_.x += 2.0 * dx / zoom_;
		clip_translation_.y -= 2.0 * dy / zoom_;
	    }
	    mouse_down_xy_ = mouse_xy_;
	}
    }

    void LabelingViewerApp::scroll_callback(double xoffset, double yoffset) {
	geo_argused(xoffset);
	double dy = -40.0*double(yoffset) / double(get_height());
	zoom_ *= (1.0 + dy);
    }

    bool LabelingViewerApp::save(const std::string& filename) {
        Logger::warn("GLUP")
	    << "Could not save " << filename << std::endl;
        Logger::warn("GLUP")
	    << "LabelingViewerApp::save() needs to be overloaded"
	    << std::endl;
        return false;
    }
    
    bool LabelingViewerApp::load(const std::string& filename) {

        // based on ext/geogram/src/lib/geogram_gfx/gui/simple_mesh_application.cpp load()

        if(!FileSystem::is_file(filename)) {
			fmt::println(Logger::out("I/O"),"{} is not a file",filename); Logger::err("I/O").flush();
			return false;
        }

		if(FileSystem::extension(filename)=="txt") {

			if(current_state_ == empty) {
				fmt::println(Logger::err("I/O"),"You need to import a triangle mesh before importing a labeling"); Logger::err("I/O").flush();
				return false;
			}

			if(!load_labeling(filename,mesh_,LABELING_ATTRIBUTE_NAME)) {
				// Should the labeling be removed ?
				// If a labeling was already displayed, it should be restored...
				mesh_gfx_.clear_custom_drawings();
				current_state_ = triangle_mesh;
				return false;
			}

			update_static_labeling_graph();

			current_state_ = labeling;

			return true;
		}

        mesh_gfx_.set_mesh(nullptr);

        mesh_.clear(true,false); // keep_attributes=true is very important. else runtime error when re-importing files
        
        if(GEO::CmdLine::get_arg_bool("single_precision")) {
            mesh_.vertices.set_single_precision();
        }
        
        MeshIOFlags flags;
        if(CmdLine::get_arg_bool("attributes")) {
            flags.set_attribute(MESH_FACET_REGION);
            flags.set_attribute(MESH_CELL_REGION);            
        } 
        if(!mesh_load(filename, mesh_, flags)) {
            return false;
        }

		if(!mesh_.facets.are_simplices() || mesh_.cells.nb()!=0) { // check if it is a triangle mesh
			fmt::println(Logger::err("I/O"),"This file does not contain a triangle mesh. Only surface triangle meshes are supported.");  Logger::err("I/O").flush();
			current_state_ = empty;
			// Instead of clearing the mesh, we should load the mesh elsewhere, check the new mesh, then replace the displayed mesh
			return false;
		}

        if(
            FileSystem::extension(filename) == "obj6" ||
            FileSystem::extension(filename) == "tet6"
        ) {
            Logger::out("Vorpaview")
                << "Displaying mesh animation." << std::endl;

	    start_animation();
            
            mesh_gfx_.set_animate(true);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(mesh_, xyzmin, xyzmax, true);
            set_region_of_interest(
                xyzmin[0], xyzmin[1], xyzmin[2],
                xyzmax[0], xyzmax[1], xyzmax[2]
            );
        } else {
            mesh_gfx_.set_animate(false);            
            mesh_.vertices.set_dimension(3);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(mesh_, xyzmin, xyzmax, false);
            set_region_of_interest(
                xyzmin[0], xyzmin[1], xyzmin[2],
                xyzmax[0], xyzmax[1], xyzmax[2]
            );
        }

        show_vertices_ = (mesh_.facets.nb() == 0);
        mesh_gfx_.set_mesh(&mesh_);

		current_state_ = triangle_mesh;

	    current_file_ = filename;
        return true;
    }
    
    bool LabelingViewerApp::can_load(const std::string& filename) {
        std::string extensions_str = supported_read_file_extensions();
        if(extensions_str == "") {
            return false;
        }
        if(extensions_str == "*") {
            return true;
        }
        std::string extension = FileSystem::extension(filename);
        std::vector<std::string> extensions;
        String::split_string(extensions_str, ';', extensions);
        for(index_t i=0; i<extensions.size(); ++i) {
            if(extensions[i] == extension) {
                return true;
            }
        }
        return false;
    }
    
    std::string LabelingViewerApp::supported_read_file_extensions() {
        return "";
    }

    std::string LabelingViewerApp::supported_write_file_extensions() {
        return "";
    }

    ImTextureID LabelingViewerApp::convert_to_ImTextureID(
	GLuint gl_texture_id_in
    ) {
        // It is not correct to directly cast a GLuint into a void*
        // (generates warnings), therefore I'm using a union.
        union {
            GLuint gl_texture_id;
            ImTextureID imgui_texture_id;
        };
        imgui_texture_id = nullptr;
        gl_texture_id = gl_texture_id_in;
        return imgui_texture_id;
    }

    void LabelingViewerApp::GL_initialize() {
	Application::GL_initialize();
        glGenTextures(1, &geogram_logo_texture_);
        glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
        glBindTexture(GL_TEXTURE_2D, geogram_logo_texture_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2Dxpm(geogram_logo_256_xpm);
	init_colormaps();
	std::string keys = CmdLine::get_arg("gfx:keypress");
	for(index_t i=0; i<index_t(keys.length()); ++i) {
	    char_callback((unsigned int)(keys[i]));
	}
    }
    
    void LabelingViewerApp::GL_terminate() {
        for(index_t i=0; i<colormaps_.size(); ++i) {
            if(colormaps_[i].texture != 0) {
                glDeleteTextures(1, &colormaps_[i].texture);
		colormaps_[i].texture = 0;
            }
        }
        if(geogram_logo_texture_ != 0) {
            glDeleteTextures(1, &geogram_logo_texture_);
	    geogram_logo_texture_ = 0;
        }
	Application::GL_terminate();
    }


    void LabelingViewerApp::browse(const std::string& path, bool subdirs) {
        std::vector<std::string> files;
        FileSystem::get_directory_entries(path,files);
        
        for(index_t i=0; i<files.size(); ++i) {
            if(FileSystem::is_directory(files[i]) && subdirs) {
                if(ImGui::BeginMenu(path_to_label(path_,files[i]))) {
                    browse(files[i]);
                    ImGui::EndMenu();
                }
            } else {
                if(can_load(files[i])) {
                    if(ImGui::MenuItem(path_to_label(path_,files[i]))) {
                        load(files[i]);
                    }
                }
            }
        }
    }
    
    void LabelingViewerApp::geogram_initialize(int argc, char** argv) {
	GEO::Application::geogram_initialize(argc, argv);
        if(filenames().size() == 1 &&
	   FileSystem::is_directory(filenames()[0])
	) {
            path_ = filenames()[0];
        } else if(filenames().size() > 0) {
            for(index_t i=0; i<filenames().size(); ++i) {
                load(filenames()[i]);
            }
            if(filenames().size() > 0) {
                path_ = FileSystem::dir_name(filenames()[filenames().size()-1]);
            }
        } else {
            path_ = FileSystem::documents_directory();
        }
        Logger::instance()->register_client(console_);
        Progress::set_client(status_bar_);
	set_default_layout();
	if(phone_screen_) {
	    object_properties_visible_ = false;
	    viewer_properties_visible_ = false;
	}
    }

    void LabelingViewerApp::init_colormap(
        const std::string& name, const char** xpm_data
    ) {
        colormaps_.push_back(ColormapInfo());
        colormaps_.rbegin()->name = name;
        glGenTextures(1, &colormaps_.rbegin()->texture);
        glBindTexture(GL_TEXTURE_2D, colormaps_.rbegin()->texture);
        glTexImage2Dxpm(xpm_data);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
        );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void LabelingViewerApp::init_colormaps() {
        for(int i = 0; i <= LAST_COLORMAP; i++) {
			init_colormap(colormap_name[i],colormap_xpm[i]);
		}
    }

    bool LabelingViewerApp::exec_command(const char* command) {
#ifdef GEOGRAM_WITH_LUA	
	if(luaL_dostring(lua_state_,command)) {
	    adjust_lua_glup_state(lua_state_);
	    const char* msg = lua_tostring(lua_state_,-1);
	    const char* msg2 = strchr(msg,']');
	    if(msg2 != nullptr) {
		msg = msg2+2;
	    }
	    Logger::err("LUA") << "line " << msg << std::endl;
	    lua_error_occured_ = true;
	} else {
	    lua_error_occured_ = false;
	}
	return !lua_error_occured_;
#else
	geo_argused(command);
	Logger::err("LUA") << "Compiled without LUA support"
			   << std::endl;
	return false;
#endif	
    }

    void LabelingViewerApp::ImGui_initialize() {
	if(phone_screen_ && !ImGui_firsttime_init_) {
	    console_visible_ = false;
	    // Tooltips do not play well with touch screens.
	    ImGui::DisableTooltips();
#ifndef GEO_OS_ANDROID
	    set_default_layout();
#endif	    
	}
	CmdLine::set_arg("gui:style","Light");
	Application::ImGui_initialize();
    }

    void LabelingViewerApp::drop_callback(int nb, const char** f) {
	for(int i=0; i<nb; ++i) {
	    load(f[i]);
	}
    }

    void LabelingViewerApp::get_bbox(
        const Mesh& M_in, double* xyzmin, double* xyzmax, bool animate
    ) {
        geo_assert(M_in.vertices.dimension() >= index_t(animate ? 6 : 3));
        for(index_t c = 0; c < 3; c++) {
            xyzmin[c] = Numeric::max_float64();
            xyzmax[c] = Numeric::min_float64();
        }

        for(index_t v = 0; v < M_in.vertices.nb(); ++v) {
            if(M_in.vertices.single_precision()) {
                const float* p = M_in.vertices.single_precision_point_ptr(v);
                for(coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = std::min(xyzmin[c], double(p[c]));
                    xyzmax[c] = std::max(xyzmax[c], double(p[c]));
                    if(animate) {
                        xyzmin[c] = std::min(xyzmin[c], double(p[c+3]));
                        xyzmax[c] = std::max(xyzmax[c], double(p[c+3]));
                    }
                }
            } else {
                const double* p = M_in.vertices.point_ptr(v);
                for(coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = std::min(xyzmin[c], p[c]);
                    xyzmax[c] = std::max(xyzmax[c], p[c]);
                    if(animate) {
                        xyzmin[c] = std::min(xyzmin[c], p[c+3]);
                        xyzmax[c] = std::max(xyzmax[c], p[c+3]);
                    }
                }
            }
        }
    }

	void LabelingViewerApp::update_static_labeling_graph() {

		// compute charts, boundaries and corners of the labeling
		static_labeling_graph.fill_from(mesh_,LABELING_ATTRIBUTE_NAME,false);

		fmt::println(Logger::out("I/O"),"There are {} charts, {} corners and {} boundaries in this labeling.",static_labeling_graph.nb_charts(),static_labeling_graph.nb_corners(),static_labeling_graph.nb_boundaries());  Logger::out("I/O").flush();

		static_labeling_graph.dump_to_file("StaticLabelingGraph.txt");

		mesh_gfx_.clear_custom_drawings();

		corner_group_index = mesh_gfx_.new_custom_points_group(corners_color_,true);
		X_boundaries_group_index = mesh_gfx_.new_custom_edges_group(labeling_colors_[0],false); // axis X -> color of label 0 = +X
		Y_boundaries_group_index = mesh_gfx_.new_custom_edges_group(labeling_colors_[2],false); // axis Y -> color of label 2 = +Y
		Z_boundaries_group_index = mesh_gfx_.new_custom_edges_group(labeling_colors_[4],false); // axis Z -> color of label 4 = +Z
		invalid_boundaries_group_index = mesh_gfx_.new_custom_edges_group(invalid_boundaries_color_,false);
		valid_but_axisless_boundaries_group_index = mesh_gfx_.new_custom_edges_group(valid_boundaries_color_,false);

		// store the coordinates of each corner in mesh_gfx_
		// to be able to draw them
		for(std::size_t i = 0; i < static_labeling_graph.nb_corners(); ++i) {
			const double* coordinates = mesh_.vertices.point_ptr(
				static_labeling_graph.corners[i].vertex
			);
			mesh_gfx_.add_custom_point(corner_group_index,coordinates[0], coordinates[1], coordinates[2]);
		}

		// store the coordinates of each boundary edge in mesh_gfx_
		// to be able to draw them
		std::size_t group_index;
		for(std::size_t i = 0; i < static_labeling_graph.nb_boundaries(); ++i) {
			const Boundary& boundary = static_labeling_graph.boundaries[i];
			for(const auto& be : boundary.halfedges) { // for each boundary edge of this boundary
				const double* coordinates_first_point = halfedge_vertex_from(mesh_,be).data();
				const double* coordinates_second_point = halfedge_vertex_to(mesh_,be).data();
				switch(boundary.axis) {
					case -1:
						// May or may not be invalid. We must look at static_labeling_graph.invalid_boundaries
						if(VECTOR_CONTAINS(static_labeling_graph.invalid_boundaries,i)) {
							group_index = invalid_boundaries_group_index;
						}
						else {
							group_index = valid_but_axisless_boundaries_group_index;
						}
						break;
					case 0:
						group_index = X_boundaries_group_index;
						break;
					case 1:
						group_index = Y_boundaries_group_index;
						break;
					case 2:
						group_index = Z_boundaries_group_index;
						break;
					default:
						geo_assert_not_reached;
				}
				mesh_gfx_.add_custom_edge(
					group_index,
					coordinates_first_point[0], coordinates_first_point[1], coordinates_first_point[2], // first point
					coordinates_second_point[0], coordinates_second_point[1], coordinates_second_point[2] // second point
				);
			}
		}

		mesh_gfx_.set_custom_points_group_visibility(corner_group_index,true);
		mesh_gfx_.set_custom_edges_group_visibility(X_boundaries_group_index,true);
		mesh_gfx_.set_custom_edges_group_visibility(Y_boundaries_group_index,true);
		mesh_gfx_.set_custom_edges_group_visibility(Z_boundaries_group_index,true);
	}