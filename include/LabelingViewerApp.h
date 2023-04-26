// copied from ext/geogram/src/lib/geogram_gfx/simple_application.h
//
// # Changelog - https://keepachangelog.com
//
// ## Added
//
// - inclusion of <geogram_gfx/mesh/mesh_gfx.h>
// - get_bbox(), copied from ext/geogram/src/lib/geogram_gfx/gui/simple_mesh_application.h
// - draw_points() draw_edges() draw_surface() draw_volume(), copied from ext/geogram/src/lib/geogram_gfx/gui/simple_mesh_application.h
// - mesh_, anim_speed_, anim_time_, show_vertices_, show_vertices_selection_, vertices_size_, vertices_color_,
//   vertices_transparency_, show_surface_, show_surface_sides_, show_mesh_, mesh_width_, mesh_color_, show_surface_borders_,
//   surface_color_, surface_color_2_, show_volume_, cells_shrink_, volume_color_, show_colored_cells_, show_hexes_, show_connectors_,
//   show_attributes_, current_colormap_texture_, attribute_, attribute_subelements_, attribute_name_, attribute_min_, attribute_max_ attributes
// - LABELING_ATTRIBUTE_NAME macro
// - show_labeling_, labeling_colors_ and axis_colors_ attributes
// - inclusion of CustomMeshGfx.h
// - mesh_gfx_ attribute of type CustomMeshGfx
//
// ## Changed
// 
// - class SimpleApplication renamed to LabelingViewerApp
// - TextEditor is now GEO::TextEditor
// - the constructor has no arguments anymore. The application name is hard-coded.
// - load() and draw_scene() no longer virtual
//
// ## Removed
//
// - namespace GEO wrapping the code
// - inclusion of <geogram_gfx/mesh/mesh_gfx.h>

#pragma once

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/gui/application.h>
#include <geogram_gfx/gui/status_bar.h>
#include <geogram_gfx/gui/console.h>
#include <geogram_gfx/gui/text_editor.h>
#include <geogram_gfx/gui/command.h>
#include <geogram_gfx/gui/arc_ball.h>
#include <geogram_gfx/full_screen_effects/full_screen_effect.h>
#include <geogram_gfx/ImGui_ext/imgui_ext.h>
#include <geogram_gfx/ImGui_ext/icon_font.h>

#include <map>
#include <functional>

#include "CustomMeshGfx.h"   // for CustomMeshGfx

#define LABELING_ATTRIBUTE_NAME "label"

struct lua_State;

using namespace GEO;

    class GEOGRAM_GFX_API LabelingViewerApp : public Application {
      public:

	/**
	 * \brief LabelingViewerApp constructor.
	 */
	LabelingViewerApp();

	/**
	 * \brief LabelingViewerApp destructor.
	 */
	~LabelingViewerApp() override;
	
	/**
	 * \copydoc GEO::Application::draw_gui()
	 */
	void draw_gui() override;

	/**
	 * \copydoc GEO::Application::draw_graphics()
	 */
	void draw_graphics() override;

        /**
         * \brief Saves the current content to a file.
         * \details Baseclass implementation does nothing. Derived classes
         *  may overload this function.
         * \retval true if the file could be sucessfully saved.
         * \retval false otherwise
         */
        virtual bool save(const std::string& filename);
        
        /**
         * \brief Loads a file.
         * \retval true if the file could be sucessfully loaded
         * \retval false otherwise
         */
        bool load(const std::string& filename);

	/**
	 * \brief Gets the text editor.
	 * \return a reference to the text editor.
	 */
	GEO::TextEditor& text_editor() {
	    return text_editor_;
	}

	/**
	 * \brief Shows the text editor.
	 */
	void show_text_editor() {
	    text_editor_visible_ = true;
	}

	/**
	 * \brief Hides the text editor.
	 */
	void hide_text_editor() {
	    text_editor_visible_ = false;
	}

	/**
	 * \brief Shows the console.
	 */
	void show_console() {
	    console_visible_ = true;
	}

	/**
	 * \brief Hides the console.
	 */
	void hide_console() {
	    console_visible_ = false;
	}

	/**
	 * \brief Restores default viewing parameters.
	 */
	void home();

	/**
	 * \copydoc GEO::Application::set_style()
	 */
	void set_style(const std::string& style) override;

	/**
	 * \brief Sets the region of interest
	 * \details This defines the default target of the camera
	 * \param[in] xmin , ymin , zmin , xmax , ymax , zmax the
	 *   bounds of the region of interest.
	 */
	void set_region_of_interest(
	    double xmin, double ymin, double zmin,
	    double xmax, double ymax, double zmax
	);

	/**
	 * \brief Gets the region of interest
	 * \see set_region_of_interest()
	 * \param[out] xmin , ymin , zmin , xmax , ymax , zmax the
	 *   bounds of the region of interest.
	 */
	void get_region_of_interest(
	    double& xmin, double& ymin, double& zmin,
	    double& xmax, double& ymax, double& zmax
	) const;
	
	void zoom_in() {
	    zoom_ *= 1.1;
	}

	void zoom_out() {
	    zoom_ /= 1.1;
	}

	void set_clipping(bool x) {
	    clipping_ = x;
	}

	void set_lighting(bool x) {
	    lighting_ = x;
	}

	void set_background_color(const vec4f& color) {
	    background_color_ = color;
	}

	virtual bool exec_command(const char* command);
	
	static LabelingViewerApp* instance() {
	    return dynamic_cast<LabelingViewerApp*>(
		GEO::Application::instance()
	    );
	}

	/**
	 * \brief Projects a point from model space to window coordinates.
	 * \param[in] p the point in model space coordinates.
	 * \return the point in window coordinates, that is
	 *    [0,width-1] x [0,height-1]
	 */
	vec3 project(const vec3& p);

	/**
	 * \brief Unprojects a 3d point from window coordinates to model space.
	 * \param[in] p the 3d point in window coordinates, that is
	 *   [0,width-1] x [0,height-1]
	 * \return the point in model space coordinates.
	 */
	vec3 unproject(const vec3& p);

	/**
	 * \brief Unprojects a 2d point from window coordinates to model space.
	 * \param[in] p the 2d point in screen space coordinates.
	 * \return the 2d point in model space coordinates.
	 */
	vec2 unproject_2d(const vec2& p);


	/**
	 * \copydoc Application::drop_callback()
	 */
	void drop_callback(int nb, const char** f) override;

    /**
     * \brief Gets the bounding box of a mesh animation.
     * \details In animated mode, the mesh animation is stored as 
     *  a mesh with 6d coordinates, that correspond to the geometric 
     *  location at the vertices at time 0 and at time 1.
     * \param[in] M_in the mesh
     * \param[out] xyzmin a pointer to the three minimum coordinates
     * \param[out] xyzmax a pointer to the three maximum coordinates
     * \param[in] animate true if displaying a mesh animation
     */
    void get_bbox(
        const Mesh& M_in, double* xyzmin, double* xyzmax, bool animate
    );
	
      protected:

	/**
	 * \brief Declares a function to be triggered when a key is pressed.
	 * \param[in] key the key ("a" for a, "F1" for F1)
	 * \param[in] cb the function to be called
	 * \param[in] help an optional help string
	 */
	void add_key_func(
	    const std::string& key, std::function<void()> cb,
	    const char* help = nullptr
	);

	/**
	 * \brief Declares a boolean to be toggled when a key is pressed.
	 * \param[in] key the key ("a" for a, "F1" for F1)
	 * \param[in] p_val a pointer to the boolean
	 * \param[in] help an optional help string
	 */
	void add_key_toggle(
	    const std::string& key, bool* p_val,
	    const char* help = nullptr
	);

	/**
	 * \copydoc GEO::Application::char_callback()
	 */
	void char_callback(unsigned int c) override;
	
	/**
	 * \copydoc GEO::Application::key_callback()
	 */
        void key_callback(int key, int scancode, int action, int mods) override;

        /**
	 * \copydoc GEO::Application::mouse_button_callback()
	 */
	void mouse_button_callback(
	    int button, int action, int mods, int source
	) override;

        /**
	 * \copydoc GEO::Application::cursor_pos_callback()
	 */
	void cursor_pos_callback(double x, double y, int source) override;

        /**
	 * \copydoc GEO::Application::scroll_callback()
	 */
	void scroll_callback(double xoffset, double yoffset) override;
	
	/**
	 * \brief Setups OpenGL for scene drawing.
	 */
	virtual void draw_scene_begin();

	/**
	 * \brief Draws the scene.
	 */
	void draw_scene();

	/**
	 * \brief Cleanups OpenGL after scene drawing.
	 */
	virtual void draw_scene_end();

    void draw_points();
    void draw_edges();
    void draw_surface();
    void draw_volume();
	
        /**
         * \brief Draws the viewer properties window frame and contents.
         */
        virtual void draw_viewer_properties_window();
        
        /**
         * \brief Draws the contents of viewer properties window.
         */
        virtual void draw_viewer_properties();

        /**
         * \brief Draw the object properties window frame and contents.
         */
        virtual void draw_object_properties_window();
        
        /**
         * \brief Draws the contents of the object properties window.
         */
        virtual void draw_object_properties();


        /**
         * \brief Draws the active command window if any.
         */
        virtual void draw_command_window();
	
        /**
         * \brief Draws the console.
         */
        virtual void draw_console();

	
	/**
         * \brief Draws the menu bar.
         */
        virtual void draw_menu_bar();


        /**
         * \brief Draws the load menu and browser.
         */
        virtual void draw_load_menu();

        /**
         * \brief Draws the save menu.
         */
        virtual void draw_save_menu();

	/**
	 * \brief Draws other file operation menu.
	 * \details Default implementation does nothing.
	 *  It can be overloaded to add other menu
	 *  items in the file menu.
	 */
	virtual void draw_fileops_menu();
	
        /**
         * \brief Draws the about box in the file menu.
         */
        virtual void draw_about();

        /**
         * \brief Draws help info (accelarators)
         */
        virtual void draw_help();
	
        /**
         * \brief Draws the windows menu.
         */
        virtual void draw_windows_menu();

        /**
         * \brief Draws the application menus.
         * \details Meant to be overloaded by derived classes.
         */
        virtual void draw_application_menus();

        /**
         * \brief Draws the application icons on the menubar.
         * \details Meant to be overloaded by derived classes.
         */
	virtual void draw_application_icons();
	
	/**
	 * \copydoc Application::post_draw()
	 */
        void post_draw() override;
	
        /**
         * \brief Tests whether a file can be loaded.
         * \details This function can be used to filter the files displayed
         *  in the "Load..." menu. Baseclass implementation always return true.
         *  Derived classes may overload it and return false for files with
         *  unknown extensions.
         */  
        virtual bool can_load(const std::string& filename);
	
        /**
         * \brief Gets the list of supported file extensions for reading.
         * \details This function may be olverloaded by derived class. Base
         *  class implementation returns "". If this function returns "", then
         *  no "Load..." option is displayed in the "File" menu. 
         * \return The semi-colon separated list of supported file extensions,
         *  or "*" if all file extensions are supported.
         */
        virtual std::string supported_read_file_extensions(); 

        /**
         * \brief Gets the list of supported file extensions for writing.
         * \details This function may be olverloaded by derived class. Base
         *  class implementation returns "". If this function returns "", then
         *  no "Save..." option is displayed in the "File" menu.  
         *   If it returns a colon-separated list of extensions, then the
         *  "Save..." option displays a list of possible file names for each
         *  supported extension.
         * \return The semi-colon separated list of supported file extensions,
         *  or "*" if all file extensions are supported.
         */
        virtual std::string supported_write_file_extensions(); 

        /**
         * \brief Converts an OpenGL texture ID into an ImGUI texture ID.
         * \param[in] gl_texture_id_in the OpenGL texture ID
         * \return the corresponding ImGUI texture ID
         */
        ImTextureID convert_to_ImTextureID(GLuint gl_texture_id_in);

	/**
	 * \copydoc GEO::Application::GL_initialize()
	 */
	void GL_initialize() override;

	/**
	 * \copydoc GEO::Application::GL_terminate()
	 */
	void GL_terminate() override;

        /**
         * \brief Recursively browses a directory and generates
         *  menu items.
         * \param[in] path the path to be browsed
         * \param[in] subdirs if true, browse subdirectories as well
         */
        void browse(const std::string& path, bool subdirs=false);

	/**
	 * \copydoc GEO::Application::geogram_initialize()
	 */
	void geogram_initialize(int argc, char** argv) override;

	/**
	 * \brief Sets the default filename used to save
	 *  the current file.
	 * \param[in] filename the default filename.
	 */
	void set_default_filename(const std::string& filename) {
	    strcpy(filename_, filename.c_str());
	}

        /**
         * \brief Initializes a new colormap from name and xpm data.
         * \details This function can be called only once the OpenGL
         *  context is ready, for instance in the init_graphics() function.
         * \param[in] name the name of the colormap
         * \param[in] xpm_data the image data of the colormap.
         */
        void init_colormap(const std::string& name, const char** xpm_data);

        /**
         * \brief Initializes all the default colormaps.
         * \details This function can be called only once the OpenGL
         *  context is ready, for instance in the init_graphics() function.
         */
        void init_colormaps();

	/**
	 * \copydoc Application::ImGui_initialize()
	 */
	void ImGui_initialize() override;

	void set_2d() {
	    three_D_ = false;
	}

	void set_3d() {
	    three_D_ = true;
	}

	void set_default_layout();

        void resize(index_t w, index_t h, index_t fb_w, index_t fb_h) override;

	virtual const char* default_layout() const;
	virtual const char* default_layout_android_vertical() const;
	virtual const char* default_layout_android_horizontal() const;		
	
      protected:
	bool lighting_;
	bool edit_light_;
	bool clipping_;
        GLUPclipMode clip_mode_;
	bool edit_clip_;
	bool fixed_clip_;
	GLenum effect_;
	vec4f background_color_;
	
        bool viewer_properties_visible_;
        bool object_properties_visible_;
        bool console_visible_;
	bool text_editor_visible_;
	bool use_text_editor_;

	Box roi_;
	double roi_radius_;
	vec3    object_translation_;
	ArcBall object_rotation_;
	ArcBall light_rotation_;
	ArcBall clip_rotation_;
	vec3    clip_translation_;
	bool    three_D_;
	double  zoom_;
	double  zoom_down_; /**< Zoom when mouse down. */

	bool props_pinned_;
	
	enum MouseOp {
	    MOUSE_NOOP, MOUSE_ROTATE, MOUSE_TRANSLATE, MOUSE_ZOOM
	} mouse_op_;
	
	enum MouseTarget {
	    MOUSE_NOTARGET, MOUSE_OBJECT, MOUSE_LIGHT, MOUSE_CLIP
	} mouse_target_;

	vec2 mouse_down_xy_; // in [-1,1] x [-1,1]
	vec2 mouse_xy_;      // in [-1,1] x [-1,1]

	// Current transform, for picking
	GLint viewport_[4];
	mat4 modelview_transpose_;
	mat4 project_transpose_;
	
        std::string path_;	
	std::string current_file_;
        char filename_[geo_imgui_string_length]; // Buffer for file dialog.
        GLuint geogram_logo_texture_;

        Console_var console_;
        StatusBar_var status_bar_;
	GEO::TextEditor text_editor_;

	std::map< std::string, std::function<void()> > key_funcs_;
	std::map< std::string, std::string > key_funcs_help_;

        struct ColormapInfo {
            ColormapInfo() : texture(0) {
            }
            GLuint texture;
            std::string name;
        };

        vector<ColormapInfo> colormaps_;
	FullScreenEffectImpl_var full_screen_effect_;

	lua_State* lua_state_;
	bool lua_error_occured_;

    // copied from ext/geogram/src/lib/geogram_gfx/gui/simple_mesh_application.h
    Mesh mesh_;
    CustomMeshGfx mesh_gfx_;

    float anim_speed_;
    float anim_time_;

    bool show_vertices_;
    bool show_vertices_selection_;
    float vertices_size_;
    vec4f vertices_color_;
    float vertices_transparency_;

    bool show_surface_;
    bool show_surface_sides_;   
    bool show_mesh_;
    float mesh_width_;
    vec4f mesh_color_;

    bool show_surface_borders_;
    vec4f surface_color_;
    vec4f surface_color_2_;

    bool show_volume_;
    float cells_shrink_;
    vec4f volume_color_;
    bool show_colored_cells_;
    bool show_hexes_;
    bool show_connectors_;

    bool show_attributes_;
    GLuint current_colormap_texture_;
    std::string       attribute_;
    MeshElementsFlags attribute_subelements_;
    std::string       attribute_name_;
    float             attribute_min_;
    float             attribute_max_;

	// added
	bool show_labeling_;
	float labeling_colors_[6][4];
	float axis_colors_[3][4];
    };