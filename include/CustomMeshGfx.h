// # Changelog - https://keepachangelog.com
//
// Modifications of Geogram's source files
//   ext/geogram/src/lib/geogram_gfx/mesh/mesh_gfx.h
// & ext/geogram/src/lib/geogram_gfx/mesh/mesh_gfx.cpp
//
// ## Rendering of boundary edges
//
// ### Added
//
// - new variable : custom_edges_
// - new functions : add_custom_edge() and draw_custom_edges()
// - inclusion of <utility> and geometry.h
//
// ## Rendering of corners
//
// ### Added
//
// - new variables : custom_points_, custom_points_color_
// - new functions : add_custom_point(), set_custom_points_color()
// - inclusion of <array>
//
// ### Changed
//
// - draw_vertices() also draws the custom points
// - set_mesh() clears the custom points
//
// ## Rendering of facet labels
//
// ### Added
//
// - new variables : facets_colors_by_int_attribute_, value_to_color_, int_attribute_
// - facets_colors_by_int_attribute_ is initialized in the constructor
// - new functions : bind_int_attribute_value_to_color(), set_facets_colors_by_int_attribute(),
//   unset_facets_colors_by_int_attribute(), draw_triangles_immediate_by_int_attrib()
// - inclusion of <map>
//
// ### Changed
//
// - draw_triangles() calls draw_triangles_immediate_by_int_attrib() when facets_colors_by_int_attribute_ is true
//
// ## To have an independent copy of MeshGfx
//
// ### Added
//
// - using namespace GEO
//
// ### Changed
//
// - replace MeshGfx by CustomMeshGfx
// - Geogram's include guard replaced by #pragma once
//
// ### Removed
//
// - namespace GEO wrapping the code
// - inclusion of <geogram_gfx/mesh/mesh_gfx.h>

#pragma once

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/GLUP/GLUP_private.h>
#include <geogram/mesh/mesh.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <map>
#include <array>
#include <utility> // for std::pair

#include "geometry.h" // for vec3d, VecngCompare 

/**
 * \file include/CustomMeshGfx.h
 * \brief A class to display a mesh using OpenGL/GLUP.
 */

using namespace GEO;

    class MeshGfxImpl;
    
    /**
     * \brief Draws a mesh using OpenGL.
     */
    class GEOGRAM_GFX_API CustomMeshGfx {
    public:

        /**
         * \brief CustomMeshGfx constructor.
         */
        CustomMeshGfx();

        /**
         * \brief CustomMeshGfx destructor.
         */
        ~CustomMeshGfx();

        /**
         * \brief Forbids CustomMeshGfx copy..
         */
        CustomMeshGfx(const CustomMeshGfx& rhs) = delete;
        
        /**
         * \brief Forbids CustomMeshGfx copy..
         */
        CustomMeshGfx& operator=(const CustomMeshGfx& rhs) = delete;

        /**
         * \brief Draws the vertices of the mesh.
         * \details If a vertices selection is set, only the
         *  vertices in the selection are drawn
         * \see set_vertices_selection()
         */
        void draw_vertices();

        /**
         * \brief Draws the edges of the mesh.
         */
        void draw_edges();

        void draw_custom_edges();
        
        /**
         * \brief Draws the surfacic part of the mesh.
         */
        void draw_surface();

        /**
         * \brief Draws the borders of the surfacic 
         *  part of the mesh.
         */
        void draw_surface_borders();
        
        /**
         * \brief Draws the volumetric part of the mesh.
         */
        void draw_volume();
        
        /**
         * \brief Gets the mesh visibility flag.
         * \details The mesh visibility flags specifies
         *  whether mesh edges should be drawn. The used
         *  color can be specified by set_mesh_color()
         * \retval true if mesh edges should be displayed
         * \retval false otherwise
         */
        bool get_show_mesh() const {
            return show_mesh_;
        }

        /**
         * \brief Sets the mesh visibility flag.
         * \param[in] x the new value of the mesh visibility flag.
         * \details The mesh visibility flags specifies
         *  whether mesh edges should be drawn. The used
         *  color can be specified by set_mesh_color()
         * \note For now, regions are only implemented for
         *   triangulated meshes and tetrahedralized meshes
         *   (not implemented yet for hybrid surfacic and 
         *    volumetric meshes).
         */
        void set_show_mesh(bool x) {
            show_mesh_ = x;
        }

        /**
         * \brief Gets the mesh width
         * \details The mesh width is taken into account 
         *   when the mesh visibility flag is set 
         *   (by set_show_mesh()), when drawing facets
         *   and cells.
         * \return the mesh width
         */
        index_t get_mesh_width() const {
            return mesh_width_;
        }

        /**
         * \brief Sets the mesh width
         * \details The mesh width is taken into account 
         *   when the mesh visibility flag is set 
         *   (by set_show_mesh()), when drawing facets
         *   and cells.
         * \param[in] x the mesh width (minimum is 1)
         */
        void set_mesh_width(index_t x) {
            mesh_width_ = x;
        }

        /**
         * \brief Gets the mesh border width
         * \details The mesh border width is the one used
         *   by draw_surface_borders()
         * \return the mesh border width
         */
        index_t get_mesh_border_width() const {
            return mesh_border_width_;
        }

        /**
         * \brief Sets the mesh border width
         * \details The mesh border width is the one used
         *   by draw_surface_borders()
         * \param[in] x the mesh width (minimum is 1)
         */
        void set_mesh_border_width(index_t x) {
            mesh_border_width_ = x;
        }
        
        /**
         * \brief Gets the cells shrink coefficient.
         * \details The cells shrink coefficient is used 
         *  to display cells slighly smaller than what they are.
         *  Cells shrinking is only supported in GLSL mode.
         * \return the cells shrink coefficient, betwe 0.0 (no shrink) 
         *  and 1.0 (full shrink)
         */
        double get_shrink() const {
            return shrink_;
        }

        /**
         * \brief Sets the cells shrink coefficient.
         * \details The cells shrink coefficient is used 
         *  to display cells slighly smaller than what they are.
         *  Cells shrinking is only supported in GLSL mode.
         * \param[in] x the cells shrink coefficient, betwe 0.0 (no shrink) 
         *  and 1.0 (full shrink)
         */
        void set_shrink(double x) {
            shrink_ = x;
            if(shrink_ < 0.0) {
                shrink_ = 0.0;
            } else if(shrink_ > 1.0) {
                shrink_ = 1.0;
            }
        }


        /**
         * \brief Gets the animate flag
         * \details When animate mode is activated and the
         *  mesh has 6d vertices, then an animation is displayed.
         *  The first three coordinates correspond to the vertex
         *  position at initial time (t=0). The last three coordinates 
         *  correspond to the vertex position at final time (t=1).
         * \retval true if animation is used
         * \retval false otherwise
         * \see get_time(), set_time()
         */
        bool get_animate() const {
            return animate_;
        }

        /**
         * \brief Gets the animate flag
         * \details When animate mode is activated and the
         *  mesh has 6d vertices, then an animation is displayed.
         *  The first three coordinates correspond to the vertex
         *  position at initial time (t=0). The last three coordinates 
         *  correspond to the vertex position at final time (t=1).
         * \param[in] x true if animation should be used, false otherwise
         * \see get_time(), set_time()
         */
        void set_animate(bool x) {
            animate_ = x;
        }
        
        /**
         * \brief Gets the time of the animation.
         * \details Used if animate mode is set.
         * \return the time of the animation, betwe 0.0 (initial) 
         *  and 1.0 (final)
         * \see get_animate(), set_animate()
         */
        double get_time() const {
            return time_;
        }

        /**
         * \brief Gets the time of the animation.
         * \details Used if animate mode is set.
         * \param[in] x the time of the animation, betwe 0.0 (initial) 
         *  and 1.0 (final)
         * \see get_animate(), set_animate()
         */
        void set_time(double x) {
            time_ = x;
            if(time_ < 0.0) {
                time_ = 0.0;
            } else if(time_ > 1.0) {
                time_ = 1.0;
            }
        }
        
        /**
         * \brief Gets the cell visibility flag.
         * \details It is possible to specify cell visibility
         *  flags for each individual cell type.
         * \param[in] type one of MESH_TET, MESH_HEX, MESH_PRISM, MESH_PYRAMID
         * \retval true if the cells of \p type should be displayed
         * \retval false otherwise
         */
        bool get_draw_cells(MeshCellType type) const {
            return draw_cells_[type];
        }


        /**
         * \brief Sets the cell visibility flag.
         * \details It is possible to specify cell visibility
         *  flags for each individual cell type.
         * \param[in] type one of MESH_TET, MESH_HEX, MESH_PRISM, MESH_PYRAMID
         * \param[in] x true if mesh cells of type \p type should be displayed,
         *  false otherwise.
         */
        void set_draw_cells(MeshCellType type, bool x) {
            draw_cells_[type] = x;
        }

        /**
         * \brief Sets the points color
         * \details Specifies the color used to display points
         * \param[in] r , g , b , a the components of the points color,
         *  in (0.0 .. 1.0)
         * \see draw_points()
         */
        void set_points_color(float r, float g, float b, float a=1.0f) {
            points_color_[0] = r;
            points_color_[1] = g;
            points_color_[2] = b;
	    points_color_[3] = a;
        }

        /**
         * \brief Gets the points color
         * \param[out] r , g , b , a the components of the points color,
         *  in (0.0 .. 1.0)
         * \see draw_points()
         */
        void get_points_color(float& r, float& g, float& b, float& a) const {
            r = points_color_[0];
            g = points_color_[1];
            b = points_color_[2];
	    a = points_color_[3];
        }

        /**
         * \brief Sets the point size
         * \param[in] x the point size (minimum 1)
         * \see draw_points()
         */
        void set_points_size(float x) {
            points_size_ = x;
        }

        /**
         * \brief Gets the point size
         * \return the point size
         * \see draw_points()
         */
        float get_points_size() const {
            return points_size_;
        }
        
        /**
         * \brief Sets the mesh color
         * \details Specifies the mesh color to be used if 
         *  mesh edges should be displayed.
         * \param[in] r , g , b , a the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void set_mesh_color(float r, float g, float b, float a=1.0f) {
            mesh_color_[0] = r;
            mesh_color_[1] = g;
            mesh_color_[2] = b;
	    mesh_color_[3] = a;
        }

        /**
         * \brief Gets the mesh color
         * \param[out] r , g , b , a the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void get_mesh_color(float& r, float& g, float& b, float& a) const {
            r = mesh_color_[0];
            g = mesh_color_[1];
            b = mesh_color_[2];
	    a = mesh_color_[3];
        }
        
        /**
         * \brief Sets the surface color
         * \details Specifies the color used to display the
         *  surfacic part of the mesh. It specifies the color 
         *  of both frontfacing and backfacing faces.
         * \param[in] r , g , b , a the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see draw_surface(), set_backface_surface_color()
         */
        void set_surface_color(float r, float g, float b, float a=1.0f) {
            surface_color_[0] = r;
            surface_color_[1] = g;
            surface_color_[2] = b;
	    surface_color_[3] = a;
            backface_surface_color_[0] = r;
            backface_surface_color_[1] = g;
            backface_surface_color_[2] = b;
            backface_surface_color_[3] = a;	    
        }

        /**
         * \brief Gets the surface color
         * \param[out] r , g , b , a the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see draw_surface()
         */
        void get_surface_color(float& r, float& g, float& b, float& a) const {
            r = surface_color_[0];
            g = surface_color_[1];
            b = surface_color_[2];
	    a = surface_color_[3];
        }
        
        /**
         * \brief Sets the surface color for backfacing faces.
         * \details Specifies the color used to display the
         *  backfaces of the surfacic part of the mesh. 
         * \param[in] r , g , b , a the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void set_backface_surface_color(
	    float r, float g, float b, float a=1.0f
	) {
            backface_surface_color_[0] = r;
            backface_surface_color_[1] = g;
            backface_surface_color_[2] = b;
	    backface_surface_color_[3] = a;
        }

        /**
         * \brief Sets the color used to display mesh cells.
         * \param[in] r , g , b , a the components of the cells color,
         *  in (0.0 .. 1.0)
         * \see set_cells_colors_by_type(), draw_volume()
         */
        void set_cells_color(float r, float g, float b, float a=1.0f) {
            for(index_t i=0; i<MESH_NB_CELL_TYPES; ++i) {
                cells_color_[i][0] = r;
                cells_color_[i][1] = g;
                cells_color_[i][2] = b;
		cells_color_[i][3] = a;
            }
        }

        /**
         * \brief Gets the cells color
         * \param[out] r , g , b , a the components of the cells color,
         *  in (0.0 .. 1.0)
         * \see set_cells_colors_by_type(), draw_volume()
         */
        void get_cells_color(float& r, float& g, float& b, float& a) const {
            r = cells_color_[0][0];
            g = cells_color_[0][1];
            b = cells_color_[0][2];
	    a = cells_color_[0][3];
        }

        /**
         * \brief Sets a different color for each mesh cell type
         * \details it uses the following colors:
         *  - tets: red
         *  - hexes: white
         *  - prisms: green
         *  - pyramids: blue
         * \see set_cells_color(), draw_volume()
         */
        void set_cells_colors_by_type() {
            cells_colors_by_type_ = true;
            set_cells_color(MESH_TET,       1.0f, 0.0f, 0.0f);
            set_cells_color(MESH_HEX,       0.9f, 0.9f, 0.9f);
            set_cells_color(MESH_PRISM,     0.0f, 1.0f, 0.0f);
            set_cells_color(MESH_PYRAMID,   0.0f, 0.0f, 1.0f);
            set_cells_color(MESH_CONNECTOR, 1.0f, 0.8f, 0.0f);            
        }

        // Bind an int attribute value to a color
        // Used to color a mesh according to an int attribute
        // Active after calling set_facets_colors_by_int_attribute()
        void bind_int_attribute_value_to_color(index_t attribute_value, const float* rgba) {
            value_to_color_[attribute_value] = rgba;
        }

        // Ask to color the mesh according to the facet attribute 'name'
        // Return true if ok ('name' is a valid int facet attribute)
        bool set_facets_colors_by_int_attribute(const std::string& name) {
            // based on set_scalar_attribute()

            if( (facets_colors_by_int_attribute_ == true) &&
                (attribute_name_ == name) &&
                (int_attribute_.is_bound()) ) {
                // already set, nothing to do
                return true;
            }

            if(int_attribute_.is_bound()) int_attribute_.unbind();

            if(!int_attribute_.bind_if_is_defined(mesh_->facets.attributes(),name)) { // try to bind int_attribute_ to the 'name' facet attribute of mesh_
                fmt::println(Logger::err("MeshGfx"),"No {} facet attribute",name); Logger::out("MeshGfx") << std::flush;
                unset_facets_colors_by_int_attribute();
                return false;
            }

            geo_assert(int_attribute_.is_bound());

            facets_colors_by_int_attribute_ = true;
            attribute_name_ = name;
            return true;
        }

        // Cancel settings of set_facets_colors_by_int_attribute()
        void unset_facets_colors_by_int_attribute() {
            facets_colors_by_int_attribute_ = false;
            attribute_name_.clear();
        }
        
        /**
         * \brief Gets the lighing flag
         * \retval true if lighting should be used
         * \retval false otherwise
         */
        bool get_lighting() const {
            return lighting_;
        }

        /**
         * \brief Sets the lighting flag
         * \param[in] x true if lighting should be used, false
         *  otherwise.
         */
        void set_lighting(bool x) {
            lighting_ = x;
        }

        /**
         * \brief Sets the mesh
         * \param[in] M a pointer to the mesh that should be
         *  displayed.
         */
        void set_mesh(const Mesh* M);
        
        /**
         * \brief Gets the mesh
         * \return a pointer to the mesh that will be displayed.
         */
        const Mesh* mesh() const {
            return mesh_;
        }

        /**
         * \brief Sets picking mode.
         * \details If picking mode is MESH_NONE, then normal drawing
         *  is activated, else the color is replaced with the index of
         *  the elements.
         * \param[in] what a bitwise or ('|') combination of 
         *  MESH_VERTICES, MESH_EDGES, MESH_FACETS, MESH_CELLS,
         *  or MESH_NONE if picking mode should be deactivated
         * \note Picking mode is currently only implemented with
         *  GLSL support at least v1.5, and only works for points,
         *  segments, pure triangle meshes and pure tetrahedral meshes
         *  (facet picking for polygonal meshes and cell picking for 
         *   hybrid meshes are not implemented yet).
         */
        void set_picking_mode(MeshElementsFlags what) {
            picking_mode_ = what;
        }

        /**
         * \brief Gets the current picking mode.
         * \return a bitwise or ('|') combination of 
         *  MESH_VERTICES, MESH_EDGES, MESH_FACETS, MESH_CELLS,
         *  or MESH_NONE if picking mode is deactivated
         */
        MeshElementsFlags get_picking_mode() const {
            return picking_mode_;
        }
        

        /**
         * \brief Sets the object-wide picking id.
         * \details When the object-wide picking id is set,
         *  it is used to draw all primitives in the picking buffer.
         * \param[in] id the object-wide picking id, or index_t(-1) if
         *  normal per-primitive picking should be used.
         */
        void set_object_picking_id(index_t id) {
            object_picking_id_ = id;
        }

        /**
         * \brief Gets the object-wide picking id.
         * \return the object-wide picking id, or index_t(-1) if
         *  normal per-primitive picking should be used.
         */
        index_t get_object_picking_id() const {
            return object_picking_id_;
        }

        /**
         * \brief Sets the vertices selection.
         * \details If set, the vertices selection is used to determine
         *  which vertices should be drawn by draw_vertices().
         * \param[in] name the name of a Property<bool> bound to
         *  the vertices of the mesh.
         * \see draw_vertices()
         */
        void set_vertices_selection(const std::string& name) {
            vertices_selection_ = name;
        }

        /**
         * \brief Gets the vertices selection.
         * \return the name of the vertices selection, or an
         *  empty string if no vertices selection is set.
         */
        const std::string& get_vertices_selection() const {
            return vertices_selection_;
        }

        /**
         * \brief Sets the parameters for displaying a
         *  scalar attribute using texture mapping.
         * \param[in] subelements one of MESH_VERTICES, MESH_FACETS,
         *  MESH_FACET_CORNERS, MESH_CELLS, MESH_CELL_CORNERS,
         *  MESH_CELL_FACETS
         * \param[in] name name of the attribute with an optional index,
         *   for instance, "foobar[5]" refers to the 5th coordinate of
         *   the "foobar" vector attribute.
         * \param[in] attr_min value of the attribute that is bound to 
         *   the leftmost color in the colormap
         * \param[in] attr_max value of the attribute that is bound to
         *   the rightmost color in the colormap
         * \param[in] colormap_texture the texture to be used to display
         *   the attribute colors
	 * \param[in] repeat the number of times the colorramp should be
	 *   repeated within the specified range.
         */
        void set_scalar_attribute(
            MeshElementsFlags subelements,
            const std::string& name,
            double attr_min, double attr_max,
            GLuint colormap_texture,
            index_t repeat = 1
        );

	/**
	 * \brief Sets the parameters for texture mapping.
	 * \param[in] subelements the subelements that have texture
	 *  coordinates.
	 * \param[in] attribute_name the name of the attribute that 
	 *  has the texture coordinates. Can be a 2d or 3d vector attribute.
	 * \param[in] texture the texture.
	 * \param[in] repeat the number of times the texture should be repeated
	 *  in the unit square in texture space.
	 */
	void set_texturing(
	    MeshElementsFlags subelements,
	    const std::string& attribute_name,
	    GLuint texture,
	    index_t texture_dim,
	    index_t repeat = 1
	);

	
        /**
         * \brief Unsets scalar attribute display.
         */
        void unset_scalar_attribute() {
            attribute_subelements_ = MESH_NONE;
            attribute_min_ = 0.0;
            attribute_max_ = 0.0;
            attribute_texture_ = 0;
            attribute_repeat_ = 1;
	    attribute_dim_ = 0;
        }

        void clear_custom_drawings() {
            custom_points_.clear();
            custom_edges_.clear();
        }

        /**
         * \brief Sets primitive filtering
         * \details Primitive filtering deactivates display of some
         *  primitives based on their id and an attribute. The attribute
         *  is of type Numeric::uint8. If the attribute value is zero, then
         *  the corresponding primitive is not displayed.
         * \param[in] subelements one of 
         *  MESH_VERTICES, MESH_FACETS, MESH_CELLS, MESH_ALL_ELEMENTS
         * \param[in] name the name of the attribute with the
         *  filter. Default is "filter"
         */
        void set_filter(
            MeshElementsFlags subelements,
            const std::string& name = "filter"
        );

        /**
         * \brief Unset primitive filtering
         */
        void unset_filters();

        void set_custom_points_color(float r, float g, float b, float a) {
            custom_points_color_[0] = r;
            custom_points_color_[1] = g;
            custom_points_color_[2] = b;
            custom_points_color_[3] = a;
        }

        void add_custom_point(double x, double y, double z) {
            custom_points_.push_back({x,y,z});
        }

        void add_custom_edge(float r, float g, float b, float a, double x1, double y1, double z1, double x2, double y2, double z2) {
            custom_edges_[{r,g,b,a}].push_back(
                std::make_pair<vec3d,vec3d>({x1,y1,z1},{x2,y2,z2})
            );
        }
        
    protected:
        
        void draw_vertices_array();
        void draw_vertices_immediate_plain();
        void draw_vertices_immediate_attrib();
        void draw_vertices_selection();
        
        void draw_edges_array();
        void draw_edges_immediate_plain();
        void draw_edges_immediate_attrib();

        void draw_triangles();
        void draw_triangles_array();
        void draw_triangles_immediate_plain();
        void draw_triangles_immediate_attrib();
        void draw_triangles_immediate_by_int_attrib();

        void draw_quads();
        void draw_quads_array();
        void draw_quads_immediate_plain();
        void draw_quads_immediate_attrib();
        
        void draw_triangles_and_quads();
        void draw_triangles_and_quads_array();
        void draw_triangles_and_quads_immediate_plain();
        void draw_triangles_and_quads_immediate_attrib();

        void draw_polygons();
        void draw_polygons_plain();
        void draw_polygons_attrib();

        void draw_tets();
        void draw_tets_array();
        void draw_tets_immediate_plain();
        void draw_tets_immediate_attrib();
        
        void draw_hybrid();
        void draw_hybrid_array();
        void draw_hybrid_immediate_plain();
        void draw_hybrid_immediate_attrib();
        
        void set_cells_color(MeshCellType type, float r, float g, float b) {
            cells_color_[type][0] = r;
            cells_color_[type][1] = g;
            cells_color_[type][2] = b;            
        }

	void draw_attribute_as_tex_coord(index_t element) {
	    if(picking_mode_ == MESH_NONE) {
		switch(attribute_dim_) {
		    case 1:
			glupPrivateTexCoord1d(scalar_attribute_[element]);
			break;
		    case 2:
			glupPrivateTexCoord2d(
			    tex_coord_attribute_[0][element],
			    tex_coord_attribute_[1][element]
			);
			break;
		    case 3:
			glupPrivateTexCoord3d(
			    tex_coord_attribute_[0][element],
			    tex_coord_attribute_[1][element],
			    tex_coord_attribute_[2][element]
			);
			break;
		}
	    }
	}
	
        void draw_vertex_with_attribute(index_t vertex) {
            if(attribute_subelements_ == MESH_VERTICES) {
		draw_attribute_as_tex_coord(vertex);
            }
            draw_vertex(vertex);
        }

        void draw_surface_vertex_with_attribute(
            index_t vertex, index_t facet, index_t corner
        ) {
	    if(attribute_subelements_ == MESH_VERTICES) {
		draw_attribute_as_tex_coord(vertex);
	    } else if(attribute_subelements_ == MESH_FACETS) {
		draw_attribute_as_tex_coord(facet);
	    } else if(attribute_subelements_ == MESH_FACET_CORNERS) {
		draw_attribute_as_tex_coord(corner);
	    }
            draw_vertex(vertex);
        }

        void draw_volume_vertex_with_attribute(
            index_t vertex, index_t cell, index_t cell_corner
        ) {
	    if(attribute_subelements_ == MESH_VERTICES) {
		draw_attribute_as_tex_coord(vertex);
	    } else if(attribute_subelements_ == MESH_CELLS) {
		draw_attribute_as_tex_coord(cell);
	    } else if(attribute_subelements_ == MESH_CELL_CORNERS) {
		draw_attribute_as_tex_coord(cell_corner);
	    }
            draw_vertex(vertex);
        }

        void draw_vertex(index_t v) {
            if(do_animation_) {
                if(mesh_->vertices.single_precision()) {
                    const GLUPfloat* p =
                        mesh_->vertices.single_precision_point_ptr(v);
                    float t = float(time_);
                    float s = 1.0f - float(time_);
                    glupPrivateVertex3f(
                        s*p[0] + t*p[3],
                        s*p[1] + t*p[4],
                        s*p[2] + t*p[5]
                    );
                } else {
                    const GLUPdouble* p = mesh_->vertices.point_ptr(v);
                    double s = 1.0 - time_;
                    glupPrivateVertex3d(
                        s*p[0] + time_*p[3],
                        s*p[1] + time_*p[4],
                        s*p[2] + time_*p[5]
                    );
                }
            } else {
                if(mesh_->vertices.single_precision()) {
		    if(mesh_->vertices.dimension() < 3) {
			glupPrivateVertex2fv(
			    mesh_->vertices.single_precision_point_ptr(v)
			);
		    } else {
			glupPrivateVertex3fv(
			    mesh_->vertices.single_precision_point_ptr(v)
			);
		    }
                } else {
		    if(mesh_->vertices.dimension() < 3) {
			glupPrivateVertex2dv(
			    mesh_->vertices.point_ptr(v)
			);
		    } else {
			glupPrivateVertex3dv(
			    mesh_->vertices.point_ptr(v)
			);
		    }
                }
            }
        }

        void draw_surface_mesh_with_lines();
        
        /**
         * \brief Sets GLUP drawing parameters.
         */
        void set_GLUP_parameters();

        /**
         * \brief Sets GLUP picking mode for drawing primitives
         *  of a given type, or deactivates GLUP picking if CustomMeshGfx picking mode
         *  is deactivated.
         * \param[in] what one of 
	 *  MESH_VERTICES, MESH_EDGES, MESH_FACETS, MESH_CELLS.
         */
        void set_GLUP_picking(MeshElementsFlags what);


        /**
         * \brief Encodes an id as the current vertex color.
         * \details This is required for drawing polygons, that
         *  cannot use standard GLUP primitives and picking.
         */
        void set_GLUP_vertex_color_from_picking_id(index_t id);


        /**
         * \brief Updates the Vertex Buffer Objects and Vertex Array
         *  Objects.
         * \details The buffer objects are updated if buffer_objects_dirty_
         *  is set, then buffer_objects_dirty_ is reset. If 
         *  buffer_objects_dirty_ is not set, it checks whether the sizes
         *  of the buffer objects match the size of the mesh arrays.
         */
        void update_buffer_objects_if_needed();


        /**
         * \brief Updates the buffer objects used to display attributes.
         * \details The buffer objects are updated if 
         *  attribute_buffer_objects_dirty_
         *  is set, then attribute_buffer_objects_dirty_ is reset. If 
         *  attribute_buffer_objects_dirty_ is not set, it checks whether 
         *  the sizes of the buffer objects match the size of the mesh arrays.
         */
        void update_attribute_buffer_objects_if_needed();


        /**
         * \brief Binds the attribute buffer object to a Vertex Array
         *  Object.
         * \param[in] VAO one of vertices_VAO_, edges_VAO_, facets_VAO_
         *  or cells_VAO_. If zero, the function does nothing.
         * \pre attribute_.is_bound()
         */
        void bind_attribute_buffer_object(GLuint VAO);


        /**
         * \brief Unbinds the attribute buffer object from a Vertex
         *   Array Object.
         * \param[in] VAO one of vertices_VAO_, edges_VAO_, facets_VAO_
         *  or cells_VAO_. If zero, the function does nothing.
         */
        void unbind_attribute_buffer_object(GLuint VAO);
        
        /**
         * \brief Binds the vertices VBO to the current VAO.
         */
        void bind_vertices_VBO();

        /**
         * \brief Setups drawing for attributes.
         * \details If no attribute is bound, does nothing.
         */
        void begin_attributes();

        /**
         * \brief Deactivates drawing for attributes.
         */
        void end_attributes();


        /**
         * \brief Tests whether array mode can be used
         *  to draw a specified GLUP primitive.
         * \param[in] prim the GLUP primitive
         */
        bool can_use_array_mode(GLUPprimitive prim) const;

        void update_surface_elements();
        
        void update_volume_elements();
        
    protected:
        bool show_mesh_;
        index_t mesh_width_;
        index_t mesh_border_width_;
        double shrink_;
        bool animate_;
        double time_;
        bool draw_cells_[MESH_NB_CELL_TYPES];
        float points_size_;
        
        float points_color_[4];
        float mesh_color_[4];
        float surface_color_[4];
        float backface_surface_color_[4];
        float cells_color_[MESH_NB_CELL_TYPES][4];
        bool cells_colors_by_type_;
        bool facets_colors_by_int_attribute_; // new drawing mode

        bool lighting_;

        MeshElementsFlags picking_mode_;
        index_t object_picking_id_;

        std::string vertices_selection_;

        bool do_animation_;
        
        const Mesh* mesh_;
        bool triangles_and_quads_;
        bool quads_;
        bool has_cells_[MESH_NB_CELL_TYPES];
        
        bool buffer_objects_dirty_;
        bool attributes_buffer_objects_dirty_;
	bool long_vector_attribute_;

        GLuint vertices_VAO_;
        GLuint edges_VAO_;
        GLuint facets_VAO_;
        GLuint cells_VAO_;
        
        GLuint vertices_VBO_;
        GLuint edge_indices_VBO_;
        GLuint facet_indices_VBO_;
        GLuint cell_indices_VBO_;
        GLuint vertices_attribute_VBO_;
        
        MeshElementsFlags attribute_subelements_; // used for scalar attributes only
        std::string attribute_name_; // used for scalar attributes & for int attributes (in this case, store the name of the attribute currently in int_attribute_)
	    index_t attribute_dim_; // used for scalar attributes only
        double attribute_min_; // used for scalar attributes only
        double attribute_max_; // used for scalar attributes only
        GLuint attribute_texture_; // used for scalar attributes only
	    index_t attribute_texture_dim_; // used for scalar attributes only
        index_t attribute_repeat_; // used for scalar attributes only
        ReadOnlyScalarAttributeAdapter scalar_attribute_; // used for scalar attributes only
	    ReadOnlyScalarAttributeAdapter tex_coord_attribute_[3]; // used for scalar attributes only
        bool ES_profile_;

        std::map<index_t,const GLUPfloat*> value_to_color_; // used for int attributes only (when facets_colors_by_int_attribute_==true)
        Attribute<index_t> int_attribute_; // used for int attributes only (when facets_colors_by_int_attribute_==true)

        std::vector<std::array<double,3>> custom_points_;
        float custom_points_color_[4];

        std::map<vec4f,std::vector<std::pair<vec3d,vec3d>>,VecngCompare<4,float>> custom_edges_; // edges grouped by color

        /**
         * \brief Filters primitives based on their id and on
         *  an attribute.
         */
        struct Filter {
            Filter();
            ~Filter();

            /**
             * \brief Begins rendering with primitive filtering.
             * \param[in] attributes_manager a reference to the
             *  attributes manager where the property is stored
             * \param[in] hw_primitive_filtering if set, uses
             *  hardware primitive fitering 
             *  (GLUP_PRIMITIVES_FILTERING), else one may use
             *  the function test() instead.
             * \retval true if attribute filtering is active, that is,
             *  if attribute name is not the empty string and if an
             *  the attribute exists
             * \retval false otherwise
             */
            bool begin(
                AttributesManager& attributes_manager,
                bool hw_primitive_filtering=true
            );

            /**
             * \brief Needs to be called after rendering.
             */
            void end();

            /**
             * \brief Deallocates GPU objects.
             */
            void deallocate();

            /**
             * \brief Tests an individual primitive.
             * \retval true if the primitive should be
             *  displayed
             * \retval false otherwise
             */
            bool test(index_t primitive_id) const {
                return (
                    !attribute.is_bound() ||
                    attribute[primitive_id] != 0
                );
            }
            
            std::string attribute_name;
            Attribute<Numeric::uint8> attribute;
            GLuint VBO;
            GLuint texture;
            bool dirty;
        };

        Filter vertices_filter_;
        Filter facets_filter_;
        Filter cells_filter_;

        /**
         * \brief Tests whether hardware primitive filtering is supported.
         * \details On some higher profiles (GLUP150, GLUP440), it is possible
         *  to filter primitives in the shaders. It reduces GPU transfers by
         *  making it possible to reuse the same VBOs even when the list of
         *  objects to display changes.
         */
        bool hw_filtering_supported() const;
        
        /**
         * \brief Generic function to extract element sequences to draw
         * \details Tests the filter if set, and finds all intervals of
         *  primitives to render. For each interval, it calls
         *  glupBasePickingId(GLUPuint64(begin)) so that picking Ids are
         *  correct.
         * \param[in] elements the MeshSubElementsStore to render
         * \param[in] draw the function to be called for each (begin,end) 
         *  sequence of elements to be drawn, where begin is the index of
         *  the first element, and end one position past the index of the
         *  last element.
         */
        void draw_sequences(
            const MeshSubElementsStore& elements,
            std::function<void(index_t, index_t)> draw
        ) {
            // GLUP hardware primitive filtering.
            // Not supported in all profiles.
            // If supported, primitives are filtered by the GPU,
            // else, they are filtered here explicitly.
            const bool hw_filtering = hw_filtering_supported();
            
            Filter* filter = nullptr;
            if(&elements == &mesh_->vertices) {
                filter = &vertices_filter_;
            } else if(&elements == &mesh_->facets) {
                filter = &facets_filter_;
            } else if(&elements == &mesh_->cells) {
                filter = &cells_filter_;
            }
            geo_assert(filter != nullptr) ;
            
            filter->begin(elements.attributes(), hw_filtering);
            index_t e = 0;
            while(e < elements.nb()) {
                while(e<elements.nb() && !(hw_filtering||filter->test(e))) {
                    ++e;
                }
                if(e<elements.nb()) {
                    index_t begin = e;
                    while(e<elements.nb() && (hw_filtering||filter->test(e))) {
                        ++e;
                    }
                    index_t end = e;
                    glupBasePickingId(GLUPuint64(begin));
                    draw(begin,end);
                }
            }
            filter->end();
            glupBasePickingId(0);
        }

        /**
         * \brief Generic function to extract element sequences to draw
         * \details Tests the filter if set, and finds all intervals of
         *  primitives to render. For each interval, it calls
         *  glupBasePickingId(GLUPuint64(begin)) so that picking Ids are
         *  correct.
         * \param[in] elements the MeshSubElementsStore to render
         * \param[in] predicate a function that tests which elements should
         *  be drawn or not
         * \param[in] draw the function to be called for each (begin,end) 
         *  sequence of elements to be drawn, where begin is the index of
         *  the first element, and end one position past the index of the
         *  last element.
         */
        void draw_sequences_if(
            const MeshSubElementsStore& elements,
            std::function<bool(index_t)> predicate,
            std::function<void(index_t, index_t)> draw
        ) {
            // GLUP hardware primitive filtering.
            // Not supported in all profiles.
            // If supported, primitives are filtered by the GPU,
            // else, they are filtered here explicitly.
            const bool hw_filtering = hw_filtering_supported();
            
            Filter* filter = nullptr;
            if(&elements == &mesh_->vertices) {
                filter = &vertices_filter_;
            } else if(&elements == &mesh_->facets) {
                filter = &facets_filter_;
            } else if(&elements == &mesh_->cells) {
                filter = &cells_filter_;
            }
            geo_assert(filter != nullptr) ;
            filter->begin(elements.attributes(), hw_filtering);
            
            index_t e = 0;
            while(e < elements.nb()) {
                while(
                    e<elements.nb() && (
                        !(hw_filtering || filter->test(e)) || !predicate(e))
                ) {
                    ++e;
                }
                if(e<elements.nb()) {
                    index_t begin = e;
                    while(
                        e<elements.nb() &&
                        (hw_filtering || filter->test(e)) &&
                        predicate(e)
                    ) {
                        ++e;
                    }
                    index_t end = e;
                    glupBasePickingId(GLUPuint64(begin));
                    draw(begin,end);
                }
            }
            filter->end();
            glupBasePickingId(0);
        }
        
    };
