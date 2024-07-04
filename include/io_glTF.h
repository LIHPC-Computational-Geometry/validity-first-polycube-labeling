#pragma once

#include "tiny_gltf.h"

#include <geogram/basic/command_line_args.h>    // for import_arg_group()
#include <geogram/basic/geometry.h>             // for vec3
#include <geogram/basic/logger.h>               // for Logger
#include <geogram/basic/assert.h>               // for geo_assert*
#include <geogram/mesh/mesh.h>                  // for Mesh
#include <geogram/mesh/mesh_io.h>               // for mesh_save()
#include <geogram/mesh/mesh_halfedges.h>        // for Halfedge
#include <geogram/mesh/mesh_geometry.h>         // for get_bbox()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "geometry_halfedges.h"     // for CustomMeshHalfedges
#include "geometry.h"               // for AdjacentFacetOfVertex

#ifdef GARGANTUA
    #define ASSERT_GARGANTUA_OFF (fmt::println("Cannot export to glTF in GARGANTUA mode, indices must have 32 bits"); geo_assert_not_reached;)
#else
    #define ASSERT_GARGANTUA_OFF geo_assert(sizeof(index_t)==4); // 4 bytes = 32 bits
#endif

using GEO::index_t; // to use the FOR() macro of Geogram

/////////////////////////////
// Labeling colors as image
/////////////////////////////

// The minimal texture image would have 6 pixels, one for each label in +/-{X,Y,Z}
// But texture images should have power-of-two dimensions -> 8x2 image
//
//                        R    G    B
// 0 = +X ->       red : 1.0, 0.0, 0.0
// 1 = -X ->  dark red : 0.4, 0.0, 0.0
// 2 = +Y ->     white : 1.0, 1.0, 1.0
// 3 = -Y ->      grey : 0.4, 0.4, 0.4
// 4 = +Z ->      blue : 0.0, 0.0, 1.0
// 5 = -Z -> dark blue : 0.0, 0.0, 0.4
// 
// Texture coordinates :
//
//   0,0                            0.5,0                            1,0
//    +-------+-------+-------+-------+-------+-------+-------+-------+
//    |               |       |       |       |       |               |
//    |               |       |       |       |       |               |
//    |               |       |       |       |       |               |
//    +       +  red  + dark  + white + grey  + blue  + dark  +       +
//    |               |  red  |       |       |       |  blue         |
//    |               |       |       |       |       |               |
//    |               |       |       |       |       |               |
//    +-------+-------+-------+-------+-------+-------+-------+-------+
//   0,1                            0.5,1                            1,1
//
const vec2f LABELING_TO_TEXTURE_COORDINATES[6] = {{0.188f, 0.5f},  // = 0/8 + 3/16, 1/2 -> texture coordinate of +X
                                                  {0.312f, 0.5f},  // = 1/8 + 3/16, 1/2 -> texture coordinate of -X
                                                  {0.438f, 0.5f},  // = 2/8 + 3/16, 1/2 -> texture coordinate of Y
                                                  {0.563f, 0.5f},  // = 3/8 + 3/16, 1/2 -> texture coordinate of -Y
                                                  {0.688f, 0.5f},  // = 4/8 + 3/16, 1/2 -> texture coordinate of +Z
                                                  {0.813f, 0.5f}}; // = 5/8 + 3/16, 1/2 -> texture coordinate of -Z

// Image created with GIMP, saved to PNG, then converted to base64
// https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/Data_URLs
const std::string labeling_texture_image_PNG = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAgAAAACCAIAAADq9gq6AAAAJUlEQVQI13XBAREAIAgEsJ1GsqABIcfHwARue7iEmknSXYSzfDysLAtr+BxIUwAAAABJRU5ErkJggg==";

/////////////////////////////
// Parula colormap as image
/////////////////////////////

// it's ext/geogram/src/lib/geogram_gfx/gui/colormaps/parula.xpm
// opened with GIMP
// flipped (so that coordinate 0 is yellow and 1 is blue)
// resized to 32x2 (power-of-two dimensions)
// exported to PNG
// then converted to base64
const std::string parula_texture_image_PNG = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAACCAYAAAA5Ht7JAAAAcElEQVQY08XBSRKCMBBA0U8jIbboyp33v5pDIVSoaMyAxGP4XhPcuW5eab2lCxYJSooDS1UWemaxLGKYG8OE5dlYHpty58j4PTCtA64or3VPLkrNPRI7TBTMG4g7SurIsSV9Wk6+cHEj6m+IvyL82Q+DHTH9UiYXWAAAAABJRU5ErkJggg==";

//////////////////////////////////////
// Geogram vs glTF mesh vertices map
//////////////////////////////////////

template <typename T>
struct glTF_vertex {
    index_t Geogram_vertex;
    T value;
};

//////////////
// Functions
//////////////

inline index_t* facet_vertex_index_ptr(GEO::Mesh& M, index_t f, index_t lv) {
    return M.facet_corners.vertex_index_ptr(M.facets.corner(f,lv));
}

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M, bool with_wireframe);

void write_glTF__labeled_triangle_mesh(std::string filename, GEO::Mesh& M, const char* attribute_name, std::vector<std::vector<AdjacentFacetOfVertex>> per_vertex_adj_facets);

void write_glTF__labeled_triangle_mesh_with_polycube_animation(std::string filename, GEO::Mesh& M, GEO::Mesh& polycube, const char* attribute_name, const std::vector<std::vector<AdjacentFacetOfVertex>>& per_vertex_adj_facets);

void write_glTF__hex_mesh_surface(std::string filename, const GEO::Mesh& hex_mesh);