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

#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges
#include "geometry.h"               // for AdjacentFacetOfVertex

#ifdef GARGANTUA
    #define ASSERT_GARGANTUA_OFF (fmt::println("Cannot export to glTF in GARGANTUA mode, indices must have 32 bits"); geo_assert_not_reached;)
#else
    #define ASSERT_GARGANTUA_OFF geo_assert(sizeof(index_t)==4); // 4 bytes = 32 bits
#endif

using GEO::index_t; // to use the FOR() macro of Geogram

/////////////////////
// Labeling texture
/////////////////////

// The minimal texture image would have 6 pixels, one for each label in +/-{X,Y,Z}
// But texture images should have power-of-two dimensions -> 8x2 image
//
//                        R    G    B
// 0 = +X ->       red : 1.0, 0.0, 0.0
// 1 = -X ->  dark red : 0.6, 0.0, 0.0
// 2 = +Y ->     white : 1.0, 1.0, 1.0
// 3 = -Y ->      grey : 0.6, 0.6, 0.6
// 4 = +Z ->      blue : 0.0, 0.0, 1.0
// 5 = -Z -> dark blue : 0.0, 0.0, 0.6
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
const std::string labeling_texture_image_PNG = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAgAAAACCAIAAADq9gq6AAAAH0lEQVQI12P4z8Dwn4FhJgPD////Z86cCePNZGLAAQBunQmThyp5vwAAAABJRU5ErkJggg==";

//////////////////////////////////////
// Geogram vs glTF mesh vertices map
//////////////////////////////////////

struct glTF_vertex {
    index_t Geogram_vertex;
    index_t label;
};

//////////////
// Functions
//////////////

inline index_t* facet_vertex_index_ptr(GEO::Mesh& M, index_t f, index_t lv) {
    return M.facet_corners.vertex_index_ptr(M.facets.corner(f,lv));
}

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M, bool with_wireframe);

void write_glTF__labeled_triangle_mesh(std::string filename, GEO::Mesh& M, const char* attribute_name, const std::vector<std::vector<AdjacentFacetOfVertex>>& per_vertex_adj_facets);