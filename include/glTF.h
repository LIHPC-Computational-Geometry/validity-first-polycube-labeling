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

#ifdef GARGANTUA
    #define ASSERT_GARGANTUA_OFF (fmt::println("Cannot export to glTF in GARGANTUA mode, indices must have 32 bits"); geo_assert_not_reached;)
#else
    #define ASSERT_GARGANTUA_OFF geo_assert(sizeof(index_t)==4); // 4 bytes = 32 bits
#endif

using GEO::index_t; // to use the FOR() macro of Geogram

/////////////////////
// Labeling texture
/////////////////////

// An image with 6x1 pixels :
// pixel 0,0 : 1.0, 0.0, 0.0 (red)
// pixel 1,0 : 0.6, 0.0, 0.0 (dark red)
// pixel 2,0 : 1.0, 1.0, 1.0 (white)
// pixel 3,0 : 0.6, 0.6, 0.6 (grey)
// pixel 4,0 : 0.0, 0.0, 1.0 (blue)
// pixel 5,0 : 0.0, 0.0, 0.6 (dark blue)
// 
// Texture coordinates :
//
//   0,0                    0.5,0                    1,0
//    +-------+-------+-------+-------+-------+-------+
//    |       | dark  |       |       |       | dark  |
//    |  red  |       | white | grey  | blue  |       |
//    |       |  red  |       |       |       |  blue |
//    +-------+-------+-------+-------+-------+-------+
//   0,1                    0.5,1                    1,1
//
const vec2f LABELING_TO_TEXTURE_COORDINATES[6] = {{0.084f, 0.5f},  // = 0/6 + 1/12, 1/2 -> texture coordinate of 0 = +X
                                                  {0.250f, 0.5f},  // = 1/6 + 1/12, 1/2 -> texture coordinate of 1 = -X
                                                  {0.417f, 0.5f},  // = 2/6 + 1/12, 1/2 -> texture coordinate of 2 = +Y
                                                  {0.534f, 0.5f},  // = 3/6 + 1/12, 1/2 -> texture coordinate of 3 = -Y
                                                  {0.750f, 0.5f},  // = 4/6 + 1/12, 1/2 -> texture coordinate of 4 = +Z
                                                  {0.917f, 0.5f}}; // = 5/6 + 1/12, 1/2 -> texture coordinate of 5 = -Z

// Image created with GIMP, saved to PNG, then converted to base64
// https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/Data_URLs
const std::string labeling_texture_image_PNG = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAYAAAABCAIAAAByq0inAAABhGlDQ1BJQ0MgcHJvZmlsZQAAKJF9kT1Iw1AUhU9TRZEWETuIOGSoTlZERRy1CkWoEGqFVh1MXvoHTRqSFBdHwbXg4M9i1cHFWVcHV0EQ/AFxdnBSdJES70sKLWK88Hgf591zeO8+QKiXmWZ1jAOabpupRFzMZFfFrlcI6EMYAYzJzDLmJCkJ3/q6pz6quxjP8u/7s8JqzmJAQCSeZYZpE28QT2/aBud94ggryirxOfGoSRckfuS64vEb54LLAs+MmOnUPHGEWCy0sdLGrGhqxFPEUVXTKV/IeKxy3uKslauseU/+wlBOX1nmOq0hJLCIJUgQoaCKEsqwEaNdJ8VCis7jPv5B1y+RSyFXCYwcC6hAg+z6wf/g92yt/OSElxSKA50vjvMxDHTtAo2a43wfO07jBAg+A1d6y1+pAzOfpNdaWvQI6N0GLq5bmrIHXO4AA0+GbMquFKQl5PPA+xl9UxbovwV61ry5Nc9x+gCkaVbJG+DgEBgpUPa6z7u72+f2b09zfj85enKQ1aMGywAAAAlwSFlzAAAuIwAALiMBeKU/dgAAAAd0SU1FB+gCEQ0GMNz9VKwAAAAZdEVYdENvbW1lbnQAQ3JlYXRlZCB3aXRoIEdJTVBXgQ4XAAAAGUlEQVQI1wXBAQEAAACCIP4/cG2BEdsqRgdOxAf5SBcyAQAAAABJRU5ErkJggg==";

//////////////
// Functions
//////////////

inline index_t* facet_vertex_index_ptr(GEO::Mesh& M, index_t f, index_t lv) {
    return M.facet_corners.vertex_index_ptr(M.facets.corner(f,lv));
}

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M, bool with_wireframe);

void write_glTF__labeled_triangle_mesh(std::string filename, GEO::Mesh& M, const char* attribute_name);