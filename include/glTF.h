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

// Image created with GIMP, then converted to base64
// https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/Data_URLs

// BMP (bitmap) format
// base64 = "Qk2eAAAAAAAAAIoAAAB8AAAABgAAAAEAAAABABgAAAAAABQAAAAjLgAAIy4AAAAAAAAAAAAAAAD/AAD/AAD/AAAAAAAAAEJHUnMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAD/AACZ////mZmZ/wAAmQAAAAA="
// 212 chars -> 158 bytes + 2 of padding
const uint32_t labeling_texture_image_BMP[20] = {
    uint32_t(0x424d9e0000000000),
    uint32_t(0x00008a0000007c00),
    uint32_t(0x0000060000000100),
    uint32_t(0x0000010018000000),
    uint32_t(0x000014000000232e),
    uint32_t(0x0000232e00000000),
    uint32_t(0x0000000000000000),
    uint32_t(0xff0000ff0000ff00),
    uint32_t(0x0000000000004247),
    uint32_t(0x5273000000000000),
    uint32_t(0x0000000000000000),
    uint32_t(0x0000000000000000),
    uint32_t(0x0000000000000000),
    uint32_t(0x0000000000000000),
    uint32_t(0x0000000000000000),
    uint32_t(0x0000020000000000),
    uint32_t(0x0000000000000000),
    uint32_t(0x00000000ff000099),
    uint32_t(0xffffff999999ff00),
    uint32_t(0x0099000000000000)
};

// PNG format
// base64:
const std::string labeling_texture_image_PNG = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAYAAAABCAIAAAByq0inAAABhGlDQ1BJQ0MgcHJvZmlsZQAAKJF9kT1Iw1AUhU9TRZEWETuIOGSoTlZERRy1CkWoEGqFVh1MXvoHTRqSFBdHwbXg4M9i1cHFWVcHV0EQ/AFxdnBSdJES70sKLWK88Hgf591zeO8+QKiXmWZ1jAOabpupRFzMZFfFrlcI6EMYAYzJzDLmJCkJ3/q6pz6quxjP8u/7s8JqzmJAQCSeZYZpE28QT2/aBud94ggryirxOfGoSRckfuS64vEb54LLAs+MmOnUPHGEWCy0sdLGrGhqxFPEUVXTKV/IeKxy3uKslauseU/+wlBOX1nmOq0hJLCIJUgQoaCKEsqwEaNdJ8VCis7jPv5B1y+RSyFXCYwcC6hAg+z6wf/g92yt/OSElxSKA50vjvMxDHTtAo2a43wfO07jBAg+A1d6y1+pAzOfpNdaWvQI6N0GLq5bmrIHXO4AA0+GbMquFKQl5PPA+xl9UxbovwV61ry5Nc9x+gCkaVbJG+DgEBgpUPa6z7u72+f2b09zfj85enKQ1aMGywAAAAlwSFlzAAAuIwAALiMBeKU/dgAAAAd0SU1FB+gCEQ0GMNz9VKwAAAAZdEVYdENvbW1lbnQAQ3JlYXRlZCB3aXRoIEdJTVBXgQ4XAAAAGUlEQVQI1wXBAQEAAACCIP4/cG2BEdsqRgdOxAf5SBcyAQAAAABJRU5ErkJggg==";
// 748 chars -> 559 bytes + 1 of padding
// 89504e470d0a1a0a
// 0000000d49484452
// 0000000600000001
// 080200000072ab48
// a700000184694343
// 504943432070726f
// 66696c6500002891
// 7d913d48c3501485
// 4f53459116113b88
// 3864a84e5644451c
// b50a45a8106a8556
// 1d4c5efa074d1a92
// 141747c1b5e0e0cf
// 62d5c1c559570757
// 4110fc0171767052
// 749112ef4b0a2d62
// bcf0781fe7dd7378
// ef3e40a897996675
// 8c039a6e9ba9445c
// cc6457c5ae5708e8
// 4318018cc9cc32e6
// 242909dffabaa73e
// aabb18cff2effbb3
// c26ace624040249e
// 658669136f104f6f
// da06e77de2082bca
// 2af139f1a8491724
// 7ee4bae2f11be782
// cb02cf8c98e9d43c
// 7184582cb4b1d2c6
// ac686ac453c45155
// d3295fc878ac72de
// e2ac95abac794ffe
// c2504e5f59e63aad
// 2124b088254810a1
// a08a12cab011a35d
// 27c5428acee33efe
// 41d72f914b215709
// 8c1c0ba84083ecfa
// c1ffe0f76cadfce4
// 8497148a039d2f8e
// f3310c74ed028d9a
// e37c1f3b4ee30408
// 3e03577acb5fa903
// 339fa4d75a5af408
// e8dd062eae5b9ab2
// 075cee00034f866c
// caae14a425e4f3c0
// fb197d5316e8bf05
// 7ad6bcb935cf71fa
// 00a46956c91be0e0
// 10182950f6bacfbb
// bbdbe7f66f4f737e
// 3f397a7290d5a306
// cb00000009704859
// 7300002e2300002e
// 230178a53f760000
// 000774494d4507e8
// 02110d0630dcfd54
// ac00000019744558
// 74436f6d6d656e74
// 0043726561746564
// 2077697468204749
// 4d5057810e170000
// 00194944415408d7
// 05c1010100000082
// 20fe3f706d8111db
// 2a46074ec407f948
// 1732010000000049
// 454e44ae42608200

//////////////
// Functions
//////////////

inline index_t* facet_vertex_index_ptr(GEO::Mesh& M, index_t f, index_t lv) {
    return M.facet_corners.vertex_index_ptr(M.facets.corner(f,lv));
}

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M, bool with_wireframe);

void write_glTF__labeled_triangle_mesh(std::string filename, GEO::Mesh& M, const char* attribute_name);