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

inline index_t* facet_vertex_index_ptr(GEO::Mesh& M, index_t f, index_t lv) {
    return M.facet_corners.vertex_index_ptr(M.facets.corner(f,lv));
}

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M, bool with_wireframe);