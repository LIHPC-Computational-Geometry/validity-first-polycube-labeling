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
#include <set>
#include <vector>

#include "glTF.h"
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M) {
    ASSERT_GARGANTUA_OFF;
    
    //////////////
    // Get edges
    //////////////

    std::set<std::set<index_t>> edges_set;
    GEO::CustomMeshHalfedges mesh_he(M);
    GEO::MeshHalfedges::Halfedge H;
    index_t H_v0 = index_t(-1);
    index_t H_v1 = index_t(-1);
    FOR(f,M.facets.nb()) { // for each facet
        H.facet = f;
        FOR(lv,3) { // for each local vertex of f
            H.corner = M.facets.corner(f,lv);
            H_v0 = halfedge_vertex_index_from(M,H);
            H_v1 = halfedge_vertex_index_to(M,H);
            edges_set.insert({H_v0,H_v1});
        }
    }
    // we need to have an index associated to each edge
    // -> convert to vector
    std::vector<std::set<index_t>> edges(edges_set.begin(),edges_set.end());

    /////////////////////////
    // Compute bounding box
    /////////////////////////

    // /!\ warning: vertices coordinates must be double precision
    std::vector<double> bounding_box_min(3,0.0), bounding_box_max(3,0.0);
    GEO::get_bbox(M, bounding_box_min.data(), bounding_box_max.data());

    //////////////////////////////////////////
    // Convert vertices coordinates to float
    //////////////////////////////////////////

    M.vertices.set_single_precision();

    ////////////////////////
    // Create a glTF model
    ////////////////////////

    tinygltf::Model m;

    ///////////////////////
    // Create 2 materials
    ///////////////////////

    m.materials.resize(2);
    const size_t MATERIAL_0_GREY = 0;
    const size_t MATERIAL_1_BLACK = 1;
    tinygltf::Material& material_0_grey = m.materials[MATERIAL_0_GREY];
    tinygltf::Material& material_1_black = m.materials[MATERIAL_1_BLACK];

    material_0_grey.pbrMetallicRoughness.baseColorFactor = {0.8f, 0.8f, 0.8f, 1.0f};
    material_0_grey.doubleSided = true;

    material_1_black.pbrMetallicRoughness.baseColorFactor = {0.0f, 0.0f, 0.0f, 1.0f};
    material_1_black.doubleSided = true;

    ////////////////////
    // Create a buffer
    ////////////////////

    m.buffers.resize(1);
    const size_t BUFFER_0 = 0;
    tinygltf::Buffer& buffer_0 = m.buffers[BUFFER_0];

    // resize buffer to 4 bytes * 3 floating points * #vertices
    //                + 4 bytes * 3 indices * #triangles
    //                + 4 bytes * 2 indices * #edges
    size_t buffer_coordinates_start     = 0;                                                // byte index
    size_t buffer_coordinates_length    = 4*3*M.vertices.nb();                              // number of bytes
    size_t buffer_triangles_start       = buffer_coordinates_length;                        // byte index
    size_t buffer_triangles_length      = 4*3*M.facets.nb();                                // number of bytes
    size_t buffer_edges_start           = buffer_triangles_start+buffer_triangles_length;   // byte index
    size_t buffer_edges_length          = 4*2*edges.size();                                 // number of bytes
    buffer_0.data.resize(buffer_coordinates_length+buffer_triangles_length+buffer_edges_length);

    // write vertices (= 3D coordinates)
    FOR(v,M.vertices.nb()) { // for each vertex index
        memcpy(buffer_0.data.data() + v*12, static_cast<void*>(M.vertices.single_precision_point_ptr(v)), 12); // write x,y,z from float* getter
    }
    // write triangles (= vertex indices)
    FOR(f,M.facets.nb()) { // for each facet (triangle) index
        memcpy(buffer_0.data.data() + buffer_triangles_start + f*12 + 0, facet_vertex_index_ptr(M,f,0),4);
        memcpy(buffer_0.data.data() + buffer_triangles_start + f*12 + 4, facet_vertex_index_ptr(M,f,1),4);
        memcpy(buffer_0.data.data() + buffer_triangles_start + f*12 + 8, facet_vertex_index_ptr(M,f,2),4);
    }
    // write edges (= vertex indices)
    FOR(e,edges.size()) {
        memcpy(buffer_0.data.data() + buffer_edges_start + 8*e + 0, &(*edges[e].begin()),  4);
        memcpy(buffer_0.data.data() + buffer_edges_start + 8*e + 4, &(*edges[e].rbegin()), 4);
    }    

    //////////////////////////
    // Create 3 buffer views
    //////////////////////////

    m.bufferViews.resize(3);
    const size_t BUFFERVIEW_0_VERTICES_COORDINATES = 0;
    const size_t BUFFERVIEW_1_TRIANGLES_VERTICES = 1;
    const size_t BUFFERVIEW_2_EDGES_VERTICES = 2;
    tinygltf::BufferView& bufferview_0_vertices_coordinates = m.bufferViews[BUFFERVIEW_0_VERTICES_COORDINATES];
    tinygltf::BufferView& bufferview_1_triangles_vertices = m.bufferViews[BUFFERVIEW_1_TRIANGLES_VERTICES];
    tinygltf::BufferView& bufferview_2_edges_vertices = m.bufferViews[BUFFERVIEW_2_EDGES_VERTICES];

    // 1st chunk of the buffer 0 : points = 3D coordinates (type ARRAY_BUFFER)
    bufferview_0_vertices_coordinates.buffer = BUFFER_0;
    bufferview_0_vertices_coordinates.byteOffset = buffer_coordinates_start;
    bufferview_0_vertices_coordinates.byteLength = buffer_coordinates_length;
    bufferview_0_vertices_coordinates.target = TINYGLTF_TARGET_ARRAY_BUFFER;

    // 2nd chunk of the bufffer 0 : triangles = vertex indices (type ELEMENT_ARRAY_BUFFER)
    bufferview_1_triangles_vertices.buffer = BUFFER_0;
    bufferview_1_triangles_vertices.byteOffset = buffer_triangles_start;
    bufferview_1_triangles_vertices.byteLength = buffer_triangles_length;
    bufferview_1_triangles_vertices.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    // 3rd chunk of the buffer 0 : edges = vertex indices (type ELEMENT_ARRAY_BUFFER)
    bufferview_2_edges_vertices.buffer = BUFFER_0;
    bufferview_2_edges_vertices.byteOffset = buffer_edges_start;
    bufferview_2_edges_vertices.byteLength = buffer_edges_length;
    bufferview_2_edges_vertices.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    ///////////////////////
    // Create 3 accessors
    ///////////////////////

    m.accessors.resize(3);
    const size_t ACCESSOR_0_VERTICES_COORDINATES = 0;
    const size_t ACCESSOR_1_TRIANGLES_VERTICES = 1;
    const size_t ACCESSOR_2_EDGES_VERTICES = 2;
    tinygltf::Accessor& accessor_0_vertices_coordinates = m.accessors[ACCESSOR_0_VERTICES_COORDINATES];
    tinygltf::Accessor& accessor_1_triangles_vertices = m.accessors[ACCESSOR_1_TRIANGLES_VERTICES];
    tinygltf::Accessor& accessor_2_edges_vertices = m.accessors[ACCESSOR_2_EDGES_VERTICES];

    // layout description of buffer view 0
    accessor_0_vertices_coordinates.bufferView = BUFFERVIEW_0_VERTICES_COORDINATES;
    accessor_0_vertices_coordinates.byteOffset = 0;
    accessor_0_vertices_coordinates.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_0_vertices_coordinates.count = M.vertices.nb(); // how many points
    accessor_0_vertices_coordinates.type = TINYGLTF_TYPE_VEC3;
    accessor_0_vertices_coordinates.maxValues = bounding_box_max; // = [max_x, max_y, max_z]
    accessor_0_vertices_coordinates.minValues = bounding_box_min; // = [min_x, min_y, min_z]

    // layout description of buffer view 1
    accessor_1_triangles_vertices.bufferView = BUFFERVIEW_1_TRIANGLES_VERTICES;
    accessor_1_triangles_vertices.byteOffset = 0;
    accessor_1_triangles_vertices.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT; // 32 bits per index
    accessor_1_triangles_vertices.count = M.facets.nb()*3; // how many indices
    accessor_1_triangles_vertices.type = TINYGLTF_TYPE_SCALAR;
    // range of indices value : [ 0 : #vertices-1 ]
    accessor_1_triangles_vertices.maxValues.push_back(M.vertices.nb()-1); // max value of all indices in the buffer view
    accessor_1_triangles_vertices.minValues.push_back(0);

    // layout description of buffer view 2
    accessor_2_edges_vertices.bufferView = BUFFERVIEW_2_EDGES_VERTICES;
    accessor_2_edges_vertices.byteOffset = 0;
    accessor_2_edges_vertices.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT; // 32 bits per index
    accessor_2_edges_vertices.count = edges.size()*2; // how many indices
    accessor_2_edges_vertices.type = TINYGLTF_TYPE_SCALAR;
    // range of indices value : [ 0 : #vertices-1 ]
    accessor_2_edges_vertices.maxValues.push_back(M.vertices.nb()-1); // max value of all indices in the buffer view
    accessor_2_edges_vertices.minValues.push_back(0);

    ////////////////////////////////////
    // Create a mesh with 2 primitives
    ////////////////////////////////////

    m.meshes.resize(1);
    const size_t MESH_0 = 0;
    tinygltf::Mesh& mesh_0 = m.meshes[MESH_0];

    mesh_0.primitives.resize(2);
    const size_t PRIMITIVE_0_TRIANGLES = 0;
    const size_t PRIMITIVE_1_EDGES = 1;
    tinygltf::Primitive& primitive_0_triangles = mesh_0.primitives[PRIMITIVE_0_TRIANGLES];
    tinygltf::Primitive& primitive_1_edges = mesh_0.primitives[PRIMITIVE_1_EDGES];

    primitive_0_triangles.indices = ACCESSOR_1_TRIANGLES_VERTICES;
    primitive_0_triangles.attributes["POSITION"] = ACCESSOR_0_VERTICES_COORDINATES;
    primitive_0_triangles.material = MATERIAL_0_GREY;
    primitive_0_triangles.mode = TINYGLTF_MODE_TRIANGLES;

    primitive_1_edges.indices = ACCESSOR_2_EDGES_VERTICES;
    primitive_1_edges.attributes["POSITION"] = ACCESSOR_0_VERTICES_COORDINATES;
    primitive_1_edges.material = MATERIAL_1_BLACK;
    primitive_1_edges.mode = TINYGLTF_MODE_LINE;

    //////////////////
    // Create a node
    //////////////////
    
    m.nodes.resize(1);
    const size_t NODE_0 = 0;

    m.nodes[NODE_0].mesh = MESH_0; // link NODE_0 to MESH_0

    ///////////////////
    // Create a scene
    ///////////////////

    m.scenes.resize(1);
    const size_t SCENE_0 = 0;

    m.scenes[SCENE_0].nodes.push_back(NODE_0); // link SCENE_0 to NODE_0

    ///////////////////////
    // Describe the asset
    ///////////////////////

    m.asset.version = "2.0"; // required
    m.asset.generator = "tinygltf";
    
    ///////////////////////
    // Save model to file
    ///////////////////////

    tinygltf::TinyGLTF gltf;
    std::string output_filename = filename + ".glb";
    fmt::println(Logger::out("glTF"),"Writing {}...",output_filename); Logger::out("glTF").flush();
    gltf.WriteGltfSceneToFile(&m, output_filename,
                            true, // embedImages
                            true, // embedBuffers
                            true, // pretty print
                            true); // write binary

}