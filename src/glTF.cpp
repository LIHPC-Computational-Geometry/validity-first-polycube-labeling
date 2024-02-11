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

void write_glTF__triangle_mesh(std::string filename, GEO::Mesh& M, bool with_wireframe) {
    ASSERT_GARGANTUA_OFF;
    
    ////////////////////////
    // Get wireframe edges
    ////////////////////////

    std::set<std::set<index_t>> wireframe_edges_as_set;
    if(with_wireframe) {
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
                wireframe_edges_as_set.insert({H_v0,H_v1});
            }
        }
    }
    // we need to have an index associated to each wireframe edge
    // -> convert to vector
    std::vector<std::set<index_t>> wireframe_edges_as_vector(wireframe_edges_as_set.begin(),wireframe_edges_as_set.end());

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

    //////////////////////////////////////////////
    // Create 1 material (+1 with the wireframe)
    //////////////////////////////////////////////

    m.materials.resize(1);
    const size_t MATERIAL_0_GREY = 0;
    const size_t MATERIAL_1_BLACK = 1;
    tinygltf::Material& material_0_grey = m.materials[MATERIAL_0_GREY];
    
    material_0_grey.pbrMetallicRoughness.baseColorFactor = {0.8f, 0.8f, 0.8f, 1.0f};
    material_0_grey.doubleSided = true;

    if(with_wireframe) {
        m.materials.resize(2);
        tinygltf::Material& material_1_black = m.materials[MATERIAL_1_BLACK];
        material_1_black.pbrMetallicRoughness.baseColorFactor = {0.0f, 0.0f, 0.0f, 1.0f};
        material_1_black.doubleSided = true;
    }

    ////////////////////
    // Create a buffer
    ////////////////////

    m.buffers.resize(1);
    const size_t BUFFER_0 = 0;
    tinygltf::Buffer& buffer_0 = m.buffers[BUFFER_0];

    // resize buffer to 4 bytes * 3 floating points * #vertices
    //                + 4 bytes * 3 indices * #triangles
    //                + 4 bytes * 2 indices * #edges
    size_t buffer_coordinates_start         = 0;                                                // byte index
    size_t buffer_coordinates_length        = 4*3*M.vertices.nb();                              // number of bytes
    size_t buffer_triangles_start           = buffer_coordinates_length;                        // byte index
    size_t buffer_triangles_length          = 4*3*M.facets.nb();                                // number of bytes
    size_t buffer_wireframe_edges_start     = buffer_triangles_start+buffer_triangles_length;   // byte index
    size_t buffer_wireframe_edges_length    = 4*2*wireframe_edges_as_vector.size();             // number of bytes
    if(with_wireframe) {
        buffer_0.data.resize(buffer_coordinates_length+buffer_triangles_length+buffer_wireframe_edges_length);
    }
    else {
        buffer_0.data.resize(buffer_coordinates_length+buffer_triangles_length);
    }

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
    if(with_wireframe) {
        // write wireframe edges (= vertex indices)
        FOR(e,wireframe_edges_as_vector.size()) {
            memcpy(buffer_0.data.data() + buffer_wireframe_edges_start + 8*e + 0, &(*wireframe_edges_as_vector[e].begin()),  4);
            memcpy(buffer_0.data.data() + buffer_wireframe_edges_start + 8*e + 4, &(*wireframe_edges_as_vector[e].rbegin()), 4);
        }
    }    

    //////////////////////////////////////////////////
    // Create 2 buffer views (+1 with the wireframe)
    //////////////////////////////////////////////////

    m.bufferViews.resize(2);
    const size_t BUFFERVIEW_0_VERTICES_COORDINATES = 0;
    const size_t BUFFERVIEW_1_TRIANGLES_VERTICES = 1;
    const size_t BUFFERVIEW_2_WIREFRAME_EDGES_VERTICES = 2;
    tinygltf::BufferView& bufferview_0_vertices_coordinates = m.bufferViews[BUFFERVIEW_0_VERTICES_COORDINATES];
    tinygltf::BufferView& bufferview_1_triangles_vertices = m.bufferViews[BUFFERVIEW_1_TRIANGLES_VERTICES];
    

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

    if(with_wireframe) {
        m.bufferViews.resize(3);
        tinygltf::BufferView& bufferview_2_wireframe_edges_vertices = m.bufferViews[BUFFERVIEW_2_WIREFRAME_EDGES_VERTICES];
        // 3rd chunk of the buffer 0 : edges = vertex indices (type ELEMENT_ARRAY_BUFFER)
        bufferview_2_wireframe_edges_vertices.buffer = BUFFER_0;
        bufferview_2_wireframe_edges_vertices.byteOffset = buffer_wireframe_edges_start;
        bufferview_2_wireframe_edges_vertices.byteLength = buffer_wireframe_edges_length;
        bufferview_2_wireframe_edges_vertices.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
    }

    ///////////////////////////////////////////////
    // Create 2 accessors (+1 with the wireframe)
    ///////////////////////////////////////////////

    m.accessors.resize(2);
    const size_t ACCESSOR_0_VERTICES_COORDINATES = 0;
    const size_t ACCESSOR_1_TRIANGLES_VERTICES = 1;
    const size_t ACCESSOR_2_WIREFRAME_EDGES_VERTICES = 2;
    tinygltf::Accessor& accessor_0_vertices_coordinates = m.accessors[ACCESSOR_0_VERTICES_COORDINATES];
    tinygltf::Accessor& accessor_1_triangles_vertices = m.accessors[ACCESSOR_1_TRIANGLES_VERTICES];

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

    if(with_wireframe) {
        m.accessors.resize(3);
        tinygltf::Accessor& accessor_2_wireframe_edges_vertices = m.accessors[ACCESSOR_2_WIREFRAME_EDGES_VERTICES];
        // layout description of buffer view 2
        accessor_2_wireframe_edges_vertices.bufferView = BUFFERVIEW_2_WIREFRAME_EDGES_VERTICES;
        accessor_2_wireframe_edges_vertices.byteOffset = 0;
        accessor_2_wireframe_edges_vertices.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT; // 32 bits per index
        accessor_2_wireframe_edges_vertices.count = wireframe_edges_as_vector.size()*2; // how many indices
        accessor_2_wireframe_edges_vertices.type = TINYGLTF_TYPE_SCALAR;
        // range of indices value : [ 0 : #vertices-1 ]
        accessor_2_wireframe_edges_vertices.maxValues.push_back(M.vertices.nb()-1); // max value of all indices in the buffer view
        accessor_2_wireframe_edges_vertices.minValues.push_back(0);
    }

    ///////////////////////////////////////////////////////////
    // Create a mesh with 1 primitive (+1 with the wireframe)
    ///////////////////////////////////////////////////////////

    m.meshes.resize(1);
    const size_t MESH_0 = 0;
    tinygltf::Mesh& mesh_0 = m.meshes[MESH_0];

    mesh_0.primitives.resize(1);
    const size_t PRIMITIVE_0_TRIANGLES = 0;
    const size_t PRIMITIVE_1_WIREFRAME_EDGES = 1;
    tinygltf::Primitive& primitive_0_triangles = mesh_0.primitives[PRIMITIVE_0_TRIANGLES];

    primitive_0_triangles.indices = ACCESSOR_1_TRIANGLES_VERTICES;
    primitive_0_triangles.attributes["POSITION"] = ACCESSOR_0_VERTICES_COORDINATES;
    primitive_0_triangles.material = MATERIAL_0_GREY;
    primitive_0_triangles.mode = TINYGLTF_MODE_TRIANGLES;

    if(with_wireframe) {
        mesh_0.primitives.resize(2);
        tinygltf::Primitive& primitive_1_wireframe_edges = mesh_0.primitives[PRIMITIVE_1_WIREFRAME_EDGES];
        primitive_1_wireframe_edges.indices = ACCESSOR_2_WIREFRAME_EDGES_VERTICES;
        primitive_1_wireframe_edges.attributes["POSITION"] = ACCESSOR_0_VERTICES_COORDINATES;
        primitive_1_wireframe_edges.material = MATERIAL_1_BLACK;
        primitive_1_wireframe_edges.mode = TINYGLTF_MODE_LINE;
    }

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

void write_glTF__labeled_triangle_mesh(std::string filename, GEO::Mesh& M, const char* attribute_name) {

    // for now, only write a simple texture
    // https://github.com/KhronosGroup/glTF-Tutorials/blob/main/gltfTutorial/gltfTutorial_013_SimpleTexture.md

    ////////////////////////
    // Create a glTF model
    ////////////////////////

    tinygltf::Model m;

    ////////////////////
    // Create an image
    ////////////////////

    m.images.resize(1);
    const size_t IMAGE_0 = 0;
    tinygltf::Image& image_0 = m.images[IMAGE_0];

    image_0.uri = "testTexture.png";
    image_0.mimeType = "image/png"; // should not be required if uri is specified

    /////////////////////
    // Create a sampler
    /////////////////////

    m.samplers.resize(1);
    const size_t SAMPLER_0 = 0;
    tinygltf::Sampler& sampler_0 = m.samplers[SAMPLER_0];

    sampler_0.magFilter = 9729; // ?
    sampler_0.minFilter = 9987; // ?
    sampler_0.wrapS = 33648; // ?
    sampler_0.wrapT = 33648; // ?

    /////////////////////
    // Create a texture
    /////////////////////

    m.textures.resize(1);
    const size_t TEXTURE_0 = 0;
    tinygltf::Texture& texture_0 = m.textures[TEXTURE_0];

    texture_0.sampler = SAMPLER_0;
    texture_0.source = 0; // IMAGE_0 ?

    //////////////////////
    // Create 1 material 
    //////////////////////

    m.materials.resize(1);
    const size_t MATERIAL_0 = 0;
    tinygltf::Material& material_0 = m.materials[MATERIAL_0];
    
    material_0.pbrMetallicRoughness.baseColorTexture.index = TEXTURE_0;
    material_0.pbrMetallicRoughness.metallicFactor = 0.0;
    material_0.pbrMetallicRoughness.roughnessFactor = 1.0;

    ////////////////////
    // Create a buffer
    ////////////////////

    m.buffers.resize(1);
    const size_t BUFFER_0 = 0;
    tinygltf::Buffer& buffer_0 = m.buffers[BUFFER_0];
    //
    //       base64 : AAABAAIAAQADAAIAAAAAAAAAAAAAAAAAAACAPwAAAAAAAAAAAAAAAAAAgD8AAAAAAACAPwAAgD8AAAAAAAAAAAAAgD8AAAAAAACAPwAAgD8AAAAAAAAAAAAAAAAAAAAAAACAPwAAAAAAAAAA (from gltfTutorial_013)
    // base16 (hex) : 0000010002000100030002000000000000000000000000000000803F0000000000000000000000000000803F000000000000803F0000803F00000000000000000000803F000000000000803F0000803F000000000000000000000000000000000000803F0000000000000000
    //
    // BUFFERVIEW_0 -> first 12 bytes
    // 000001000200010003000200
    // BUFFERVIEW_1 -> next 96 bytes
    // with a byteStride of 12 https://github.com/KhronosGroup/glTF-Tutorials/blob/main/gltfTutorial/gltfTutorial_005_BuffersBufferViewsAccessors.md#data-interleaving
    // 0000000000000000000000000000803F0000000000000000000000000000803F000000000000803F0000803F00000000000000000000803F000000000000803F0000803F000000000000000000000000000000000000803F0000000000000000
    // 
    // ACCESSOR_0 -> triangles = vertex indices
    //   triangle 0 (f0) : 00 00 01 00 02 00 -> vertices 0x0000, 0x0001 and 0x0002
    //   triangle 1 (f1) : 01 00 03 00 02 00 -> vertices 0x0001, 0x0003 and 0x0002
    // ACCESSOR_1 -> vertices = 3D coordinates = 3 float per vertex = 12 bytes per vertex
    //     00 00 00 00 00 00 00 00 00 00 00 00 -> v0 position : 0x00000000, 0x00000000, 0x00000000 -> {0.0f, 0.0f, 0.0f}
    //     00 00 80 3f 00 00 00 00 00 00 00 00 -> v1 position : 0x3f800000, 0x00000000, 0x00000000 -> {1.0f, 0.0f, 0.0f}
    //     00 00 00 00 00 00 80 3f 00 00 00 00 -> v2 position : 0x00000000, 0x3f800000, 0x00000000 -> {0.0f, 1.0f, 0.0f}
    //     00 00 80 3f 00 00 80 3f 00 00 00 00 -> v3 position : 0x3f800000, 0x3f800000, 0x00000000 -> {1.0f, 1.0f, 0.0f}
    // ACCESSOR_2 -> vertices texture = 2D coordinates. byteStride of 12 in BUFFERVIEW_1 -> a texture coordinate every 12 bytes, even if 2 floats have only 8 bytes
    //     00 00 00 00 00 00 80 3f 00 00 00 00 -> v0 texture coord. : 0x00000000, 0x3f800000, 0x00000000 -> {0.0f, 1.0f} + padding
    //     00 00 80 3f 00 00 80 3f 00 00 00 00 -> v1 texture coord. : 0x3f800000, 0x3f800000, 0x00000000 -> {1.0f, 1.0f} + padding
    //     00 00 00 00 00 00 00 00 00 00 00 00 -> v2 texture coord. : 0x00000000, 0x00000000, 0x00000000 -> {0.0f, 0.0f} + padding
    //     00 00 80 3f 00 00 00 00 00 00 00 00 -> v3 texture coord. : 0x3f800000, 0x00000000, 0x00000000 -> {1.0f, 0.0f} + padding
    //
    //       v2          v3                            0,0          1,0 
    //       +------------+                             +------------+
    //       | \\         |                             |            |
    //       |   \\   f1  |                             |            |
    //       |     \\     |      texture coordinates :  |            |
    //       |  f0   \\   |                             |            |
    // Y     |         \\ |                             |            |
    //       +------------+                             +------------+
    // ^     v0          v1                            0,1          1,1
    // |
    // o-- >  X

    buffer_0.uri = "data:application/gltf-buffer;base64,AAABAAIAAQADAAIAAAAAAAAAAAAAAAAAAACAPwAAAAAAAAAAAAAAAAAAgD8AAAAAAACAPwAAgD8AAAAAAAAAAAAAgD8AAAAAAACAPwAAgD8AAAAAAAAAAAAAAAAAAAAAAACAPwAAAAAAAAAA";
    // buffer_0.byteLength = 108;

    //////////////////////////
    // Create 2 buffer views
    //////////////////////////

    m.bufferViews.resize(2);
    const size_t BUFFERVIEW_0 = 0;
    const size_t BUFFERVIEW_1 = 1;
    tinygltf::BufferView& bufferview_0 = m.bufferViews[BUFFERVIEW_0];
    tinygltf::BufferView& bufferview_1 = m.bufferViews[BUFFERVIEW_1];
    
    bufferview_0.buffer = BUFFER_0;
    bufferview_0.byteOffset = 0;
    bufferview_0.byteLength = 12;
    bufferview_0.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    bufferview_1.buffer = BUFFER_0;
    bufferview_1.byteOffset = 12;
    bufferview_1.byteLength = 96;
    bufferview_1.byteStride = 12;
    bufferview_1.target = TINYGLTF_TARGET_ARRAY_BUFFER;

    ///////////////////////
    // Create 3 accessors 
    ///////////////////////

    m.accessors.resize(3);
    const size_t ACCESSOR_0 = 0;
    const size_t ACCESSOR_1 = 1;
    const size_t ACCESSOR_2 = 2;
    tinygltf::Accessor& accessor_0 = m.accessors[ACCESSOR_0];
    tinygltf::Accessor& accessor_1 = m.accessors[ACCESSOR_1];
    tinygltf::Accessor& accessor_2 = m.accessors[ACCESSOR_2];

    accessor_0.bufferView = BUFFERVIEW_0;
    accessor_0.byteOffset = 0;
    accessor_0.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT;
    accessor_0.count = 6;
    accessor_0.type = TINYGLTF_TYPE_SCALAR;
    accessor_0.maxValues = {3};
    accessor_0.minValues = {0};

    accessor_1.bufferView = BUFFERVIEW_1;
    accessor_1.byteOffset = 0;
    accessor_1.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_1.count = 4;
    accessor_1.type = TINYGLTF_TYPE_VEC3;
    accessor_1.maxValues = {1.0, 1.0, 0.0};
    accessor_1.minValues = {0.0, 0.0, 0.0};

    accessor_2.bufferView = BUFFERVIEW_1;
    accessor_2.byteOffset = 48;
    accessor_2.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_2.count = 4;
    accessor_2.type = TINYGLTF_TYPE_VEC2;
    accessor_2.maxValues = {1.0, 1.0};
    accessor_2.minValues = {0.0, 0.0};

    ///////////////////////////////////
    // Create a mesh with 1 primitive
    ///////////////////////////////////

    m.meshes.resize(1);
    const size_t MESH_0 = 0;
    tinygltf::Mesh& mesh_0 = m.meshes[MESH_0];

    mesh_0.primitives.resize(1);
    const size_t PRIMITIVE_0 = 0;
    tinygltf::Primitive& primitive_0 = mesh_0.primitives[PRIMITIVE_0];

    primitive_0.indices = ACCESSOR_0;
    primitive_0.attributes["POSITION"] = ACCESSOR_1;
    primitive_0.attributes["TEXCOORD_0"] = ACCESSOR_2;
    primitive_0.material = MATERIAL_0;
    primitive_0.mode = TINYGLTF_MODE_TRIANGLES; // default value in spec, but tinygltf writes -1 if not specified

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
    fmt::println(Logger::out("glTF"),"Writing SimpleTexture.gltf..."); Logger::out("glTF").flush();
    gltf.WriteGltfSceneToFile(&m, "SimpleTexture.gltf",
                            true, // embedImages
                            true, // embedBuffers
                            true, // pretty print
                            false); // write binary
}