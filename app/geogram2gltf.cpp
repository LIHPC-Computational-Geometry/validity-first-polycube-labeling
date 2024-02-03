// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "tiny_gltf.h"

#include <geogram/basic/command_line_args.h>    // for import_arg_group()
#include <geogram/basic/geometry.h>             // for vec3
#include <geogram/basic/logger.h>               // for Logger
#include <geogram/mesh/mesh.h>                  // for Mesh
#include <geogram/mesh/mesh_io.h>               // for mesh_save()
#include <geogram/mesh/mesh_halfedges.h>        // for Halfedge

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

using GEO::index_t; // to use the FOR() macro of Geogram

inline index_t* facet_vertex_index_ptr(GEO::Mesh& M, index_t f, index_t lv) {
    return M.facet_corners.vertex_index_ptr(M.facets.corner(f,lv));
}

int main(int argc, char **argv)
{
    // based on
    // - https://github.com/syoyo/tinygltf/issues/63
    // - https://gitlab.com/dodgyville/pygltflib#a-simple-mesh
    // - https://gist.github.com/dov/2f730851aaaf84b85539bba5ed5fb268
    // - https://gitlab.com/dodgyville/pygltflib#create-a-mesh-convert-to-bytes-convert-back-to-mesh

    GEO::initialize();
    GEO::CmdLine::import_arg_group("sys"); // declares sys:compression_level, needed by mesh_save() for .geogram files

    #ifdef GARGANTUA
        fmt::println(Logger::err("I/O"),"Cannot export to glTF in GARGANTUA mode, indices must have 32 bits"); Logger::err("I/O").flush();
        return 1;
    #endif

    geo_assert(sizeof(index_t)==4); // 4 bytes = 32 bits

    // create a simple Geogram mesh (a box)

    GEO::Mesh box;
    box.vertices.create_vertices(8);
    box.vertices.point(0) = GEO::vec3(-0.5, -0.5,  0.5);
    box.vertices.point(1) = GEO::vec3( 0.5, -0.5,  0.5);
    box.vertices.point(2) = GEO::vec3(-0.5,  0.5,  0.5);
    box.vertices.point(3) = GEO::vec3( 0.5,  0.5,  0.5);
    box.vertices.point(4) = GEO::vec3( 0.5, -0.5, -0.5);
    box.vertices.point(5) = GEO::vec3(-0.5, -0.5, -0.5);
    box.vertices.point(6) = GEO::vec3( 0.5,  0.5, -0.5);
    box.vertices.point(7) = GEO::vec3(-0.5,  0.5, -0.5);
    box.vertices.set_single_precision();
    box.facets.create_triangles(12);
    // triangle 0
    box.facets.set_vertex(0,0,0);
    box.facets.set_vertex(0,1,1);
    box.facets.set_vertex(0,2,2);
    // triangle 1
    box.facets.set_vertex(1,0,3);
    box.facets.set_vertex(1,1,2);
    box.facets.set_vertex(1,2,1);
    // triangle 2
    box.facets.set_vertex(2,0,1);
    box.facets.set_vertex(2,1,0);
    box.facets.set_vertex(2,2,4);
    // triangle 3
    box.facets.set_vertex(3,0,5);
    box.facets.set_vertex(3,1,4);
    box.facets.set_vertex(3,2,0);
    // triangle 4
    box.facets.set_vertex(4,0,3);
    box.facets.set_vertex(4,1,1);
    box.facets.set_vertex(4,2,6);
    // triangle 5
    box.facets.set_vertex(5,0,4);
    box.facets.set_vertex(5,1,6);
    box.facets.set_vertex(5,2,1);
    // triangle 6
    box.facets.set_vertex(6,0,2);
    box.facets.set_vertex(6,1,3);
    box.facets.set_vertex(6,2,7);
    // triangle 7
    box.facets.set_vertex(7,0,6);
    box.facets.set_vertex(7,1,7);
    box.facets.set_vertex(7,2,3);
    // triangle 8
    box.facets.set_vertex(8,0,0);
    box.facets.set_vertex(8,1,2);
    box.facets.set_vertex(8,2,5);
    // triangle 9
    box.facets.set_vertex(9,0,7);
    box.facets.set_vertex(9,1,5);
    box.facets.set_vertex(9,2,2);
    // triangle 10
    box.facets.set_vertex(10,0,5);
    box.facets.set_vertex(10,1,7);
    box.facets.set_vertex(10,2,4);
    // triangle 11
    box.facets.set_vertex(11,0,6);
    box.facets.set_vertex(11,1,4);
    box.facets.set_vertex(11,2,7);
    // compute facets adjacency
    box.facets.connect();
    // fill edges. for unicity, only add halfedges H where v0 < v1
    GEO::CustomMeshHalfedges mesh_he(box);
    GEO::MeshHalfedges::Halfedge H;
    index_t H_v0 = index_t(-1);
    index_t H_v1 = index_t(-1);
    FOR(f,box.facets.nb()) { // for each facet
        H.facet = f;
        FOR(lv,3) { // for each local vertex of f
            H.corner = box.facets.corner(f,lv);
            H_v0 = halfedge_vertex_index_from(box,H);
            H_v1 = halfedge_vertex_index_to(box,H);
            if(H_v0 < H_v1) {
                box.edges.create_edge(H_v0,H_v1);
            }
        }
    }

    GEO::mesh_save(box,"box.geogram");

    ////////////////////////
    // Create a glTF model
    ////////////////////////

    tinygltf::Model m;

    ///////////////////////
    // Create 2 materials
    ///////////////////////

    m.materials.resize(2);
    const size_t MATERIAL_0_RED = 0;
    const size_t MATERIAL_1_BLACK = 1;
    tinygltf::Material& material_0_red = m.materials[MATERIAL_0_RED];
    tinygltf::Material& material_1_black = m.materials[MATERIAL_1_BLACK];

    material_0_red.pbrMetallicRoughness.baseColorFactor = {1.0f, 0.9f, 0.9f, 1.0f};
    material_0_red.doubleSided = true;

    material_1_black.pbrMetallicRoughness.baseColorFactor = {0.0f, 0.0f, 0.0f, 1.0f};
    material_1_black.doubleSided = true;

    ////////////////////
    // Create a buffer
    ////////////////////

    m.buffers.resize(1);
    const size_t BUFFER_0 = 0;
    tinygltf::Buffer& buffer_0 = m.buffers[BUFFER_0];

    // resize buffer to 4 bytes * 3 indices * #triangles
    //                + 4 bytes * 3 floating points * #vertices
    //                + 4 bytes * 2 indices * #edges
    size_t buffer_indices_start = 0; // byte index
    size_t buffer_indices_length = 4*3*box.facets.nb(); // number of bytes
    size_t buffer_coordinates_start = buffer_indices_length; // byte index. just after the indices chunk, no padding
    size_t buffer_coordinates_length = 4*3*box.vertices.nb(); // number of bytes
    size_t buffer_edges_start = buffer_coordinates_start+buffer_coordinates_length; // byte index
    size_t buffer_edges_length = 4*2*box.edges.nb(); // number of bytes
    buffer_0.data.resize(buffer_indices_length+buffer_coordinates_length+buffer_edges_length);

    // write triangles (= vertex indices)
    FOR(f,box.facets.nb()) { // for each facet (triangle) index
        memcpy(buffer_0.data.data() + f*12 + 0, facet_vertex_index_ptr(box,f,0),4);
        memcpy(buffer_0.data.data() + f*12 + 4, facet_vertex_index_ptr(box,f,1),4);
        memcpy(buffer_0.data.data() + f*12 + 8, facet_vertex_index_ptr(box,f,2),4);
    }
    // write vertices (= 3D coordinates)
    FOR(v,box.vertices.nb()) { // for each vertex index
        memcpy(buffer_0.data.data()+buffer_coordinates_start+v*12+0,static_cast<void*>(box.vertices.single_precision_point_ptr(v)),12); // write x,y,z from float* getter
    }
    // write edges (= vertex indices)
    FOR(e,box.edges.nb()) {
        memcpy(buffer_0.data.data()+buffer_edges_start + 8*e + 0, box.edges.vertex_index_ptr(2*e+0), 4);
        memcpy(buffer_0.data.data()+buffer_edges_start + 8*e + 4, box.edges.vertex_index_ptr(2*e+1), 4);
    }    

    //////////////////////////
    // Create 3 buffer views
    //////////////////////////

    m.bufferViews.resize(3);
    const size_t BUFFERVIEW_0_TRIANGLES_VERTICES = 0;
    const size_t BUFFERVIEW_1_VERTICES_COORDINATES = 1;
    const size_t BUFFERVIEW_2_EDGES_VERTICES = 2;
    tinygltf::BufferView& bufferview_0_triangles_vertices = m.bufferViews[BUFFERVIEW_0_TRIANGLES_VERTICES];
    tinygltf::BufferView& bufferview_1_vertices_coordinates = m.bufferViews[BUFFERVIEW_1_VERTICES_COORDINATES];
    tinygltf::BufferView& bufferview_2_edges_vertices = m.bufferViews[BUFFERVIEW_2_EDGES_VERTICES];

    // 1st chunk of the buffer 0 : triangles = vertex indices (type ELEMENT_ARRAY_BUFFER)
    bufferview_0_triangles_vertices.buffer = BUFFER_0;
    bufferview_0_triangles_vertices.byteOffset = buffer_indices_start;
    bufferview_0_triangles_vertices.byteLength = buffer_indices_length;
    bufferview_0_triangles_vertices.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    // 2nd chunk of the bufffer 0 : points = 3D coordinates (type ARRAY_BUFFER)
    bufferview_1_vertices_coordinates.buffer = BUFFER_0;
    bufferview_1_vertices_coordinates.byteOffset = buffer_coordinates_start;
    bufferview_1_vertices_coordinates.byteLength = buffer_coordinates_length;
    bufferview_1_vertices_coordinates.target = TINYGLTF_TARGET_ARRAY_BUFFER;

    // 3rd chunk of the buffer 0 : indices of edge vertices (type ELEMENT_ARRAY_BUFFER)
    bufferview_2_edges_vertices.buffer = BUFFER_0;
    bufferview_2_edges_vertices.byteOffset = buffer_edges_start;
    bufferview_2_edges_vertices.byteLength = buffer_edges_length;
    bufferview_2_edges_vertices.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    ///////////////////////
    // Create 3 accessors
    ///////////////////////

    m.accessors.resize(3);
    const size_t ACCESSOR_0_TRIANGLES_VERTICES = 0;
    const size_t ACCESSOR_1_VERTICES_COORDINATES = 1;
    const size_t ACCESSOR_2_EDGES_VERTICES = 2;
    tinygltf::Accessor& accessor_0_triangles_vertices = m.accessors[ACCESSOR_0_TRIANGLES_VERTICES];
    tinygltf::Accessor& accessor_1_vertices_coordinates = m.accessors[ACCESSOR_1_VERTICES_COORDINATES];
    tinygltf::Accessor& accessor_2_edges_vertices = m.accessors[ACCESSOR_2_EDGES_VERTICES];

    // layout description of buffer view 0
    accessor_0_triangles_vertices.bufferView = BUFFERVIEW_0_TRIANGLES_VERTICES;
    accessor_0_triangles_vertices.byteOffset = 0;
    accessor_0_triangles_vertices.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT; // 32 bits per index
    accessor_0_triangles_vertices.count = box.facets.nb()*3; // how many indices
    accessor_0_triangles_vertices.type = TINYGLTF_TYPE_SCALAR;
    // range of indices : [ 0 : box.facets.nb()*3-1 ]
    accessor_0_triangles_vertices.maxValues.push_back(box.facets.nb()*3-1);
    accessor_0_triangles_vertices.minValues.push_back(0);

    // layout description of buffer view 1
    accessor_1_vertices_coordinates.bufferView = BUFFERVIEW_1_VERTICES_COORDINATES;
    accessor_1_vertices_coordinates.byteOffset = 0;
    accessor_1_vertices_coordinates.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_1_vertices_coordinates.count = box.vertices.nb(); // how many points
    accessor_1_vertices_coordinates.type = TINYGLTF_TYPE_VEC3;
    // range of coordinate values
    accessor_1_vertices_coordinates.maxValues = { 0.5,  0.5,  0.5};
    accessor_1_vertices_coordinates.minValues = {-0.5, -0.5, -0.5};

    // layout description of buffer view 2
    accessor_2_edges_vertices.bufferView = BUFFERVIEW_2_EDGES_VERTICES;
    accessor_2_edges_vertices.byteOffset = 0;
    accessor_2_edges_vertices.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT; // 32 bits per index
    accessor_2_edges_vertices.count = box.edges.nb(); // how many indices
    accessor_2_edges_vertices.type = TINYGLTF_TYPE_SCALAR;
    // range of indices : [ 0 : 1 ]
    accessor_2_edges_vertices.maxValues.push_back(box.edges.nb()-1);
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

    primitive_0_triangles.indices = ACCESSOR_0_TRIANGLES_VERTICES;
    primitive_0_triangles.attributes["POSITION"] = ACCESSOR_1_VERTICES_COORDINATES;
    primitive_0_triangles.material = MATERIAL_0_RED;
    primitive_0_triangles.mode = TINYGLTF_MODE_TRIANGLES;

    primitive_1_edges.indices = ACCESSOR_2_EDGES_VERTICES;
    primitive_1_edges.attributes["POSITION"] = ACCESSOR_1_VERTICES_COORDINATES;
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
    gltf.WriteGltfSceneToFile(&m, "box.glb",
                            true, // embedImages
                            true, // embedBuffers
                            true, // pretty print
                            true); // write binary
  
  exit(0);
}