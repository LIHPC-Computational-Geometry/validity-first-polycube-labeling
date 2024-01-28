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

#include <fmt/core.h>
#include <fmt/ostream.h>

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

    GEO::mesh_save(box,"box.geogram");

    // Create a glTF model with a single mesh and save it as a gltf file

    tinygltf::Model m;
    tinygltf::Scene scene;
    tinygltf::Mesh mesh;
    tinygltf::Primitive primitive;
    tinygltf::Node node;
    tinygltf::Buffer buffer; // the raw data buffer storing points & triangles
    tinygltf::BufferView bufferView1;
    tinygltf::BufferView bufferView2;
    tinygltf::Accessor accessor1;
    tinygltf::Accessor accessor2;
    tinygltf::Asset asset;

    // resize `buffer` to 4 bytes * 3 indices * #triangles
    //                  + 4 bytes * 3 floating points * #vertices
    size_t buffer_indices_start = 0; // byte index
    size_t buffer_indices_length = 4*3*box.facets.nb(); // number of bytes
    size_t buffer_coordinates_start = buffer_indices_length; // byte index. just after the indices chunk, no padding
    size_t buffer_coordinates_length = 4*3*box.vertices.nb(); // number of bytes
    buffer.data.resize(buffer_indices_length+buffer_coordinates_length);

    // write triangles (= vertex indices)
    FOR(f,box.facets.nb()) { // for each facet (triangle) index
        memcpy(buffer.data.data() + f*12 + 0, facet_vertex_index_ptr(box,f,0),4);
        memcpy(buffer.data.data() + f*12 + 4, facet_vertex_index_ptr(box,f,1),4);
        memcpy(buffer.data.data() + f*12 + 8, facet_vertex_index_ptr(box,f,2),4);
    }
    // write vertices (= 3D coordinates)
    FOR(v,box.vertices.nb()) { // for each vertex index
        memcpy(buffer.data.data()+buffer_coordinates_start+v*12+0,static_cast<void*>(box.vertices.single_precision_point_ptr(v)),12); // write x,y,z from float* getter
    }

    // 1st chunk of the buffer : triangles = vertex indices (type ELEMENT_ARRAY_BUFFER)
    bufferView1.buffer = 0;
    bufferView1.byteOffset=buffer_indices_start;
    bufferView1.byteLength=buffer_indices_length;
    bufferView1.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    // 2nd chunk of the bufffer : points = 3D coordinates (type ARRAY_BUFFER)
    bufferView2.buffer = 0;
    bufferView2.byteOffset=buffer_coordinates_start;
    bufferView2.byteLength=buffer_coordinates_length;
    bufferView2.target = TINYGLTF_TARGET_ARRAY_BUFFER;

    // layout description of bufferView1
    accessor1.bufferView = 0;
    accessor1.byteOffset = 0;
    accessor1.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT; // 32 bits per index
    accessor1.count = box.facets.nb()*3; // how many indices
    accessor1.type = TINYGLTF_TYPE_SCALAR;
    // range of indices : [ 0 : box.facets.nb()*3-1 ]
    accessor1.maxValues.push_back(box.facets.nb()*3-1);
    accessor1.minValues.push_back(0);

    // layout description of bufferView2
    accessor2.bufferView = 1;
    accessor2.byteOffset = 0;
    accessor2.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor2.count = box.vertices.nb(); // how many points
    accessor2.type = TINYGLTF_TYPE_VEC3;
    // range of coordinate values
    accessor2.maxValues = { 0.5,  0.5,  0.5};
    accessor2.minValues = {-0.5, -0.5, -0.5};

    // Build the mesh primitive and add it to the mesh
    primitive.indices = 0;                 // The index of the accessor for the vertex indices (`accessor1`)
    primitive.attributes["POSITION"] = 1;  // The index of the accessor for positions (`accessor2`)
    primitive.material = 0;
    primitive.mode = TINYGLTF_MODE_TRIANGLES;
    mesh.primitives.push_back(primitive);

    // Other tie ups
    node.mesh = 0;
    scene.nodes.push_back(0); // Default scene

    // Define the asset. The version is required
    asset.version = "2.0";
    asset.generator = "tinygltf";

    // Now all that remains is to tie back all the loose objects into the
    // our single model.
    m.scenes.push_back(scene);
    m.meshes.push_back(mesh);
    m.nodes.push_back(node);
    m.buffers.push_back(buffer);
    m.bufferViews.push_back(bufferView1);
    m.bufferViews.push_back(bufferView2);
    m.accessors.push_back(accessor1);
    m.accessors.push_back(accessor2);
    m.asset = asset;

    // Create a simple material
    tinygltf::Material mat;
    mat.pbrMetallicRoughness.baseColorFactor = {1.0f, 0.9f, 0.9f, 1.0f};  
    mat.doubleSided = true;
    m.materials.push_back(mat);

    // Save it to a file
    tinygltf::TinyGLTF gltf;
    gltf.WriteGltfSceneToFile(&m, "box.glb",
                            true, // embedImages
                            true, // embedBuffers
                            true, // pretty print
                            true); // write binary
  
  exit(0);
}