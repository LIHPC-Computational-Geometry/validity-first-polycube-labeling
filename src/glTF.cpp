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
#include "containers.h"             // for index_of_last()
#include "geometry.h"               // for AdjacentFacetOfVertex
#include "hex_mesh.h"               // for compute_scaled_jacobian()
#include "dump_mesh.h"              // for dump_edges()

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
    fmt::println(Logger::out("glTF"),"Writing {}...",filename); Logger::out("glTF").flush();
    gltf.WriteGltfSceneToFile(&m, filename,
                            true, // embedImages
                            true, // embedBuffers
                            true, // pretty print
                            true); // write binary

}

void write_glTF__labeled_triangle_mesh(std::string filename, GEO::Mesh& M, const char* attribute_name, const std::vector<std::vector<AdjacentFacetOfVertex>>& per_vertex_adj_facets) {
    ASSERT_GARGANTUA_OFF;
    geo_assert(per_vertex_adj_facets.size() == M.vertices.nb());

    // retreive per facet label
    GEO::Attribute<index_t> label(M.facets.attributes(),attribute_name);
    // create a vector that will enumerate vertices in the glTF mesh
    std::vector<glTF_vertex<index_t>> glTF_vertices; // for each vertex in the glTF mesh, store the corresponding vertex in M and the label (future texture coordinate) of this glTF vertex
    // create per facet corner attribute to store the glTF mesh vertex index (can be different of the index in M)
    GEO::Attribute<index_t> facet_corner_to_glTF_vertex(M.facet_corners.attributes(),"glTF_vertex");
    
    std::vector<AdjacentFacetOfVertex>::const_iterator it;
    index_t label_around_vertex = index_t(-1);
    bool same_label_all_around_vertex = true;
    index_t facet_corner = index_t(-1);
    FOR(v,M.vertices.nb()) { // for each vertex of M
        same_label_all_around_vertex = true;
        // get label of first facet around vertex v
        label_around_vertex = label[per_vertex_adj_facets[v][0].facet_index];
        // parse remaining facets around v
        it = per_vertex_adj_facets[v].begin();
        it++; // start with facet at index 1 (we already have the label of facet 0)
        do {
            if(label[it->facet_index] != label_around_vertex) {
                // facets around vertex v don't have the same label
                same_label_all_around_vertex = false;
                break;
            }
            it++;
        } while (it != per_vertex_adj_facets[v].end());

        if(same_label_all_around_vertex) {
            // only one vertex in the glTF mesh in necessary to represent v^th vertex of M
            // create new glTF vertex and link it to M vertex
            glTF_vertices.push_back(glTF_vertex{v,label_around_vertex});
            // link M facet corners to glTF vertex
            for(const auto& adj_facets : per_vertex_adj_facets[v]) {
                facet_corner = M.facets.corner(adj_facets.facet_index,adj_facets.local_vertex);
                facet_corner_to_glTF_vertex[facet_corner] = (index_t) index_of_last(glTF_vertices);
            }
        }
        else {
            // split this vertex in as many vertices as there are adjacent facets
            for(const auto& adj_facets : per_vertex_adj_facets[v]) {
                // create new glTF vertex and link it to M vertex
                glTF_vertices.push_back(glTF_vertex{v,label[adj_facets.facet_index]});
                // link M facet corners to glTF vertex
                facet_corner = M.facets.corner(adj_facets.facet_index,adj_facets.local_vertex);
                facet_corner_to_glTF_vertex[facet_corner] = (index_t) index_of_last(glTF_vertices);
            }
        }
    }

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

    ////////////////////
    // Create an image
    ////////////////////

    m.images.resize(1);
    const size_t IMAGE_0 = 0;
    tinygltf::Image& image_0 = m.images[IMAGE_0];

    image_0.uri = labeling_texture_image_PNG;

    /////////////////////
    // Create a sampler
    /////////////////////

    m.samplers.resize(1);
    const size_t SAMPLER_0 = 0;
    tinygltf::Sampler& sampler_0 = m.samplers[SAMPLER_0];

    sampler_0.magFilter = TINYGLTF_TEXTURE_FILTER_LINEAR;       // magnification filter
    sampler_0.minFilter = TINYGLTF_TEXTURE_FILTER_LINEAR;       // minification filter
    sampler_0.wrapS = TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT;    // s wrapping mode
    sampler_0.wrapT = TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT;    // t wrapping mode

    /////////////////////
    // Create a texture
    /////////////////////

    m.textures.resize(1);
    const size_t TEXTURE_0 = 0;
    tinygltf::Texture& texture_0 = m.textures[TEXTURE_0];

    texture_0.sampler = SAMPLER_0;
    texture_0.source = IMAGE_0;

    //////////////////////
    // Create a material 
    //////////////////////

    m.materials.resize(1);
    const size_t MATERIAL_0 = 0;
    tinygltf::Material& material_0 = m.materials[MATERIAL_0];
    
    material_0.pbrMetallicRoughness.baseColorTexture.index = TEXTURE_0;
    material_0.pbrMetallicRoughness.metallicFactor = 0.0;
    material_0.pbrMetallicRoughness.roughnessFactor = 1.0;
    material_0.doubleSided = true;

    ////////////////////
    // Create a buffer
    ////////////////////

    m.buffers.resize(1);
    const size_t BUFFER_0 = 0;
    tinygltf::Buffer& buffer_0 = m.buffers[BUFFER_0];

    // resize buffer to 4 bytes * 3 indices * #triangles (2*3*#triangles in the tutorial, but GEO::index_t is 4 bytes)
    //                + 4 bytes * 3 floating points * #vertices
    //                + 4 bytes * 3 floating points * #vertices (2 floating points + 1 of padding, per vertex -> byteStride of 12 bytes)
    size_t buffer_triangles_start               = 0;                                                                    // byte index
    size_t buffer_triangles_length              = 4*3*M.facets.nb();                                                    // number of bytes
    size_t buffer_vertices_coordinates_start    = buffer_triangles_length;                                              // byte index
    size_t buffer_vertices_coordinates_length   = 4*3*glTF_vertices.size();                                        // number of bytes
    size_t buffer_texture_coordinates_start     = buffer_vertices_coordinates_start+buffer_vertices_coordinates_length; // byte index
    size_t buffer_texture_coordinates_length    = 4*3*glTF_vertices.size();                                        // number of bytes
    buffer_0.data.resize(buffer_triangles_length+buffer_vertices_coordinates_length+buffer_texture_coordinates_length);

    // write triangles (= vertex indices)
    FOR(f,M.facets.nb()) { // for each facet (triangle) index
        FOR(lv,3) { // for each local vertex of this facet
            facet_corner = M.facets.corner(f,lv);
            memcpy(buffer_0.data.data() + f*12 + 4*lv, &facet_corner_to_glTF_vertex[facet_corner], 4);
        }
    }
    // write vertices position coordinates (= 3D coordinates)
    // & write vertices texture coordinates (= 2D coordinates aligned on 3D coordiates)
    uint16_t zeros_on_4_bytes = 0x00000000;
    FOR(glTF_v,glTF_vertices.size()) { // for each vertex index
        // write position
        memcpy(buffer_0.data.data() + buffer_vertices_coordinates_start + glTF_v*12, M.vertices.single_precision_point_ptr(glTF_vertices[glTF_v].Geogram_vertex), 12); // write x,y,z from float* getter
        // write texture
        memcpy(buffer_0.data.data() + buffer_texture_coordinates_start + glTF_v*12, LABELING_TO_TEXTURE_COORDINATES[glTF_vertices[glTF_v].value].data(), 8);
        memcpy(buffer_0.data.data() + buffer_texture_coordinates_start + glTF_v*12 + 8, &zeros_on_4_bytes, 4); // padding
    }

    //////////////////////////
    // Create 2 buffer views
    //////////////////////////

    m.bufferViews.resize(2);
    const size_t BUFFERVIEW_0 = 0;
    const size_t BUFFERVIEW_1 = 1;
    tinygltf::BufferView& bufferview_0 = m.bufferViews[BUFFERVIEW_0];
    tinygltf::BufferView& bufferview_1 = m.bufferViews[BUFFERVIEW_1];
    
    bufferview_0.buffer = BUFFER_0;
    bufferview_0.byteOffset = buffer_triangles_start;
    bufferview_0.byteLength = buffer_triangles_length;
    bufferview_0.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    bufferview_1.buffer = BUFFER_0;
    bufferview_1.byteOffset = buffer_vertices_coordinates_start;
    bufferview_1.byteLength = buffer_vertices_coordinates_length+buffer_texture_coordinates_length;
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
    accessor_0.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
    accessor_0.count = M.facets.nb()*3;
    accessor_0.type = TINYGLTF_TYPE_SCALAR;
    accessor_0.maxValues = {(double) (glTF_vertices.size()-1)};
    accessor_0.minValues = {0};

    accessor_1.bufferView = BUFFERVIEW_1;
    accessor_1.byteOffset = 0;
    accessor_1.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_1.count = glTF_vertices.size();
    accessor_1.type = TINYGLTF_TYPE_VEC3;
    accessor_1.maxValues = bounding_box_max;
    accessor_1.minValues = bounding_box_min;

    accessor_2.bufferView = BUFFERVIEW_1;
    accessor_2.byteOffset = buffer_texture_coordinates_start-bufferview_1.byteOffset; // byte offset relative to the bufferview
    accessor_2.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_2.count = glTF_vertices.size();
    accessor_2.type = TINYGLTF_TYPE_VEC2;
    accessor_2.maxValues = {LABELING_TO_TEXTURE_COORDINATES[5][0], 0.5};
    accessor_2.minValues = {LABELING_TO_TEXTURE_COORDINATES[0][0], 0.5};

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
    fmt::println(Logger::out("glTF"),"Writing {}...",filename); Logger::out("glTF").flush();
    gltf.WriteGltfSceneToFile(&m, filename,
                            true, // embedImages
                            true, // embedBuffers
                            true, // pretty print
                            true); // write binary
}

void write_glTF__hex_mesh_surface(std::string filename, const GEO::Mesh& hex_mesh) {
    // 1. Compute the per-cell Scaled Jacobian
    // 2. Extract the surface mesh
    // 3. Store the quads wireframe edges
    // 4. Triangulate the quads
    // 5. Compute per-vertex adjacent triangles
    // 6. For each vertex, create 1 glTF vertex if attribute value doesn't change much in neighborhood,
    //    else create a glTF vertex for each adjacent triangle
    // 7. Assemble the glTF asset, embed the colormap
    // 8. Write to file

    GEO::Mesh mesh;
    mesh.copy(hex_mesh,false); // don't copy attributes

    // 1.

    BasicStats SJ_stats;
    compute_scaled_jacobian(mesh,SJ_stats); // creates a "SJ" cell attribute, of type double
    GEO::Attribute<double> SJ(mesh.cells.attributes(), "SJ"); // retrieve the attribute

    // 2.

    GEO::Attribute<index_t> per_quad_cell_index(mesh.cells.attributes(),"cell_index");
    mesh.cells.compute_borders(per_quad_cell_index);
    // transfer attribute from cells to quads
    std::vector<double> SJ_on_surface_quads(mesh.facets.nb());
    FOR(f,mesh.facets.nb()) { // for each facet (quads)
        SJ_on_surface_quads[f] = SJ[per_quad_cell_index[f]];
    }
    mesh.cells.clear();
    mesh.vertices.remove_isolated();

    // 3.

    std::set<std::set<index_t>> wireframe_edges_as_set;
    GEO::CustomMeshHalfedges mesh_he(mesh);
    GEO::MeshHalfedges::Halfedge H;
    index_t H_v0 = index_t(-1);
    index_t H_v1 = index_t(-1);
    FOR(f,mesh.facets.nb()) { // for each facet
        H.facet = f;
        geo_assert(mesh.facets.nb_vertices(f)==4); // assert it's a quad
        FOR(lv,4) { // for each local vertex of f
            H.corner = mesh.facets.corner(f,lv);
            H_v0 = halfedge_vertex_index_from(mesh,H);
            H_v1 = halfedge_vertex_index_to(mesh,H);
            wireframe_edges_as_set.insert({H_v0,H_v1});
        }
    }
    // we need to have an index associated to each wireframe edge in the glTF asset
    // -> convert to vector
    std::vector<std::set<index_t>> wireframe_edges_as_vector(wireframe_edges_as_set.begin(),wireframe_edges_as_set.end());

    // 4.

    std::vector<index_t> triangle_to_old_facet; // index map before/after triangulation
    triangulate_facets(mesh,triangle_to_old_facet);
    // transfer attribute from quads to triangles
    GEO::Attribute<double> SJ_on_surface_triangles(mesh.facets.attributes(),"SJ");
    FOR(t,mesh.facets.nb()) { // for each triangle
        SJ_on_surface_triangles[t] = SJ_on_surface_quads[triangle_to_old_facet[t]];
    }

    // 5.

    std::vector<std::vector<AdjacentFacetOfVertex>> per_facet_adj_facet_corners;
    compute_adjacent_facets_of_vertices(mesh,per_facet_adj_facet_corners);

    // 6.

    // create a vector that will enumerate vertices in the glTF mesh
    std::vector<glTF_vertex<double>> glTF_vertices; // for each vertex in the glTF mesh, store the corresponding vertex in `mesh` and the SJ value (future texture coordinate) of this glTF vertex
    // create per facet corner attribute to store the glTF mesh vertex index (can be different of the index in `mesh`)
    GEO::Attribute<index_t> facet_corner_to_glTF_vertex(mesh.facet_corners.attributes(),"glTF_vertex");
    
    std::vector<AdjacentFacetOfVertex>::const_iterator it;
    double minSJ_around_vertex = index_t(-1);
    double maxSJ_around_vertex = index_t(-1);
    bool similar_SJ_around_vertex = true;
    index_t facet_corner = index_t(-1);
    FOR(v,mesh.vertices.nb()) { // for each vertex of M
        similar_SJ_around_vertex = true;
        // get label of first facet around vertex v
        minSJ_around_vertex = SJ_on_surface_triangles[per_facet_adj_facet_corners[v][0].facet_index];
        maxSJ_around_vertex = minSJ_around_vertex;
        // parse remaining facets around v
        it = per_facet_adj_facet_corners[v].begin();
        it++; // start with facet at index 1 (we already have the label of facet 0)
        do {
            minSJ_around_vertex = std::min(minSJ_around_vertex,SJ_on_surface_triangles[it->facet_index]);
            maxSJ_around_vertex = std::max(maxSJ_around_vertex,SJ_on_surface_triangles[it->facet_index]);
            if(maxSJ_around_vertex-minSJ_around_vertex > 0.05) {
                // facets around vertex v don't have a similar SJ value
                similar_SJ_around_vertex = false;
                break;
            }
            it++;
        } while (it != per_facet_adj_facet_corners[v].end());

        if(similar_SJ_around_vertex) {
            // only one vertex in the glTF mesh in necessary to represent v^th vertex of M
            // create new glTF vertex and link it to the vertex of `mesh`
            glTF_vertices.push_back(glTF_vertex{v,(maxSJ_around_vertex-minSJ_around_vertex)/2.0});
            // link `mesh` facet corners to glTF vertex
            for(const auto& adj_facets : per_facet_adj_facet_corners[v]) {
                facet_corner = mesh.facets.corner(adj_facets.facet_index,adj_facets.local_vertex);
                facet_corner_to_glTF_vertex[facet_corner] = (index_t) index_of_last(glTF_vertices);
            }
        }
        else {
            // split this vertex in as many vertices as there are adjacent facets
            for(const auto& adj_facets : per_facet_adj_facet_corners[v]) {
                // create new glTF vertex and link it to M vertex
                glTF_vertices.push_back(glTF_vertex{v,SJ_on_surface_triangles[adj_facets.facet_index]});
                // link M facet corners to glTF vertex
                facet_corner = mesh.facets.corner(adj_facets.facet_index,adj_facets.local_vertex);
                facet_corner_to_glTF_vertex[facet_corner] = (index_t) index_of_last(glTF_vertices);
            }
        }
    }
    // update `wireframe_edges_as_vector`
    index_t v0 = index_t(-1);
    index_t v1 = index_t(-1);
    FOR(i,wireframe_edges_as_vector.size()) {
        v0 = *wireframe_edges_as_vector[i].begin();
        v1 = *wireframe_edges_as_vector[i].rbegin();
        v0 = facet_corner_to_glTF_vertex[mesh.facets.corner(
            per_facet_adj_facet_corners[v0][0].facet_index,
            per_facet_adj_facet_corners[v0][0].local_vertex
        )];
        v1 = facet_corner_to_glTF_vertex[mesh.facets.corner(
            per_facet_adj_facet_corners[v1][0].facet_index,
            per_facet_adj_facet_corners[v1][0].local_vertex
        )];
        wireframe_edges_as_vector[i] = {v0,v1};
    }

    // 7.

    // Warning: vertices coordinates must be double precision for get_bbox()
    // Should be the case because of compute_scaled_jacobian()
    geo_assert(mesh.vertices.double_precision());
    std::vector<double> bounding_box_min(3,0.0), bounding_box_max(3,0.0);
    GEO::get_bbox(mesh, bounding_box_min.data(), bounding_box_max.data());

    // convert vertices coordinates to float
    mesh.vertices.set_single_precision();

    // create a glTF model
    tinygltf::Model m;

    // embed the Parula colormap
    m.images.resize(1);
    const size_t IMAGE_0 = 0;
    tinygltf::Image& image_0 = m.images[IMAGE_0];
    image_0.uri = parula_texture_image_PNG;

    // create a sampler
    m.samplers.resize(1);
    const size_t SAMPLER_0 = 0;
    tinygltf::Sampler& sampler_0 = m.samplers[SAMPLER_0];
    sampler_0.magFilter = TINYGLTF_TEXTURE_FILTER_LINEAR;   // magnification filter
    sampler_0.minFilter = TINYGLTF_TEXTURE_FILTER_LINEAR;   // minification filter
    sampler_0.wrapS = TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE;  // s wrapping mode
    sampler_0.wrapT = TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE;  // t wrapping mode

    // create a texture
    m.textures.resize(1);
    const size_t TEXTURE_0 = 0;
    tinygltf::Texture& texture_0 = m.textures[TEXTURE_0];
    texture_0.sampler = SAMPLER_0;
    texture_0.source = IMAGE_0;

    // create a material
    m.materials.resize(2);
    const size_t MATERIAL_0 = 0;
    const size_t MATERIAL_1 = 1; // for the wireframe edges
    tinygltf::Material& material_0 = m.materials[MATERIAL_0];
    tinygltf::Material& material_1 = m.materials[MATERIAL_1];
    material_0.pbrMetallicRoughness.baseColorTexture.index = TEXTURE_0;
    material_0.pbrMetallicRoughness.metallicFactor = 0.0;
    material_0.pbrMetallicRoughness.roughnessFactor = 1.0;
    material_0.doubleSided = true;
    material_1.pbrMetallicRoughness.baseColorFactor = {0.0f, 0.0f, 0.0f, 1.0f};
    material_1.doubleSided = true;

    // create a buffer
    m.buffers.resize(1);
    const size_t BUFFER_0 = 0;
    tinygltf::Buffer& buffer_0 = m.buffers[BUFFER_0];
    // resize buffer to 4 bytes * 3 indices * #triangles (GEO::index_t is 4 bytes)
    //                + 4 bytes * 3 floating points * #vertices
    //                + 4 bytes * 3 floating points * #vertices (2 floating points + 1 of padding, per vertex -> byteStride of 12 bytes)
    //                + 4 bytes * 2 indices * #edges
    size_t buffer_triangles_start               = 0;
    size_t buffer_triangles_length              = 4*3*mesh.facets.nb();
    size_t buffer_vertices_coordinates_start    = buffer_triangles_length;
    size_t buffer_vertices_coordinates_length   = 4*3*glTF_vertices.size();
    size_t buffer_texture_coordinates_start     = buffer_vertices_coordinates_start+buffer_vertices_coordinates_length;
    size_t buffer_texture_coordinates_length    = 4*3*glTF_vertices.size();
    size_t buffer_wireframe_edges_start         = buffer_texture_coordinates_start+buffer_texture_coordinates_length;
    size_t buffer_wireframe_edges_length        = 4*2*wireframe_edges_as_vector.size();
    buffer_0.data.resize(buffer_triangles_length+buffer_vertices_coordinates_length+buffer_texture_coordinates_length+buffer_wireframe_edges_length);

    // write triangles (= vertex indices)
    FOR(f,mesh.facets.nb()) { // for each facet (triangle) index
        FOR(lv,3) { // for each local vertex of this facet
            facet_corner = mesh.facets.corner(f,lv);
            memcpy(buffer_0.data.data() + f*12 + 4*lv, &facet_corner_to_glTF_vertex[facet_corner], 4);
        }
    }

    // write vertices position coordinates (= 3D coordinates)
    // & write vertices texture coordinates (= 2D coordinates aligned on 3D coordinates)
    uint16_t zeros_on_4_bytes = 0x00000000;
    FOR(glTF_v,glTF_vertices.size()) { // for each vertex index
        // write position
        memcpy(buffer_0.data.data() + buffer_vertices_coordinates_start + glTF_v*12, mesh.vertices.single_precision_point_ptr(glTF_vertices[glTF_v].Geogram_vertex), 12); // write x,y,z from float* getter
        // write texture
        memcpy(buffer_0.data.data() + buffer_texture_coordinates_start + glTF_v*12,  GEO::vec2f((float) glTF_vertices[glTF_v].value,0.5f).data(), 8);
        memcpy(buffer_0.data.data() + buffer_texture_coordinates_start + glTF_v*12 + 8, &zeros_on_4_bytes, 4); // padding
    }

    // write wireframe edges (= vertex indices)
    FOR(e,wireframe_edges_as_vector.size()) {
        memcpy(buffer_0.data.data() + buffer_wireframe_edges_start + 8*e + 0, &(*wireframe_edges_as_vector[e].begin()),  4);
        memcpy(buffer_0.data.data() + buffer_wireframe_edges_start + 8*e + 4, &(*wireframe_edges_as_vector[e].rbegin()), 4);
    }

    // create 3 buffer views
    m.bufferViews.resize(3);
    const size_t BUFFERVIEW_0 = 0; // vertex indices of triangles
    const size_t BUFFERVIEW_1 = 1; // 3D coordinates of vertices + 2D texture coordinates of vertices
    const size_t BUFFERVIEW_2 = 2; // vertex indices of wireframe edges
    tinygltf::BufferView& bufferview_0 = m.bufferViews[BUFFERVIEW_0];
    tinygltf::BufferView& bufferview_1 = m.bufferViews[BUFFERVIEW_1];
    tinygltf::BufferView& bufferview_2 = m.bufferViews[BUFFERVIEW_2];
    bufferview_0.buffer = BUFFER_0;
    bufferview_0.byteOffset = buffer_triangles_start;
    bufferview_0.byteLength = buffer_triangles_length;
    bufferview_0.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
    bufferview_1.buffer = BUFFER_0;
    bufferview_1.byteOffset = buffer_vertices_coordinates_start;
    bufferview_1.byteLength = buffer_vertices_coordinates_length+buffer_texture_coordinates_length;
    bufferview_1.byteStride = 12;
    bufferview_1.target = TINYGLTF_TARGET_ARRAY_BUFFER;
    bufferview_2.buffer = BUFFER_0;
    bufferview_2.byteOffset = buffer_wireframe_edges_start;
    bufferview_2.byteLength = buffer_wireframe_edges_length;
    bufferview_2.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;

    // create 4 accessors
    m.accessors.resize(4);
    const size_t ACCESSOR_0 = 0; // vertex indices of triangles
    const size_t ACCESSOR_1 = 1; // 3D coordinates of vertices
    const size_t ACCESSOR_2 = 2; // 2D texture coordinates of vertices
    const size_t ACCESSOR_3 = 3; // vertex indices of wireframe edges
    tinygltf::Accessor& accessor_0 = m.accessors[ACCESSOR_0];
    tinygltf::Accessor& accessor_1 = m.accessors[ACCESSOR_1];
    tinygltf::Accessor& accessor_2 = m.accessors[ACCESSOR_2];
    tinygltf::Accessor& accessor_3 = m.accessors[ACCESSOR_3];
    accessor_0.bufferView = BUFFERVIEW_0;
    accessor_0.byteOffset = 0;
    accessor_0.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
    accessor_0.count = mesh.facets.nb()*3;
    accessor_0.type = TINYGLTF_TYPE_SCALAR;
    accessor_0.maxValues = {(double) (glTF_vertices.size()-1)};
    accessor_0.minValues = {0};
    accessor_1.bufferView = BUFFERVIEW_1;
    accessor_1.byteOffset = 0;
    accessor_1.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_1.count = glTF_vertices.size();
    accessor_1.type = TINYGLTF_TYPE_VEC3;
    accessor_1.maxValues = bounding_box_max;
    accessor_1.minValues = bounding_box_min;
    accessor_2.bufferView = BUFFERVIEW_1;
    accessor_2.byteOffset = buffer_texture_coordinates_start-bufferview_1.byteOffset; // byte offset relative to the bufferview
    accessor_2.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    accessor_2.count = glTF_vertices.size();
    accessor_2.type = TINYGLTF_TYPE_VEC2;
    accessor_2.maxValues = {SJ_stats.max(), 0.5};
    accessor_2.minValues = {SJ_stats.min(), 0.5};
    accessor_3.bufferView = BUFFERVIEW_2;
    accessor_3.byteOffset = 0;
    accessor_3.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
    accessor_3.count = wireframe_edges_as_vector.size()*2;
    accessor_3.type = TINYGLTF_TYPE_SCALAR;
    accessor_3.maxValues.push_back((double) (glTF_vertices.size()-1));
    accessor_3.minValues.push_back(0);

    // create a mesh with 2 primitives

    m.meshes.resize(1);
    const size_t MESH_0 = 0;
    tinygltf::Mesh& mesh_0 = m.meshes[MESH_0];
    mesh_0.primitives.resize(2);
    const size_t PRIMITIVE_0 = 0;
    const size_t PRIMITIVE_1 = 1;
    tinygltf::Primitive& primitive_0 = mesh_0.primitives[PRIMITIVE_0];
    tinygltf::Primitive& primitive_1 = mesh_0.primitives[PRIMITIVE_1];
    primitive_0.indices = ACCESSOR_0;
    primitive_0.attributes["POSITION"] = ACCESSOR_1;
    primitive_0.attributes["TEXCOORD_0"] = ACCESSOR_2;
    primitive_0.material = MATERIAL_0;
    primitive_0.mode = TINYGLTF_MODE_TRIANGLES; // default value in spec, but tinygltf writes -1 if not specified
    primitive_1.indices = ACCESSOR_3;
    primitive_1.attributes["POSITION"] = ACCESSOR_1;
    primitive_1.material = MATERIAL_1;
    primitive_1.mode = TINYGLTF_MODE_LINE;

    // create a node
    m.nodes.resize(1);
    const size_t NODE_0 = 0;
    m.nodes[NODE_0].mesh = MESH_0; // link NODE_0 to MESH_0

    // create a scene
    m.scenes.resize(1);
    const size_t SCENE_0 = 0;
    m.scenes[SCENE_0].nodes.push_back(NODE_0); // link SCENE_0 to NODE_0

    // describe the asset
    m.asset.version = "2.0"; // required
    m.asset.generator = "tinygltf";

    // 8.

    tinygltf::TinyGLTF gltf;
    fmt::println(Logger::out("glTF"),"Writing {}...",filename); Logger::out("glTF").flush();
    gltf.WriteGltfSceneToFile(
        &m,
        filename,
        true, // embedImages
        true, // embedBuffers
        true, // pretty print
        true  // write binary
    );
}