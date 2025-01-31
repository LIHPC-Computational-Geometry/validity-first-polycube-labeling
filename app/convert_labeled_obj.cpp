// PolyCut (http://www.cs.ubc.ca/labs/imager/tr/2013/polycut/)
// does not export a .txt file for the output labeling but a colored .obj
// Following Evocube (see https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/app/synthesize_data.cpp#L139)
// we parse the .obj for the `usemtl` keywords and map {"red", "darkred", "green", "darkgreen", "blue", "darkblue"} to [0:5]
// (`usemtl` reference material names defined in an external .mtl file, see https://en.wikipedia.org/wiki/Wavefront_.obj_file#Reference_materials)
// and assemble a GEO::Mesh manually (because GEO::mesh_load() does not like .obj with vertex texture coordinates & vertex normals attached to facets https://en.wikipedia.org/wiki/Wavefront_.obj_file#Face_elements)
// to export a color-less .obj
//
// Usage:
//   ./bin/convert_labeled_obj input.obj output.txt
// or
//   ./bin/convert_labeled_obj input.obj output.txt colorless_obj=output.obj
//
// Input format:
// mtllib segmentation.mtl                // Reference the materials definition
// v -5.000000 5.000000 4.052147          // First vertex (xyz coordinates)
// v 5.000000 4.584275 -0.006765          // ...
// ...                                    // ...
// v 5.000000 4.235790 -0.751739          // Last vertex
// vn 1.000000 0.000000 0.000000          // First vertex normal (xyz vector)
// vn 1.000000 0.000000 0.000000          // ...
// ...                                    // ...
// vn 1.000000 0.000000 0.000000          // Last vertex normal
// usemtl green                           // Material of the first facet
// f 13/13/13 1/1/1 54/54/54              // First facet (v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 with vt for textures and vn for normals, we only read v1, v2 and v3 the vertex indices)
// usemtl darkred                         // ...
// f 1/1/1 13/13/13 2761/2761/2761        // ...
// ...                                    // ...
// usemtl blue                            // ...
// f 54/54/54 1/1/1 1304/1304/1304        // Last facet

#include <geogram/mesh/mesh.h> // for Mesh
#include <geogram/mesh/mesh_io.h> // for mesh_load(), mesh_save()
#include <geogram/basic/command_line.h> // for CmdLine::parse()
#include <geogram/basic/file_system.h> // for FileSystem::is_file()
#include <geogram/basic/string.h> // for String::split_string()

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/os.h> // for fmt::output_file()

#include <vector>
#include <map>
#include <fstream>
#include <string>

using namespace GEO;

const std::map<std::string,index_t> MATERIAL_TO_LABEL = {
    {"red",         0}, // +X
    {"darkred",     1}, // -X
    {"green",       2}, // +Y
    {"darkgreen",   3}, // -Y
    {"blue",        4}, // +Z
    {"darkblue",    5}  // -Z
};

int main(int argc, char** argv) {
    GEO::initialize();

    CmdLine::declare_arg(
        "colorless_obj",
        "",
        "path to the output, color-less, .obj mesh"
    );

    std::vector<std::string> filenames;
    if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_obj output_labeling.txt"
    )) {
		return 1;
	}

    std::string input_obj_filepath = filenames[0];
    std::string output_labeling_filepath = filenames[1];
    std::string colorless_obj_filepath = CmdLine::get_arg("colorless_obj");

    if(!FileSystem::is_file(input_obj_filepath)) {
        fmt::println(Logger::err("I/O"),"{} does not exist",input_obj_filepath); Logger::err("I/O").flush();
        return 1;
    }

    Mesh M;
    M.vertices.set_dimension(3);
    M.vertices.set_double_precision();
    vec3 coordinates;
    vec3i indices;

    std::fstream file;
    std::string word;
    std::vector<std::string> word_split;
    std::vector<index_t> labeling;
    file.open(input_obj_filepath.c_str());
    if(!file.good()) {
        fmt::println(Logger::err("I/O"),"Cannot open {}",input_obj_filepath); Logger::err("I/O").flush();
        return 1;
    }
    while (file >> word) { // parse the file word by word
        if (word == "mtllib") {
            file >> word; // read the next word = the materials file & ignore it
        }
        else if (word == "v") {
            file >> coordinates.x;
            file >> coordinates.y;
            file >> coordinates.z;
            M.vertices.create_vertex(coordinates.data());
        }
        else if (word == "vn") {
            file >> coordinates.x;
            file >> coordinates.y;
            file >> coordinates.z;
            // do not store vertex normals
        }
        else if (word == "usemtl") { // found a reference to a material
            file >> word; // read the next word = the material name
            labeling.push_back(MATERIAL_TO_LABEL.at(word)); // store the label in [0:5] corresponding to this material
        }
        else if (word == "f") {
            // read first vertex index of this facet
            file >> word;
            String::split_string(word,'/',word_split);
            indices.x = std::stoi(word_split[0])-1;
            // ignore word_split[1] which is a texture coordinate
            // ignore word_split[2] which is a vertex normal
            word_split.clear();
            // read second vertex index of this facet
            file >> word;
            String::split_string(word,'/',word_split);
            indices.y = std::stoi(word_split[0])-1;
            // ignore word_split[1] which is a texture coordinate
            // ignore word_split[2] which is a vertex normal
            word_split.clear();
            // read third vertex index of this facet
            file >> word;
            String::split_string(word,'/',word_split);
            indices.z = std::stoi(word_split[0])-1;
            // ignore word_split[1] which is a texture coordinate
            // ignore word_split[2] which is a vertex normal
            word_split.clear();
            // store facet in the GEO::Mesh
            M.facets.create_triangle(
                (index_t) indices.x,
                (index_t) indices.y,
                (index_t) indices.z
            );
        }
        else {
            fmt::println(Logger::err("I/O"),"Unexpected word '{}' in {}",word,input_obj_filepath); Logger::err("I/O").flush();
        }
    }
    file.close();

    fmt::println(Logger::out("I/O"),"Found {} vertices",M.vertices.nb()); Logger::out("I/O").flush();
    fmt::println(Logger::out("I/O"),"Found {} triangles",M.facets.nb()); Logger::out("I/O").flush();
    fmt::println(Logger::out("I/O"),"Found {} labels",labeling.size()); Logger::out("I/O").flush();

    geo_assert(M.facets.nb() == labeling.size());

    if(!colorless_obj_filepath.empty()) {
        // expand '~' to $HOME
        if(colorless_obj_filepath[0] == '~') {
            colorless_obj_filepath.replace(0, 1, std::string(getenv("HOME")));
        }
        // write the color-less .obj
        if(!mesh_save(M,colorless_obj_filepath)) {
            fmt::println(Logger::err("I/O"),"Unable to write {}",colorless_obj_filepath); Logger::err("I/O").flush();
        }
    }

    // write labeling as .txt
    fmt::println(Logger::out("I/O"),"Saving file {}...",output_labeling_filepath); Logger::out("I/O").flush();
    auto out = fmt::output_file(output_labeling_filepath);
    FOR(f,labeling.size()) { // for each facet (= label)
        out.print("{}\n",labeling[f]);
    }
}