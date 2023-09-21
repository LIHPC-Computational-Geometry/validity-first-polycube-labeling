// OpenVolumeMesh ASCII format: (tet-mesh)
// https://www.graphics.rwth-aachen.de/media/openvolumemesh_static/Documentation/OpenVolumeMesh-Doc-Latest/file_format.html
// To my understanding, see below
// ------------------
// | OVM ASCII
// | Vertices
// | <nb_vertices>
// | <vertex[0].x> <vertex[0].y> <vertex[0].z>
// | ...
// | Edges
// | <nb_edges>
// | <edge[0].vertex0> <edge[0].vertex1>
// | Faces
// | 3 <face[0].edge0> <face[0].edge1> <face[0].edge2>
// | ...
// | Polyhedra
// | <nb_polyhedra>
// | 4 <polyhedron[0].face0> <polyhedron[0].face1> <polyhedron[0].face2> <polyhedron[0].face3>
// | ...
// ------------------

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>

using namespace GEO;

int main(int argc, const char** argv) {

    if (argc<2) {
        fmt::println("Missing arguments");
        fmt::println("./extract_surface <input.ext/.ovm> <output.ovm/.mesh>");
        return 1;
    }

    GEO::initialize();

    std::string input = argv[1];
    std::string output = argv[2];

    geo_assert(FileSystem::is_file(input));

    if ( (FileSystem::extension(input)  != "ovm") &&
         (FileSystem::extension(output) == "ovm") ) {

        fmt::println(Logger::out("I/O"),"Converting to OpenVolumeMesh"); Logger::out("I/O").flush();

        Mesh input_mesh;
        if(!GEO::mesh_load(input,input_mesh)) {
            fmt::println(Logger::err("I/O"),"Unable to open {}",input); Logger::err("I/O").flush();
            return 1;
        }

        geo_assert(input_mesh.cells.nb() != 0); // must be a volumetric mesh
        geo_assert(input_mesh.cells.are_simplices()); // must be a tetrahedral mesh

        // open output file, create variables
        auto out = fmt::output_file(output);
        index_t max_elem = 0;
        // index_t current_index = 0;
        vec3 coordinates = {0.0,0.0,0.0};

        // write header
        out.print("OVM ASCII\n");

        // write vertices
        max_elem = input_mesh.vertices.nb();
        out.print("Vertices\n");
        out.print("{}\n",max_elem);
        FOR(v,max_elem) {
            coordinates = input_mesh.vertices.point(v);
            out.print("{} {} {}\n",coordinates.x, coordinates.y, coordinates.z);
        }

        // write edges
        max_elem = input_mesh.edges.nb();
        out.print("Edges\n");
        out.print("{}\n",max_elem);
        FOR(e,max_elem) {
            out.print("{} {}\n",input_mesh.edges.vertex(e,0), input_mesh.edges.vertex(e,1));
        }

        // write faces
        /*max_elem = input_mesh.facets.nb();
        out.print("Faces\n");
        out.print("{}\n",max_elem);
        FOR(f,max_elem) {
            out.print("3 {} {} {}\n", ??? ); // TODO get edges around current facet
        }*/

        // TODO write cells

        out.close();

    }
    else if ( (FileSystem::extension(input)  == "ovm") &&
              (FileSystem::extension(output) != "ovm") ) {
        fmt::println(Logger::out("I/O"),"Converting from OpenVolumeMesh"); Logger::out("I/O").flush();
        fmt::println(Logger::err("I/O"),"Not implemented"); Logger::err("I/O").flush();
        return 1;
    }
    else {
        fmt::println(Logger::err("I/O"),"Neither the input nor the ouput is a .ovm"); Logger::err("I/O").flush();
        return 1;
    }

    return 0;
}