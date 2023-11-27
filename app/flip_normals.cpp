// Flip the normals of a triangle mesh
// that is, outwards -> inwards or inwards -> outwards (assuming all have the same direction)

#include <geogram/mesh/mesh.h>  // for GEO::Mesh
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>

#include <vector>
#include <string>
#include <set>

using namespace GEO;

int main(int argc, char** argv) {

    std::vector<std::string> filenames;
	GEO::initialize();
	if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh output_mesh"
		))
	{
        fmt::println(Logger::err("I/O"),"Usage should be\n{} input_mesh output_mesh",argv[0]); Logger::err("I/O").flush();
		return 1;
	}

    Mesh M;
    geo_assert(FileSystem::is_file(filenames[0]));

    if(!GEO::mesh_load(filenames[0],M)) {
        fmt::println("Unable to open {}",filenames[0]);
        return 1;
    }
    geo_assert(M.cells.nb() == 0); // must be a surface mesh
    geo_assert(M.facets.are_simplices()); // must be a triangle mesh

    index_t tmp_vertex_index = index_t(-1);
    index_t tmp_facet_index = index_t(-1);
    FOR(f,M.facets.nb()) { // for each facet
        //
        // f = index of current facet
        // lv = local vertex in {0,1,2}
        // le = local edge in {0,1,2}
        // v = vertex index
        // af = adjacent facet
        //
        //                         vA
        //  ------------------------------------------------
        //   \                     / \                     /
        //    \                   /lv0\                   /
        //     \                 /     \                 /
        //      \     afK       /       \      afI      /
        //       \             /         \             /
        //        \           /           \           /
        //         \         /le2       le0\         /
        //          \       /      f        \       /
        //           \     /                 \     /
        //            \   /                   \   /
        //             \ /lv2     le1       lv1\ /
        //           vC ------------------------- vB
        //               \                     /
        //                \                   /
        //                 \                 /
        //                  \               /
        //                   \             /
        //                    \    afJ    /
        //                     \         /
        //                      \       /
        //                       \     /
        //                        \   /
        //                         \ /
        //
        // local vertices are clockwise -> right hand rule -> facet normal is inwards
        // swap lv1 and lv2 to flip the normal, and swap adjacent facet accordingly
        //
        //                         vA
        //  ------------------------------------------------
        //   \                     / \                     /
        //    \                   /lv0\                   /
        //     \                 /     \                 /
        //      \     afK       /       \      afI      /
        //       \             /         \             /
        //        \           /           \           /
        //         \         /le0       le2\         /
        //          \       /      f        \       /
        //           \     /                 \     /
        //            \   /                   \   /
        //             \ /lv1     le1       lv2\ /
        //           vC ------------------------- vB
        //               \                     /
        //                \                   /
        //                 \                 /
        //                  \               /
        //                   \             /
        //                    \    afJ    /
        //                     \         /
        //                      \       /
        //                       \     /
        //                        \   /
        //                         \ /
        //
        // local vertices are counterclockwise -> right hand rule -> facet normal is outwards
        //
        tmp_vertex_index = M.facets.vertex(f,1); // copy vB
        tmp_facet_index = M.facets.adjacent(f,0); // copy afI
        M.facets.set_vertex(f,1,M.facets.vertex(f,2));// at lv1 is no longer vB but vC
        M.facets.set_adjacent(f,0,M.facets.adjacent(f,2)); // beyond le0 is no longer afI but afK
        M.facets.set_vertex(f,2,tmp_vertex_index); // at lv2 is no longer vC but vB
        M.facets.set_adjacent(f,2,tmp_facet_index); // beyond le2 is no longer afK but afI
    }

    mesh_save(M,filenames[1]);

    return 0;
}