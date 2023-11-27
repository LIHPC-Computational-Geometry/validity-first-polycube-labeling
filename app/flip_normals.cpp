// Flip the normals of a triangle mesh
// that is, outwards -> inwards or inwards -> outwards (assuming all have the same direction)
// If no output mesh filename is provided, determine the normals direction

#include <geogram/mesh/mesh.h>  // for GEO::Mesh
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/numeric.h>              // for max_float64()

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>

#include <vector>
#include <string>
#include <set>

#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

using namespace GEO;

int main(int argc, char** argv) {

    std::vector<std::string> filenames;
	GEO::initialize();
	if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"input_mesh <output_mesh>"
		))
	{
        fmt::println(Logger::err("I/O"),"Usage should be\n{} input_mesh <output_mesh>",argv[0]); Logger::err("I/O").flush();
		return 1;
	}

    Mesh M;
    geo_assert(FileSystem::is_file(filenames[0]));

    if(!GEO::mesh_load(filenames[0],M)) {
        fmt::println("Unable to open {}",filenames[0]);
        return 1;
    }
    geo_assert(M.cells.nb() == 0); // must be a surface mesh
    geo_assert(M.facets.nb() != 0); // must have facets
    geo_assert(M.facets.are_simplices()); // must be a triangle mesh

    if(filenames.size()==1) {
        // no output mesh -> determine normals direction
        // https://forums.cgsociety.org/t/check-if-mesh-is-inside-out/1688290

        // Find the vertex with the smallest X coordinate
        double smallest_X_coordinate = Numeric::max_float64();
        index_t corresponding_vertex = index_t(-1);
        FOR(v,M.vertices.nb()) { // for each vertex
            if(smallest_X_coordinate > M.vertices.point(v).x) {
                smallest_X_coordinate = M.vertices.point(v).x;
                corresponding_vertex = v;
            }
        }

        // put a point outside the mesh
        vec3 point_outside_mesh = M.vertices.point(corresponding_vertex) - vec3(1.0,0.0,0.0);

        // find a facet adjacent to the corresponding_vertex
        // and create a halfedge whose origin is corresponding_vertex
        // Is there a easier way that going through all facets?
        MeshHalfedges::Halfedge halfedge;
        CustomMeshHalfedges mesh_he(M);
        FOR(f,M.facets.nb()) { // for each facet
            if(
            (M.facets.vertex(f,0) == corresponding_vertex) ||
            (M.facets.vertex(f,1) == corresponding_vertex) ||
            (M.facets.vertex(f,2) == corresponding_vertex) ) {
                halfedge.facet = f;
                break;
            }
        }
        halfedge.corner = M.facets.corner(halfedge.facet,0); // try the halfedge at local vertex 0
        while (Geom::halfedge_vertex_index_from(M,halfedge) != corresponding_vertex) {
            mesh_he.move_to_next_around_facet(halfedge);
        }

        // compute the vertex normal of the corresponding_vertex
        vec3 normal_of_corresponding_vertex;
        MeshHalfedges::Halfedge init_halfedge = halfedge;
        do {
            normal_of_corresponding_vertex += mesh_facet_normal(M,halfedge.facet);
            mesh_he.move_to_next_around_vertex(halfedge);
        } while (halfedge != init_halfedge);

        // compute the dot product of the vertex normal and a vector going outwards the mesh (towards -X)
        double dot_product = dot(normal_of_corresponding_vertex,point_outside_mesh-M.vertices.point(corresponding_vertex));

        if(dot_product > 0) {
            fmt::println(Logger::out("normals dir."),"The facet normals are outwards"); Logger::out("normals dir.").flush();
        }
        else {
            fmt::println(Logger::out("normals dir."),"The facet normals are inwards"); Logger::out("normals dir.").flush();
        }
        
        return 0;
    }

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