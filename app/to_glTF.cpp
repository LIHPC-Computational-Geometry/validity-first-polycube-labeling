// Define these only in *one* .cc file.
// from ext/tinygltf/README.md > Loading glTF 2.0 model
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <geogram/basic/command_line.h>			// for CmdLine::parse()
#include <geogram/basic/command_line_args.h>    // for CmdLine::import_arg_group()
#include <geogram/basic/logger.h>               // for GEO::Logger::*
#include <geogram/mesh/mesh.h>                  // for Mesh
#include <geogram/mesh/mesh_io.h>               // for mesh_load()

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <vector>

#include "geometry_halfedges.h"    // for MeshHalfedgesExt
#include "io_glTF.h"
#include "labeling.h"
#include "labeling_io.h"

// To see the inside of some 3D model, remove specific facets
// #define REMOVE_FACETS_FOR_IN_VOLUME_TWIST_MODEL
// #define REMOVE_FACETS_FOR_IN_VOLUME_KNOT_MODEL

using GEO::index_t; // to use the FOR() macro of Geogram

int main(int argc, char **argv)
{
    GEO::initialize();
    CmdLine::import_arg_group("sys"); // to export .geogram files

    CmdLine::declare_arg(
        "labeling",
        "",
        "path to the labeling file"
    );

    CmdLine::declare_arg(
        "polycube",
        "",
        "path to the polycube mesh"
    );

    std::vector<std::string> filenames;
    if(!CmdLine::parse(argc,argv,filenames,"input_mesh output_filename"))	{
		return 1;
	}
    std::string labeling_filename = GEO::CmdLine::get_arg("labeling");
    std::string polycube_filename = GEO::CmdLine::get_arg("polycube");

    GEO::Mesh M;
    MeshExt M_ext(M);
    mesh_load(filenames[0],M);

    if(M.cells.nb() != 0) {
        // The user provided a volume mesh
        // Expect an hex-mesh
        geo_assert(!M.cells.are_simplices());
        FOR(c,M.cells.nb()) {
            geo_assert(M.cells.type(c) == MESH_HEX);
        }
        geo_assert(labeling_filename.empty()); // the user should not provide a labeling in this case
        write_glTF__hex_mesh_surface(filenames[1],M);
        return 0;
    }

    // ensure facet normals are outward
	if(facet_normals_are_inward(M)) {
		flip_facet_normals(M);
		fmt::println(GEO::Logger::warn("normals dir."),"Facet normals of the input mesh were inward");
		fmt::println(GEO::Logger::warn("normals dir."),"You should flip them and update the file for consistency.");
		GEO::Logger::warn("normals dir.").flush();
	}

    if(labeling_filename.empty()) {
        write_glTF__triangle_mesh(filenames[1],M_ext,true);
        return 0;
	}
    //else: the user provided a labeling

    // expand '~' to $HOME
    if(labeling_filename[0] == '~') {
        labeling_filename.replace(0, 1, std::string(getenv("HOME")));
    }

    Attribute<index_t> labeling(M.facets.attributes(),LABELING_ATTRIBUTE_NAME);
    load_labeling(labeling_filename,M,labeling);

    M_ext.adj_facet_corners.recompute();

    if(polycube_filename.empty()) {

        auto remove_facet_fx = [](const Mesh& M, const Attribute<index_t>& label, index_t facet_index) { // return true if the facet at `facet_index` should be remove before glTF export
        #ifdef REMOVE_FACETS_FOR_IN_VOLUME_TWIST_MODEL
            // Delete some facets for the new 'in-volume_twist' model (https://github.com/LIHPC-Computational-Geometry/nightmare_of_polycubes)
            // Facets whose:
            // - label is +Z (=4)
            // and
            // - vertices are near 9000 in the Z axis
            const double Z_axis_threshold = 9000*0.99; // = remove 1% to the delta Z
            return (
                (label[facet_index] == 4) &&
                (mesh_vertex(M,M.facets.vertex(facet_index,0)).z > Z_axis_threshold) && 
                (mesh_vertex(M,M.facets.vertex(facet_index,1)).z > Z_axis_threshold) && 
                (mesh_vertex(M,M.facets.vertex(facet_index,2)).z > Z_axis_threshold)
            );
        #elif REMOVE_FACETS_FOR_IN_VOLUME_KNOT_MODEL
            // Delete some facets for the new 'in-volume_knot' model (https://github.com/LIHPC-Computational-Geometry/nightmare_of_polycubes)
            // TODO
        #else
            geo_argused(M);
            geo_argused(label);
            geo_argused(facet_index);
            return false; // remove no facets
        #endif
        };
        
        write_glTF__labeled_triangle_mesh(filenames[1],M_ext,labeling,remove_facet_fx);
        return 0;
    }
    //else: the user provided a polycube

    // expand '~' to $HOME
    if(polycube_filename[0] == '~') {
        polycube_filename.replace(0, 1, std::string(getenv("HOME")));
    }

    GEO::Mesh polycube;
    mesh_load(polycube_filename,polycube);

    write_glTF__labeled_triangle_mesh_with_polycube_animation(filenames[1],M_ext,polycube,labeling);

    return 0;
}