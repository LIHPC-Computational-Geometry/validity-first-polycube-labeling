#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/command_line_args.h>	// for CmdLine::import_arg_group()
#include <geogram/basic/command_line.h>			// for CmdLine::parse()

#include "LabelingViewerApp.h"
#include "labeling.h"

#include <vector>
#include <string>

using namespace GEO;

int main(int argc, char** argv) {

    std::vector<std::string> filenames;
	initialize();
    CmdLine::import_arg_group("sys"); // declares sys:compression_level, needed by mesh_save() for .geogram files
	if(!CmdLine::parse(
		argc,
		argv,
		filenames,
		"<input_surface_mesh> <input_surface_labeling> <output_geogram_file>"
		))
	{
		return 1;
	}

    if(filenames.size()<3) { // if <output_geogram_file> is not provided -> launch GUI (default)
        LabelingViewerApp app;
        app.start(argc,argv);
        return 0;
    }
    // else : write a .geogram mesh with the labeling as attribute, to be visualized with Graphite

    Mesh M;
	if(!mesh_load(filenames[0],M)) {
		fmt::println(Logger::err("I/O"),"Unable to load mesh from {}",filenames[0]);
		return 1;
	}

    load_labeling(filenames[1],M,LABELING_ATTRIBUTE_NAME);

    if(!mesh_save(M,filenames[2])) {
		fmt::println(Logger::err("I/O"),"Unable to save mesh to {}",filenames[2]);
		return 1;
	}

    return 0;
}