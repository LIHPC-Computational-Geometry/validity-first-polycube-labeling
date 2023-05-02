#include <fmt/core.h>
#include <fmt/ostream.h>

#include <iostream>

#include "HexMeshViewerApp.h"

int main(int argc, char** argv) {
    
	if(argc != 2) {
		fmt::println(std::cerr,"Usage should be ./hex_mesh_viewer hexmesh.mesh");
		exit(1);
	}
	
    HexMeshViewerApp app;
	app.start(argc,argv);
    return 0;
}