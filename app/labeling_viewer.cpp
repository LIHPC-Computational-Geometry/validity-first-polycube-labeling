#include <iostream>

#include "LabelingViewerApp.h"

int main(int argc, char** argv) {
    
	if(argc != 3) {
		std::cerr << "Usage should be ./labeling_viewer surface.obj labeling.txt\n";
		exit(1);
	}
	
    LabelingViewerApp app;
	app.start(argc,argv);
    return 0;
}