# Source files

- **[containers.h](containers.h)/[.cpp](../src/containers.cpp)** : additional functions & macros for std containers (contains, min, concatenation...)
- **[CustomMeshHalfedges.h](CustomMeshHalfedges.h)** : modified version of Geogram's [MeshHalfedges](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram/mesh/mesh_halfedges.h) class
- **[geometry.h](geometry.h)** : geometry-related functions (labels & axes, comparison operator for geometric objects)
- **[GraphCutLabeling.h](GraphCutLabeling.h)/[.cpp](../src/GraphCutLabeling.cpp)** : constrained labeling generation with a Graph-Cut Optimization
- **[hex_mesh.h](hex_mesh.h)/[.cpp](../src/hex_mesh.cpp)** : hexahedral meshes related functions (quality metric)
- **[labeling.h](labeling.h)/[.cpp](../src/labeling.cpp)** : operators on labelings
- **[LabelingGraph.h](LabelingGraph.h)/[.cpp](../src/LabelingGraph.cpp)** : computation of charts/boundaries/corners from a labeling
- **[SimpleMeshApplicationExt.h](SimpleMeshApplicationExt.h)** : modified version of Geogram's [SimpleMeshApplication](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram_gfx/gui/simple_mesh_application.h) class