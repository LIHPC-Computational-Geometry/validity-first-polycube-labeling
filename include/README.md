# Source files

## Containers

Additional functions/macros for containers.

- `containers_macros.h` : macros for iterators, eg `VECTOR_CONTAINS()`
- `containers_std.h` : for Standard Library containers, eg `std::vector`
- `containers_Geogram.h` : for Geogram containers, eg `GEO::vec3`

## Statistics

- `basic_stats.h` : compute min/max/sum/avg/sd

## Geometry

Axes, comparison operator for geometric objects, distance...

- `geometry.h`
- `geometry_hexahedra.h` : geometric quality of hexahedra
- `CustomMeshHalfedges.h` : extended version of Geogram's [MeshHalfedges](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram/mesh/mesh_halfedges.h) class

## Labeling

Generation of / operations on polycube labelings.

- `labeling.h` : core labeling definition and functions
- `labeling_io.h` : read/write labeling files
- `labeling_graphcuts.h` : high-level interface of [GraphCutOptimization](../ext/GraphCutOptimization/), for labeling generation
- `labeling_generators.h` : functions to generate a labeling
- `labeling_operators_on_invalidity.h` : functions to edit a labeling, targeting a valid one
- `labeling_operators_on_distortion.h` : functions to edit a labeling, targeting a low distortion one
- `LabelingGraph.h` : computation of charts/boundaries/corners from a labeling

## I/O

- `glTF.h` : export a 3D model as glTF
- `dump_mesh.h` : functions to easily write a .geogram containing a vertex/edge/boundary/facet

See also `labeling_io.h`.

## GUI

Graphical User Interface classes

- `SimpleMeshApplicationExt.h` : extended version of Geogram's [SimpleMeshApplication](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram_gfx/gui/simple_mesh_application.h) class
- `LabelingViewerApp.h` : a GUI app to display labelings on triangle meshes, based on `SimpleMeshApplicationExt`