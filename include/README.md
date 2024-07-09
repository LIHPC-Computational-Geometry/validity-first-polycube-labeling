# Source files

## Containers

Additional functions/macros for containers.

- [`containers_macros.h`](containers_macros.h) : macros for iterators, eg `VECTOR_CONTAINS()`
- [`containers_std.h`](containers_std.h) : for Standard Library containers, eg `std::vector`
- [`containers_Geogram.h`](containers_Geogram.h) : for Geogram containers, eg `GEO::vec3`

## Statistics

- [`stats.h`](stats.h) : compute min/max/sum/avg/sd

## Geometry

Axes, comparison operator for geometric objects, distance...

- [`geometry.h`](geometry.h)
- [`geometry_hexahedra.h`](geometry_hexahedra.h) : geometric quality of hexahedra
- [`geometry_distortion.h`](geometry_distortion.h) : computation of mapping distortion metrics
- [`geometry_halfedges.h`](geometry_halfedges.h) : extended version of Geogram's [MeshHalfedges](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram/mesh/mesh_halfedges.h) class

## Labeling

Generation of / operations on polycube labelings.

- [`labeling.h`](labeling.h) : core labeling definition and functions
- [`labeling_io.h`](labeling_io.h) : read/write labeling files
- [`labeling_graphcuts.h`](labeling_graphcuts.h) : high-level interface of [GraphCutOptimization](../ext/GraphCutOptimization/), for labeling generation
- [`labeling_generators.h`](labeling_generators.h) : functions to generate a labeling
- [`labeling_graph.h`](labeling_graph.h) : computation of charts/boundaries/corners from a labeling
- [`labeling_operators_on_invalidity.h`](labeling_operators_on_invalidity.h) : functions to edit a labeling, targeting a valid one
- [`labeling_operators_on_distortion.h`](labeling_operators_on_distortion.h) : functions to edit a labeling, targeting a low distortion one

See also [`gui_labeling.h`](gui_labeling.h).

## I/O

- [`io_glTF.h`](io_glTF.h) : export a 3D model as glTF
- [`io_dump.h`](io_dump.h) : functions to dump a vertex/edge/boundary/facet as .geogram file

See also [`labeling_io.h`](labeling_io.h).

## GUI

Graphical User Interface classes

- [`SimpleMeshApplicationExt.h`](SimpleMeshApplicationExt.h) : extended version of Geogram's [SimpleMeshApplication](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram_gfx/gui/simple_mesh_application.h) class
- [`gui_labeling.h`](gui_labeling.h) : a GUI app to display labelings on triangle meshes, based on `SimpleMeshApplicationExt`