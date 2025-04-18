# Applications

## Polycube labeling generation

From a 3D triangle mesh, generate a labeling which maps each triangle to one of the 6 labels $\{0 = +X, 1 = -X, 2 = +Y, 3 = -Y, 4 = +Z, 5 = -Z\}$.

- [`naive_labeling.cpp`](naive_labeling.cpp) : pick the nearest label of the triangle's normal
- [`graphcut_labeling.cpp`](graphcut_labeling.cpp) : graph-cuts optimization
- [`automatic_polycube.cpp`](automatic_polycube.cpp) : main contribution, an automatic pipeline targeting labeling validity before low mapping distortion
- [`labeling_painter.cpp`](labeling_painter.cpp) : manual editing of a labeling

## Viewers

- [`labeling_viewer.cpp`](labeling_viewer.cpp) : visualize a colored triangle mesh based on a labeling
- [`hex_mesh_viewer.cpp`](hex_mesh_viewer.cpp) : visualize a colored hexahedral mesh based on the cells' quality

## Mesh transformations

- [`extract_surface.cpp`](extract_surface.cpp) : extract the surface triangle mesh from a tetrahedral mesh
- [`flip_normals.cpp`](flip_normals.cpp) : flip all the triangles, so that their normals are in the other direction
- [`reorient_mesh.cpp`](reorient_mesh.cpp) : rotate a mesh according to the principal axes of the point cloud
- [`convert_labeled_obj.cpp`](convert_labeled_obj.cpp) : extract the labeling (as .txt) and the un-labeled mesh (as .obj) from a labeled .obj (PolyCut outputs)

## Stats and metrics computation

- [`labeling_stats.cpp`](labeling_stats.cpp) : from a triangle mesh and a labeling, write the labeling graph details
- [`mesh_stats.cpp`](mesh_stats.cpp) : write vertices/facets/cells stats of a given mesh
- [`polycube_distortion.cpp`](polycube_distortion.cpp) : write distortion metrics of the mapping between a triangle mesh and a polycube mesh

## Export

- [`to_glTF.cpp`](to_glTF.cpp) : write a 3D render (triangle mesh, labeled triangle mesh, hexahedral mesh) as glTF 2.0
- [`volume_labeling.cpp`](volume_labeling.cpp) : export a per-tetrahedron-facet labeling from a per-surface-triangle labeling and a triangle-to-tetrahedron-facet map

## Debug

- [`areas_to_tilt.cpp`](areas_to_tilt.cpp) : to tweak the sensitivity for [`../include/labeling_generators.h`](../include/labeling_generators.h) > `naive_labeling()` & `graphcut_labeling()`
- [`halfedge_movements.cpp`](halfedge_movements.cpp) : to understand the half-edges API of Geogram