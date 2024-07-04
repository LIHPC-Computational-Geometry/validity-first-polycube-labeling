# Applications

## Polycube labeling generation

From a 3D triangle mesh, generate a labeling which maps each triangle to one of the 6 labels $\{0 = +X, 1 = -X, 2 = +Y, 3 = -Y, 4 = +Z, 5 = -Z\}$.

- `naive_labeling.cpp` : pick the nearest label of the triangle's normal
- `graphcut_labeling.cpp` : graph-cuts optimization
- `automatic_polycube.cpp` : main contribution, an automatic pipeline targeting labeling validity before low mapping distortion
- `labeling_painter.cpp` : manual editing of a labeling

## Viewers

- `labeling_viewer.cpp` : visualize a colored triangle mesh based on a labeling
- `hex_mesh_viewer.cpp` : visualize a colored hexahedral mesh based on the cells' quality

## Mesh transformations

- `extract_surface.cpp` : extract the surface triangle mesh from a tetrahedral mesh
- `flip_normals.cpp` : flip all the triangles, so that their normals are in the other direction
- `reorient_mesh.cpp` : rotate a mesh according to the principal axes of the point cloud

## Stats and metrics computation

- `labeling_stats.cpp` : from a triangle mesh and a labeling, write the labeling graph details
- `mesh_stats.cpp` : write vertices/facets/cells stats of a given mesh
- `polycube_distortion.cpp` : write distortion metrics of the mapping between a triangle mesh and a polycube mesh

## Export

- `to_glTF` : write a 3D render (triangle mesh, labeled triangle mesh, hexahedral mesh) as glTF 2.0

## Debug

- `areas_to_tilt.cpp` : to analyse sensitivity for `include/labeling_generators.h` > `tweaked_naive_labeling()`
- `halfedge_movements.cpp` : to understand the half-edges API of Geogram