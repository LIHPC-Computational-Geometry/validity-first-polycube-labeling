# automatic_polycube
Automatic polycube generation for hex-meshing. WIP

## How to build

```bash
git clone --recurse-submodules https://github.com/LIHPC-Computational-Geometry/automatic_polycube.git
cd automatic_polycube
```

Open `ext/geogram/src/lib/geogram/mesh/mesh_halfedges.h` and change accessibility of `MeshHalfedges` variables from `private` to `protected`.

```diff
         void move_to_opposite(Halfedge& H) const;
 
-    private:
+    protected:
         Mesh& mesh_;
         Attribute<index_t> facet_region_;
```

Then:

```bash
# from automatic_polycube/
./configure.sh
cd build/Linux64-gcc-dynamic-Release/
make automatic_polycube
```

## Run `automatic_polycube` app

```bash
# from build/Linux64-gcc-dynamic-Release
./bin/automatic_polycube ../../data/B0/surface.obj
```

1. In the left panel, click on "Compute naive labeling"
1. If the labeling in invalid*, try "Auto fix validity", which loop over manual operators just over the button
1. If not all boundaries are monotone*, try "Auto fix monotonicity" and take a coffee. Turning points should now be removed.
1. In the menu bar, click on "Save as", choose a location and a filename, and select "txt" as extension to export the per-surface-triangle labeling.

*Not the case with the provided B0 model.

## Other applications

This project contains other applications designed for [HexMeshWorkshop](https://github.com/LIHPC-Computational-Geometry/HexMeshWorkshop). The following ones may be of interest here:
- `labeling_viewer` : like `automatic_polycube` but without labeling operators, for the purpose of visualizing existing labeling files (you can import the mesh & the labeling from the command line, or drag-and-drop them one by one)
- `labeling_painter` : an interactive app to manually transform a labeling with pencil & bucket fill tools