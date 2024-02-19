# automatic_polycube
Automatic polycube generation for hex-meshing. WIP

## Requirements

- Linux-based OS
- gcc >= 11 (other C++20 compilers should work)
- [CMake](https://cmake.org/) > 3.11
- [geogram](https://github.com/BrunoLevy/geogram) : mesh data structure, I/O and GUI (included as submodule)
- [fmt](https://github.com/fmtlib/fmt) : modern string formatting (included as submodule)
- [dbg-macro](https://github.com/sharkdp/dbg-macro) : the must-have dbg(...) macro to replace cout/printf-based debugging (included as submodule)
- [nlohmann/json](https://github.com/nlohmann/json) : JSON for modern C++ (included as submodule)
- [DisjointSet](https://www.nayuki.io/page/disjoint-set-data-structure) : disjoint-set/union-find data structure, implementation by [Nayuki](https://www.nayuki.io/) (in the source code)
- [gco-v3.0](https://vision.cs.uwaterloo.ca/code/) :  multi-label energies optimization, implementation by Olga Veksler and Andrew Delong, **for research purposes only** (in the source code)

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
mkdir build_Release
cd build_Release
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja automatic_polycube
```

Remove `-G Ninja` to use Make instead, if you don't have Ninja.

## Run `automatic_polycube` app

### With the GUI

```bash
# from automatic_polycube/build_Release
./bin/automatic_polycube ../../data/B0/surface.obj gui=true
```

You can ommit the file path to drag-and-drop the input file.

1. In the right panel, click on "Compute naive labeling"
1. If the labeling in invalid*, try "Auto fix validity", which loop over manual operators just over the button
1. If not all boundaries are monotone*, try "Auto fix monotonicity", which loop over some of the operators listed above.
1. In the menu bar, click on "Save as", choose a location and a filename, and select "txt" as extension to export the per-surface-triangle labeling.

*Not the case with the provided B0 model.

### Without the GUI

```bash
# from automatic_polycube/build_Release
./bin/automatic_polycube ../../data/B0/surface.obj ../../data/B0/labeling.txt gui=false
```

## Other applications

This project contains other applications designed for [HexMeshWorkshop](https://github.com/LIHPC-Computational-Geometry/HexMeshWorkshop). The following ones may be of interest here:
- `labeling_viewer` : like `automatic_polycube` but without labeling operators, for the purpose of visualizing existing labeling files (you can import the mesh & the labeling from the command line, or drag-and-drop them one by one)
- `labeling_painter` : an interactive app to manually transform a labeling with pencil & bucket fill tools
- `graphcut_labeling` : labeling generation from Graph-Cut optimization

## Wiki

[Centralize within the HexMeshWorkshop repo](https://github.com/LIHPC-Computational-Geometry/HexMeshWorkshop/wiki)