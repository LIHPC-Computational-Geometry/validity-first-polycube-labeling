# automatic_polycube

Automatic polycube generation for hex-meshing. WIP

## Requirements

- Linux-based OS
- gcc >= 11 (C++20)
- [CMake](https://cmake.org/) > 3.11

## Dependencies

 Name | License | Description | Inclusion
------|---------|-------------|-----------
[geogram](https://github.com/BrunoLevy/geogram) | [BSD-3-Clause](https://github.com/BrunoLevy/geogram/blob/main/LICENSE) | mesh data structure, I/O and GUI | git submodule
[Eigen](https://gitlab.com/libeigen/eigen/) | [Mozilla Public License 2.0](https://gitlab.com/libeigen/eigen/-/blob/master/COPYING.MPL2) | linear algebra, including JacobiSVD solver | [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) (CMake's FetchContent)
[fmt](https://github.com/fmtlib/fmt) | [MIT](https://github.com/fmtlib/fmt/blob/master/LICENSE) | modern string formatting | git submodule
[dbg-macro](https://github.com/sharkdp/dbg-macro) | [MIT](https://github.com/sharkdp/dbg-macro/blob/master/LICENSE) | the must-have dbg(...) macro to replace cout/printf-based debugging | git submodule
[nlohmann/json](https://github.com/nlohmann/json) | [MIT](https://github.com/nlohmann/json/blob/develop/LICENSE.MIT) | JSON for modern C++ | git submodule
[tinygltf](https://github.com/syoyo/tinygltf) | [MIT](https://github.com/syoyo/tinygltf/blob/release/LICENSE) | glTF 2.0 (the "JPEG of 3D") export | git submodule
[DisjointSet](https://www.nayuki.io/page/disjoint-set-data-structure) | [MIT](https://www.nayuki.io/page/disjoint-set-data-structure) | disjoint-set/union-find data structure, implementation by [Nayuki](https://www.nayuki.io/) | in the source code
[gco-v3.0](https://vision.cs.uwaterloo.ca/code/) | [**for research purposes only**, patented](ext/GraphCutOptimization/GCO_README.TXT) | multi-label energies optimization, implementation by Olga Veksler and Andrew Delong | in the source code

## How to build

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
# You may need to specify the platform by adding something like
# '-DVORPALINE_PLATFORM=Linux64-gcc-dynamic' to the next command
# See https://github.com/BrunoLevy/geogram/wiki#compiling
# and https://github.com/BrunoLevy/geogram/blob/main/configure.sh#L120
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
