# automatic_polycube
Automatic polycube generation for hex-meshing. WIP

## How to build

```bash
git clone --recurse-submodules https://github.com/LIHPC-Computational-Geometry/automatic_polycube.git
cd automatic_polycube
```

Open `ext/geogram/src/lib/geogram/mesh/mesh_halfedges.h` and change accessibility of `MeshHalfedges` variables from `private` to `protected`. Then:

```bash
# from automatic_polycube/
./configure.sh
cd build/Linux64-gcc-dynamic-Release/
make
```

## Run `automatic_polycube` app

```bash
# from build/Linux64-gcc-dynamic-Release
./bin/automatic_polycube ../../data/B0/surface.obj
```

## Run `hex_mesh_viewer` app

```bash
# from build/Linux64-gcc-dynamic-Release
./bin/hex_mesh_viewer ../../data/B0/hex.mesh
```