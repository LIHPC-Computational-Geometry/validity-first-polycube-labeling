# automatic_polycube
Automatic polycube generation for hex-meshing. WIP

## How to build

```bash
git clone --recurse-submodules https://github.com/LIHPC-Computational-Geometry/automatic_polycube.git
cd automatic_polycube
./configure.sh
cd build/Linux64-gcc-dynamic-Release/
make
```

## App `labeling_viewer`

```bash
# from build/Linux64-gcc-dynamic-Release
cd bin
./labeling_viewer ../../../data/B0/surface.obj ../../../data/B0/labeling.txt
```

## App `hex_mesh_viewer`

```bash
# from build/Linux64-gcc-dynamic-Release
cd bin
./hex_mesh_viewer ../../../data/B0/hex.mesh
```