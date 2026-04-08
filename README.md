# bromesh

A C++20 static library for mesh generation, manipulation, and I/O. Designed for game engines, tools, and procedural content pipelines.

## Features

| Category | Algorithms |
|---|---|
| **Isosurface** | Marching cubes, surface nets, dual contouring (QEF), transvoxel (seamless LOD) |
| **Voxel** | Greedy meshing with palette colors |
| **Primitives** | Box, sphere, cylinder, capsule, plane, torus, heightmap grid |
| **Manipulation** | Smooth normals, flat normals, tangents, simplify (quadric error), LOD chain, weld, split components |
| **Analysis** | Bounding box, manifold check, volume, convex decomposition (V-HACD), convex hull |
| **UV** | Box, planar (XY/XZ/YZ), cylindrical, spherical projection |
| **Optimization** | Vertex cache, vertex fetch, overdraw (via meshoptimizer) |
| **I/O** | OBJ read/write, STL read/write, glTF/GLB read/write, MagicaVoxel VOX read |

All algorithms produce `bromesh::MeshData` -- a flat struct with separate position, normal, UV, color, and index arrays ready for GPU upload or TypedArray transfer.

## Building

```bash
git clone --recurse-submodules https://github.com/user/bromesh.git
cd bromesh
cmake -B build
cmake --build build --config Release
```

Requires CMake 3.24+ and a C++20 compiler (MSVC 2022, GCC 12+, Clang 15+).

### Running tests

```bash
cd build
ctest --build-config Release --output-on-failure
```

## Third-party dependencies

All vendored under `third_party/` as git submodules or single headers:

| Library | Purpose | License |
|---|---|---|
| [meshoptimizer](https://github.com/zeux/meshoptimizer) | Simplification, vertex cache/fetch/overdraw optimization | MIT |
| [V-HACD](https://github.com/kmammou/v-hacd) | Approximate convex decomposition | BSD-3 |
| [tinygltf](https://github.com/syoyo/tinygltf) | glTF/GLB loading and saving | MIT |
| [par_shapes](https://github.com/prideout/par) | Parametric primitive generation | MIT |

If a submodule is missing, its features are disabled at configure time and the library still builds.

## Usage

```cpp
#include "bromesh/primitives/primitives.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/simplify.h"
#include "bromesh/uv/projection.h"
#include "bromesh/optimization/optimize.h"
#include "bromesh/io/gltf.h"

// Generate a sphere, simplify it, add UVs, optimize, and save
auto mesh = bromesh::sphere(2.0f, 32, 24);
mesh = bromesh::simplify(mesh, 0.5f);
bromesh::computeNormals(mesh);
bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical);
bromesh::optimizeVertexCache(mesh);
bromesh::optimizeVertexFetch(mesh);
bromesh::saveGLTF(mesh, "sphere.glb");
```

### Isosurface extraction

```cpp
#include "bromesh/isosurface/marching_cubes.h"

// Build a scalar field (positive = inside)
const int N = 32;
std::vector<float> field(N * N * N);
// ... fill with SDF values ...

auto mesh = bromesh::marchingCubes(field.data(), N, N, N, 0.0f, 1.0f);
```

Surface nets, dual contouring, and transvoxel follow the same pattern. Transvoxel additionally accepts an LOD level and neighbor LOD array for seamless chunk boundaries.

### Convex decomposition

```cpp
#include "bromesh/analysis/convex_decomposition.h"

bromesh::ConvexDecompParams params;
params.maxHulls = 8;
auto hulls = bromesh::convexDecomposition(mesh, params);
// Each hull is a MeshData suitable for physics collision shapes
```

## Integration

bromesh is a static library. Add it as a CMake subdirectory:

```cmake
add_subdirectory(path/to/bromesh)
target_link_libraries(your_target PRIVATE bromesh)
```

All public headers are under `include/bromesh/`.

## License

MIT
