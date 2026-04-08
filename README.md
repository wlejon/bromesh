# bromesh

A C++20 static library for mesh generation, manipulation, and I/O. Designed for game engines, tools, and procedural content pipelines.

## Features

| Category | Algorithms |
|---|---|
| **Isosurface** | Marching cubes, surface nets, dual contouring (QEF), transvoxel (seamless LOD) |
| **Voxel** | Greedy meshing with palette colors |
| **Primitives** | Box, sphere, cylinder, capsule, plane, torus, heightmap grid; parametric: geodesic sphere, Platonic solids, cone, disc, rock, trefoil knot, Klein bottle |
| **Subdivision** | Loop, Catmull-Clark, midpoint (iterative, arbitrary depth) |
| **Manipulation** | Smooth/flat normals, tangents, simplify (quadric error + attribute-aware), target triangle count decimation, LOD chain, weld, split components, merge, repair (degenerate/duplicate removal, hole filling) |
| **Smoothing** | Laplacian, Taubin (shrinkage-free) |
| **Remeshing** | Isotropic remeshing (edge split/collapse/flip + tangential relaxation) |
| **Reconstruction** | Point cloud to mesh via implicit surface estimation + marching cubes |
| **Analysis** | Bounding box, manifold check, volume, convex decomposition (V-HACD), convex hull |
| **Baking** | Ambient occlusion, mean curvature, thickness — all to vertex colors |
| **UV** | Box, planar (XY/XZ/YZ), cylindrical, spherical projection; automatic unwrapping and atlas packing (xatlas) |
| **Optimization** | Vertex cache, vertex fetch, overdraw, meshlet generation, spatial sorting, shadow index buffer, mesh encoding/compression, triangle strips |
| **Boolean/CSG** | Union, difference, intersection, plane splitting (manifold) |
| **I/O** | OBJ read/write, STL read/write, PLY read/write, glTF/GLB read/write, FBX read, MagicaVoxel VOX read |

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
| [meshoptimizer](https://github.com/zeux/meshoptimizer) | Simplification, vertex cache/fetch/overdraw optimization, encoding, meshlets | MIT |
| [V-HACD](https://github.com/kmammou/v-hacd) | Approximate convex decomposition | BSD-3 |
| [tinygltf](https://github.com/syoyo/tinygltf) | glTF/GLB loading and saving | MIT |
| [par_shapes](https://github.com/prideout/par) | Parametric primitive generation | MIT |
| [xatlas](https://github.com/jpcy/xatlas) | Automatic UV unwrapping and atlas packing | MIT |
| [manifold](https://github.com/elalish/manifold) | Boolean/CSG operations | Apache-2.0 |
| [OpenFBX](https://github.com/nem0/OpenFBX) | FBX file loading | MIT |

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

### Subdivision

```cpp
#include "bromesh/manipulation/subdivide.h"

auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
auto smooth = bromesh::subdivideLoop(mesh, 2);          // 2 iterations of Loop subdivision
auto cc     = bromesh::subdivideCatmullClark(mesh, 1);   // Catmull-Clark
auto simple = bromesh::subdivideMidpoint(mesh, 1);       // No smoothing, pure refinement
```

### Smoothing and remeshing

```cpp
#include "bromesh/manipulation/smooth.h"
#include "bromesh/manipulation/remesh.h"

bromesh::smoothTaubin(mesh, 0.5f, -0.53f, 10);  // Volume-preserving smoothing
auto remeshed = bromesh::remeshIsotropic(mesh);  // Uniform triangle quality
```

### Point cloud reconstruction

```cpp
#include "bromesh/reconstruction/reconstruct.h"

bromesh::ReconstructParams params;
params.gridResolution = 64;
auto mesh = bromesh::reconstructFromPointCloud(positions, normals, pointCount, params);
```

### Vertex color baking

```cpp
#include "bromesh/analysis/bake.h"

bromesh::bakeAmbientOcclusion(mesh, 64);   // AO into vertex colors
bromesh::bakeCurvature(mesh, 1.0f);        // Mean curvature visualization
bromesh::bakeThickness(mesh, 32);          // Thickness for SSS approximation
```

### Boolean/CSG operations

```cpp
#include "bromesh/csg/boolean.h"

auto result = bromesh::booleanDifference(cube, sphere);
auto [top, bottom] = bromesh::splitByPlane(mesh, 0, 1, 0, 0);
```

### Convex decomposition

```cpp
#include "bromesh/analysis/convex_decomposition.h"

bromesh::ConvexDecompParams params;
params.maxHulls = 8;
auto hulls = bromesh::convexDecomposition(mesh, params);
// Each hull is a MeshData suitable for physics collision shapes
```

### PLY and FBX I/O

```cpp
#include "bromesh/io/ply.h"
#include "bromesh/io/fbx.h"

auto mesh = bromesh::loadPLY("scan.ply");
bromesh::savePLY(mesh, "output.ply");

auto fbxMeshes = bromesh::loadFBX("model.fbx");  // Returns all meshes in the scene
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
