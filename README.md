# bromesh

[![CI](https://github.com/wlejon/bromesh/actions/workflows/ci.yml/badge.svg)](https://github.com/wlejon/bromesh/actions/workflows/ci.yml)
[![CodeQL](https://github.com/wlejon/bromesh/actions/workflows/codeql.yml/badge.svg)](https://github.com/wlejon/bromesh/actions/workflows/codeql.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A C++20 static library for mesh generation, manipulation, and I/O. Designed for game engines, tools, and procedural content pipelines.

bromesh is one of the engine libraries of the [bro ecosystem](https://github.com/wlejon/bro/blob/main/docs/ecosystem.md): [bro](https://github.com/wlejon/bro) links it for its 3D subsystem and exposes it to apps as `bro.mesh` and `bro.rigging` through the JavaScript binding in `src/api/` (`bromesh_api`), and [broflora](https://github.com/wlejon/broflora) builds its plant geometry with it. It depends on [bromath](https://github.com/wlejon/bromath); the binding needs [bronze](https://github.com/wlejon/bronze) and [brass](https://github.com/wlejon/brass), and brass also compiles signed-distance-field graphs to native code for the isosurface extractors.

Everything runs on the CPU. It is built and tested on Windows (MSVC), Linux (GCC and Clang) and macOS (arm64).

## Features

| Category | Algorithms |
|---|---|
| **Isosurface** | Marching cubes, surface nets, dual contouring (QEF), transvoxel (seamless LOD) |
| **Voxel** | Greedy meshing with palette colors |
| **Primitives** | Box, sphere, cylinder, capsule, plane, torus, heightmap grid; parametric: geodesic sphere, Platonic solids, cone, disc, rock, blob (rock + scale + translate) |
| **Sweep / extrusion** | Sweep a 2D profile along a 3D path (parallel-transport frames, per-ring scale/twist, mitered joints), circular `tube`, cubic-bezier sweep; triangulate 2D/planar-3D polygons (with holes, via manifold) |
| **Edit mesh** | `PolyMesh` half-edge adjacency over N-gon faces — extrude/translate face, split/flip/collapse edge, group merge to N-gons, tessellation, validation, compaction; the editable topology behind isotropic remeshing |
| **Procedural** | Foliage cards (leaf/petal, flower, blade strip), space-colonization branch trees + pipe-model thickening, L-system turtle geometry, leaf scatter, obstacle/capsule fields — the renderer surface the `broflora` sibling composes |
| **Subdivision** | Loop, Catmull-Clark, midpoint (iterative, arbitrary depth) |
| **Manipulation** | Smooth/flat normals, tangents, simplify (quadric error + attribute-aware), target triangle count decimation, LOD chain, weld, split components, merge, repair (degenerate/duplicate removal, hole filling), translate/rotate/scale/mirror/center/transform, shrinkwrap (nearest / normal-project / axis-project) |
| **Skinning** | Apply bone transforms (4 weights/vertex), morph target blending, weight normalization, closest-point skin weight transfer between meshes |
| **Rigging** | One-call auto-rig from landmarks + RigSpec (bundled humanoid/quadruped/hexapod/octopod specs, also loadable/serializable as JSON); geometric landmark detection; skeleton fitting; skin weighting via bone heat, bounded biharmonic weights (BBW, OSQP), or voxel-bind; Laplacian weight smoothing; skin validation |
| **Animation** | Pose evaluation (bind/animation/two-way and N-way weighted blend with bone masks), world & skinning matrix composition, socket resolution, two-bone IK, FABRIK, look-at IK, name-based animation retargeting (Rigify/Mixamo), procedural locomotion cycles (biped/quadruped/hexapod/octopod gaits) |
| **Smoothing** | Laplacian, Taubin (shrinkage-free) |
| **Remeshing** | Isotropic remeshing (edge split/collapse/flip + tangential relaxation) |
| **Reconstruction** | Point cloud to mesh via implicit surface estimation + marching cubes |
| **Analysis** | Bounding box, manifold check, volume, surface area, triangle areas, convex decomposition (V-HACD), convex hull, surface sampling, self-intersection detection |
| **Queries** | Raycast (closest/all/test), closest point on surface, mesh-mesh intersection test |
| **Baking** | Ambient occlusion, mean curvature, thickness — to vertex colors or UV-space textures; world-space normal maps, position maps; high-poly→low-poly transfer of tangent-space normals and AO |
| **UV** | Box, planar (XY/XZ/YZ), cylindrical, spherical projection; automatic unwrapping and atlas packing (xatlas); quality metrics (L2 stretch, area/angle distortion, packing efficiency) |
| **Optimization** | Vertex cache, vertex fetch, overdraw, meshlet generation, spatial sorting, shadow index buffer, mesh encoding/compression, triangle strips, progressive mesh (continuous LOD with serialization) |
| **Boolean/CSG** | Union, difference, intersection, plane splitting (manifold) |
| **Gaussian splats** | `GaussianSplatCloud` (SoA positions/scales/rotations/opacities/SH, degrees 0-3), 3DGS `.ply` read/write with activation on load and inverse on save, affine transform/translate/scale, opacity/crop/scale filtering, heterogeneous cloud merging, mesh surface to splats |
| **I/O** | OBJ read/write, STL read/write, PLY read/write, glTF/GLB read/write (meshes, skins, skeletons, animations, materials, embedded images), Draco decode/encode, FBX read, MagicaVoxel VOX read, 3DGS splat PLY read/write |

All mesh algorithms produce `bromesh::MeshData` -- a flat struct with separate position, normal, UV, color, tangent, and index arrays ready for GPU upload or TypedArray transfer. `GaussianSplatCloud` mirrors that layout for splat data.

## Building

```bash
git clone https://github.com/wlejon/bromesh.git
cd bromesh
cmake -B build
cmake --build build --config Release
```

Requires CMake 3.24+ and a C++20 compiler (MSVC 2022, GCC 12+, Clang 15+).
A plain clone is all it takes: every dependency is pinned to a commit in
`CMakeLists.txt` (`bro_dependency()`, `cmake/bro_deps.cmake`). A working tree
beside bromesh (`../bromath`, `../bronze`, ...) wins when present; otherwise the
head of its main branch is fetched at configure. `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>`
points one dependency anywhere else. bronze and brass compile inside bromesh's
build tree, because the JavaScript binding has to be compiled against the same
bronze as the program that loads it. `-DBROMESH_ENABLE_API=OFF` skips the
binding itself.

### Running tests

```bash
cd build
ctest --build-config Release --output-on-failure
```

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

Worked examples for every subsystem — isosurfaces, subdivision, sweeps,
PolyMesh editing, procedural foliage, auto-rigging, animation/IK, baking,
progressive LOD, CSG, splats, and more — live in the
[usage cookbook](docs/usage.md). The auto-rigging design is documented in
[docs/auto-rig-strategy.md](docs/auto-rig-strategy.md).

## Dependencies

Every dependency is either declared with `bro_dependency()` and resolved as
described under Building, or vendored under `third_party/`.

[bromath](https://github.com/wlejon/bromath) (header-only Vec/Quat/Mat, AABB,
curves, easing, `SpatialHash3D`) is the one required library. An existing
`bromath` target wins (one copy serves a whole build), then `../bromath`, then
the head of its main branch.

The rest are third-party. draco and par_shapes are vendored in tree; the
others are fetched at their pinned commits:

| Library | Purpose | License |
|---|---|---|
| [meshoptimizer](https://github.com/zeux/meshoptimizer) | Simplification, vertex cache/fetch/overdraw optimization, encoding, meshlets | MIT |
| [V-HACD](https://github.com/kmammou/v-hacd) | Approximate convex decomposition | BSD-3 |
| [tinygltf](https://github.com/syoyo/tinygltf) | glTF/GLB loading and saving | MIT |
| [par_shapes](https://github.com/prideout/par) | Parametric primitive generation | MIT |
| [xatlas](https://github.com/jpcy/xatlas) | Automatic UV unwrapping and atlas packing | MIT |
| [manifold](https://github.com/elalish/manifold) | Boolean/CSG operations and 2D/3D polygon triangulation | Apache-2.0 |
| [OpenFBX](https://github.com/nem0/OpenFBX) | FBX file loading | MIT |
| [OSQP](https://github.com/osqp/osqp) | Quadratic program solver backing bounded biharmonic weights (BBW) | Apache-2.0 |
| [draco](https://github.com/google/draco) | Mesh compression (decode/encode) | Apache-2.0 |

## Integration

bromesh is a static library. Add it as a CMake subdirectory (resolve bromath
first if your project also uses it, so one copy serves both):

```cmake
add_subdirectory(path/to/bromesh)
target_link_libraries(your_target PRIVATE bromesh)
```

All public headers are under `include/bromesh/`. The JavaScript surface
(`bromesh_api`) is documented in bro's
[mesh-api.js](https://github.com/wlejon/bro/blob/main/docs/mesh-api.js),
[mesh-io-api.js](https://github.com/wlejon/bro/blob/main/docs/mesh-io-api.js),
[mesh-plants-api.js](https://github.com/wlejon/bro/blob/main/docs/mesh-plants-api.js) and
[rigging-api.js](https://github.com/wlejon/bro/blob/main/docs/rigging-api.js).

## Tests and CI

`ctest` runs the C++ suite plus the JavaScript binding tests. CI builds and
tests a plain clone on Linux (GCC and Clang), Windows (MSVC) and macOS/arm64
against the pinned dependencies, and reports coverage of `include/bromesh/` and `src/` in each
run's summary. [CodeQL](.github/workflows/codeql.yml) analyses the library on
every push and weekly.

## License

[MIT](LICENSE)
