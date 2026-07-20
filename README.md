# bromesh

[![CI](https://github.com/wlejon/bromesh/actions/workflows/ci.yml/badge.svg)](https://github.com/wlejon/bromesh/actions/workflows/ci.yml)
[![CodeQL](https://github.com/wlejon/bromesh/actions/workflows/codeql.yml/badge.svg)](https://github.com/wlejon/bromesh/actions/workflows/codeql.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A C++20 static library for mesh generation, manipulation, and I/O. Designed for game engines, tools, and procedural content pipelines.

## Features

| Category | Algorithms |
|---|---|
| **Isosurface** | Marching cubes, surface nets, dual contouring (QEF), transvoxel (seamless LOD) |
| **Voxel** | Greedy meshing with palette colors |
| **Primitives** | Box, sphere, cylinder, capsule, plane, torus, heightmap grid; parametric: geodesic sphere, Platonic solids, cone, disc, rock, blob (rock + scale + translate), trefoil knot, Klein bottle |
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
| **Gaussian splats** | `GaussianSplatCloud` (SoA positions/scales/rotations/opacities/SH, degrees 0-3), 3DGS `.ply` read/write with activation on load and inverse on save |
| **I/O** | OBJ read/write, STL read/write, PLY read/write, glTF/GLB read/write (meshes, skins, skeletons, animations, materials, embedded images), FBX read, MagicaVoxel VOX read, 3DGS splat PLY read/write |

All mesh algorithms produce `bromesh::MeshData` -- a flat struct with separate position, normal, UV, color, tangent, and index arrays ready for GPU upload or TypedArray transfer. `GaussianSplatCloud` mirrors that layout for splat data.

## Building

```bash
git clone --recurse-submodules https://github.com/wlejon/bromesh.git
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

## Dependencies

bromesh links one sibling library, [bromath](https://github.com/wlejon/bromath)
(header-only Vec/Quat/Mat, AABB, curves, easing, `SpatialHash3D`), found at
`../bromath` relative to the repo. It is a hard dependency — configure fails if
it is missing.

Everything else is vendored under `third_party/` as git submodules or single
headers:

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

// Build a scalar field. Standard SDF convention: f < iso is inside,
// f >= iso is outside. e.g. for a sphere of radius r: f = length(p) - r.
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

### Sweep and extrusion

```cpp
#include "bromesh/manipulation/sweep.h"
#include "bromesh/manipulation/bezier_sweep.h"

// Sweep a square profile along a path with rotation-minimizing frames
std::vector<bromath::Vec2> profile = {{-1,-1},{1,-1},{1,1},{-1,1}};
std::vector<bromath::Vec3> path    = {{0,0,0},{0,2,0},{1,4,0}};
auto extruded = bromesh::sweep(profile, path);

// Circular tube (branches, vines, stems): unit circle scaled per ring
auto branch = bromesh::tube(path, /*radii=*/{0.3f, 0.2f, 0.05f});

// Sweep along a cubic-bezier polyline (4 control points per shared segment)
auto curved = bromesh::bezierSweep(controlPoints, profile);
```

### Polygon triangulation

```cpp
#include "bromesh/manipulation/polygon.h"

// 2D outline (flat x,y,...) with optional holes, triangulated via manifold
auto fill = bromesh::triangulatePolygon2D(outerXY, /*holes=*/{innerXY});

// Planar 3D polygon projected to its plane and indexed back to 3D
float n[3] = {0, 1, 0};
auto cap = bromesh::triangulatePolygon3D(outerXYZ, holesXYZ, n);
```

### Half-edge edit mesh (`PolyMesh`)

```cpp
#include "bromesh/manipulation/poly_mesh.h"

// Lift a triangle soup into N-gon faces, edit topology, tessellate back
auto pm = bromesh::PolyMesh::fromMeshData(mesh.positions, mesh.indices, triToGroup);
pm.mergeFacesByGroup();                 // coplanar tris -> single N-gon faces
float up[3] = {0, 0.5f, 0};
pm.extrudeFace(faceIdx, up);            // SketchUp-style push/pull
auto tess = pm.tessellate();            // render-ready triangles + triToFace map
```

### Point cloud reconstruction

```cpp
#include "bromesh/reconstruction/reconstruct.h"

bromesh::ReconstructParams params;
params.gridResolution = 64;
auto mesh = bromesh::reconstructFromPointCloud(positions, normals, pointCount, params);
```

### Procedural foliage

Stateless plant-shaped constructors that emit `MeshData` directly. The sibling
`broflora` ecosystem-sim library composes these as its renderer surface;
standalone callers (game backgrounds, distant LODs) use them the same way.

```cpp
#include "bromesh/procedural/plants.h"

// One-call tree: spherical attractor cloud -> space colonization -> pipe-model
// thickening -> merged branch tubes. Foliage is left to the caller.
bromesh::TreeOptions topts;
topts.canopyRadius = 3.0f;
auto t = bromesh::tree(topts);          // t.branches (MeshData), t.segments

// Low-poly leaf / petal cards with 4x4-atlas UVs, plus a radial flower
auto leaf   = bromesh::leafCard(bromesh::LeafShape::Oval);
auto bloom  = bromesh::flower(bromesh::FlowerOptions{});
```

### Transforms

```cpp
#include "bromesh/manipulation/transform.h"

bromesh::translateMesh(mesh, 0, 5, 0);          // Move up
bromesh::scaleMesh(mesh, 2.0f);                  // Uniform scale
bromesh::scaleMesh(mesh, 1, 2, 1);               // Non-uniform scale (normals auto-corrected)
bromesh::rotateMesh(mesh, 0, 1, 0, 3.14159f/4);  // 45° around Y axis
bromesh::mirrorMesh(mesh, 0);                     // Mirror across YZ plane (winding corrected)
bromesh::centerMesh(mesh);                        // Center at origin

float mat[16] = { /* column-major 4x4 */ };
bromesh::transformMesh(mesh, mat);                // Arbitrary affine transform
```

### Skinning and morph targets

```cpp
#include "bromesh/manipulation/skin.h"

bromesh::normalizeWeights(skin);                          // Clean up weights
bromesh::applySkinning(mesh, skin, poseMatrices);         // Skeletal animation
bromesh::applyMorphTarget(mesh, morphTarget, 0.5f);      // 50% blend shape
```

### Auto-rigging

One call fits a bundled `RigSpec` (humanoid, quadruped, hexapod, octopod) to a
mesh from a set of landmarks and produces a ready-to-skin skeleton + weights.
`WeightingMethod::Auto` picks bone-heat for manifold meshes and voxel-bind for
non-manifold input; BBW is opt-in (higher quality, costs a QP solve per bone).

```cpp
#include "bromesh/rigging/auto_rig.h"
#include "bromesh/rigging/landmark_detect.h"
#include "bromesh/rigging/rig_spec.h"

auto spec      = bromesh::builtinHumanoidSpec();
auto landmarks = bromesh::detectHumanoidLandmarks(mesh);  // or author by hand

bromesh::WeightingOptions wopts;
wopts.method = bromesh::WeightingMethod::Auto;             // or BoneHeat / BBW / VoxelBind
auto result  = bromesh::autoRig(mesh, spec, landmarks, wopts);
// result.skeleton, result.skin, result.missingLandmarks, result.warnings
```

Specs are also data. The bundled ones ship as JSON under `data/rig_specs/`
(`humanoid.json`, `quadruped.json`, `hexapod.json`, `octopod.json`) and a custom
skeleton template can be authored the same way:

```cpp
auto custom = bromesh::loadRigSpecFile("data/rig_specs/quadruped.json");
auto parsed = bromesh::parseRigSpecJSON(jsonText);
auto text   = bromesh::serializeRigSpecJSON(custom);
auto named  = bromesh::builtinRigSpec("hexapod");
```

### Skin weight transfer

Project skin weights from a source mesh onto a target (e.g. swappable armor
that should ride the same skeleton as a base body).

```cpp
#include "bromesh/manipulation/skin_transfer.h"

auto armorSkin = bromesh::transferSkinWeights(armorMesh, bodyMesh, bodySkin);
```

### Pose evaluation and animation

```cpp
#include "bromesh/animation/pose.h"

auto pose = bromesh::evaluateAnimation(skeleton, anim, tSeconds, /*loop=*/true);
bromesh::blendPoses(pose, upperBodyPose, 0.5f, boneMaskOrNull);

// N-way weighted blend (blend spaces: idle/walk/run, directional strafe sets).
// Weights are normalized internally; two sources reproduce blendPoses exactly.
const bromesh::Pose* sources[3] = {&idle, &walk, &run};
float weights[3] = {0.2f, 0.5f, 0.3f};
bromesh::Pose blended;
bromesh::blendPosesN(sources, weights, 3, blended, boneMaskOrNull);

std::vector<float> skinningMatrices;
bromesh::computeSkinningMatrices(skeleton, pose, skinningMatrices);
bromesh::applySkinning(mesh, skin, skinningMatrices.data());

float socket[16];
bromesh::socketWorldMatrix(skeleton, pose, "hand.R", socket);
```

### Inverse kinematics

```cpp
#include "bromesh/animation/ik.h"

float target[3] = {0.3f, 1.2f, 0.4f};
bromesh::solveTwoBoneIK(skeleton, pose, shoulder, elbow, wrist, target);
bromesh::solveFABRIK(skeleton, pose, spineChain, target);
bromesh::solveLookAt(skeleton, pose, headBone, target);
```

### Retargeting and procedural locomotion

```cpp
#include "bromesh/animation/retarget.h"
#include "bromesh/animation/locomotion.h"

auto retargeted = bromesh::retargetAnimation(srcAnim, srcSkeleton, dstSkeleton);

bromesh::LocomotionParams params;
params.cycleDuration = 1.0f;
params.strideLength  = 0.30f;
auto walk = bromesh::generateLocomotionCycle(skeleton, spec, params);
```

### Shrinkwrap

Project one mesh's vertices onto another (armor to body, cloth to form, etc.).

```cpp
#include "bromesh/manipulation/shrinkwrap.h"

bromesh::shrinkwrap(armor, body, bromesh::ShrinkwrapMode::Nearest,
                    /*maxDistance=*/0.0f, /*offset=*/0.002f);
```

### Vertex color baking

```cpp
#include "bromesh/analysis/bake.h"

bromesh::bakeAmbientOcclusion(mesh, 64);   // AO into vertex colors
bromesh::bakeCurvature(mesh, 1.0f);        // Mean curvature visualization
bromesh::bakeThickness(mesh, 32);          // Thickness for SSS approximation
```

### Texture-space baking

```cpp
#include "bromesh/analysis/bake_texture.h"

auto aoMap  = bromesh::bakeAmbientOcclusionToTexture(mesh, 512, 512);
auto curvMap = bromesh::bakeCurvatureToTexture(mesh, 512, 512);
auto nrmMap  = bromesh::bakeNormalsToTexture(mesh, 1024, 1024);
auto posMap  = bromesh::bakePositionToTexture(mesh, 1024, 1024);
// Access pixels: aoMap.at(x, y) returns float*, aoMap.pixels for raw data
```

### High-poly → low-poly transfer

```cpp
#include "bromesh/analysis/bake_transfer.h"

auto nrmMap = bromesh::bakeNormalsFromReference(lowPoly, highPoly, 1024, 1024);
auto aoMap  = bromesh::bakeAOFromReference(lowPoly, highPoly, 512, 512, /*rays=*/64);
```

### UV quality metrics

```cpp
#include "bromesh/uv/uv_metrics.h"

auto metrics = bromesh::measureUVQuality(mesh);
// metrics.avgStretch, metrics.maxStretch, metrics.avgAngleDistortion, metrics.uvSpaceUsage

auto perTri = bromesh::computeUVDistortion(mesh);  // Per-triangle stretch/area/angle
```

### Surface sampling

```cpp
#include "bromesh/analysis/sample.h"

auto points = bromesh::sampleSurface(mesh, 10000, 42);  // 10k uniform random points
float area = bromesh::computeSurfaceArea(mesh);
```

### Progressive mesh (continuous LOD)

```cpp
#include "bromesh/optimization/progressive.h"

auto pm = bromesh::buildProgressiveMesh(mesh);             // Build once
auto lod = bromesh::progressiveMeshAtRatio(pm, 0.5f);     // Extract at 50%
auto lod2 = bromesh::progressiveMeshAtTriangleCount(pm, 1000); // Or by triangle count

// Serialize for streaming (send coarse first, then refine)
auto data = bromesh::serializeProgressiveMesh(pm);
auto pm2 = bromesh::deserializeProgressiveMesh(data.data(), data.size());
```

### Raycasting and queries

```cpp
#include "bromesh/analysis/raycast.h"

float origin[3] = {0, 5, 0};
float dir[3] = {0, -1, 0};
auto hit = bromesh::raycast(mesh, origin, dir);          // Closest hit
if (hit.hit) { /* hit.position, hit.normal, hit.distance, hit.triangleIndex */ }

auto all = bromesh::raycastAll(mesh, origin, dir);       // All hits, sorted by distance
bool any  = bromesh::raycastTest(mesh, origin, dir);     // Fast boolean test

float point[3] = {2, 3, 1};
auto cp = bromesh::closestPoint(mesh, point);            // Nearest surface point
```

### Self-intersection and mesh-mesh intersection

```cpp
#include "bromesh/analysis/intersect.h"

bool bad = bromesh::hasSelfIntersections(mesh);           // Fast boolean check
auto pairs = bromesh::findSelfIntersections(mesh);        // All intersecting triangle pairs
bool overlap = bromesh::meshesIntersect(meshA, meshB);    // Cross-mesh test
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

### Gaussian splats

`GaussianSplatCloud` holds a 3D Gaussian Splat cloud in the same separate-stream
layout as `MeshData`, stored render-ready: linear scales, `[0,1]` opacities,
normalized `xyzw` rotations, and spherical-harmonic coefficients interleaved by
coefficient (`[r0 g0 b0 r1 g1 b1 ...]`, stride `shStride()`). The PLY loader and
saver convert to and from the log-scale / logit-opacity, channel-major storage
that standard 3DGS files (INRIA / PlayCanvas) use.

```cpp
#include "bromesh/io/splat_ply.h"

auto cloud = bromesh::loadSplatPLY("scene.ply");  // SH degree inferred from f_rest_*
// cloud.count(), cloud.positions, cloud.scales, cloud.rotations,
// cloud.opacities, cloud.sh, cloud.shDegree
bromesh::saveSplatPLY(cloud, "scene_out.ply");    // binary little-endian, round-trips
```

## Integration

bromesh is a static library. Add it as a CMake subdirectory:

```cmake
add_subdirectory(path/to/bromesh)
target_link_libraries(your_target PRIVATE bromesh)
```

All public headers are under `include/bromesh/`.

## License

[MIT](LICENSE)
