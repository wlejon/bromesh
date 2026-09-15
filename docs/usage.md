# Usage cookbook

Worked examples for each bromesh subsystem. Every algorithm produces or
consumes `bromesh::MeshData` (`include/bromesh/mesh_data.h`) — flat position /
normal / UV / color / tangent / index arrays ready for GPU upload or
TypedArray transfer.

- [Basics](#basics)
- [Isosurface extraction](#isosurface-extraction)
- [Subdivision](#subdivision)
- [Smoothing and remeshing](#smoothing-and-remeshing)
- [Sweep and extrusion](#sweep-and-extrusion)
- [Polygon triangulation](#polygon-triangulation)
- [Half-edge edit mesh (PolyMesh)](#half-edge-edit-mesh-polymesh)
- [Point cloud reconstruction](#point-cloud-reconstruction)
- [Procedural foliage](#procedural-foliage)
- [Transforms](#transforms)
- [Skinning and morph targets](#skinning-and-morph-targets)
- [Auto-rigging](#auto-rigging)
- [Skin weight transfer](#skin-weight-transfer)
- [Pose evaluation and animation](#pose-evaluation-and-animation)
- [Inverse kinematics](#inverse-kinematics)
- [Retargeting and procedural locomotion](#retargeting-and-procedural-locomotion)
- [Shrinkwrap](#shrinkwrap)
- [Vertex color baking](#vertex-color-baking)
- [Texture-space baking](#texture-space-baking)
- [High-poly to low-poly transfer](#high-poly-to-low-poly-transfer)
- [UV quality metrics](#uv-quality-metrics)
- [Surface sampling](#surface-sampling)
- [Progressive mesh (continuous LOD)](#progressive-mesh-continuous-lod)
- [Raycasting and queries](#raycasting-and-queries)
- [Self-intersection and mesh-mesh intersection](#self-intersection-and-mesh-mesh-intersection)
- [Boolean/CSG operations](#booleancsg-operations)
- [Convex decomposition](#convex-decomposition)
- [PLY and FBX I/O](#ply-and-fbx-io)
- [Gaussian splats](#gaussian-splats)

## Basics

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

## Isosurface extraction

```cpp
#include "bromesh/isosurface/marching_cubes.h"

// Build a scalar field. Standard SDF convention: f < iso is inside,
// f >= iso is outside. e.g. for a sphere of radius r: f = length(p) - r.
const int N = 32;
std::vector<float> field(N * N * N);
// ... fill with SDF values ...

auto mesh = bromesh::marchingCubes(field.data(), N, N, N, 0.0f, 1.0f);
```

`surfaceNets` and `dualContour` share this signature. `transvoxel` differs:
it takes a single cubic `gridSize` plus an LOD level and a six-entry neighbor
LOD array (before the iso level) for seamless chunk boundaries:

```cpp
#include "bromesh/isosurface/transvoxel.h"

int neighborLods[6] = {1, 1, 0, 0, 1, 1};
auto chunk = bromesh::transvoxel(field.data(), N, /*lod=*/0, neighborLods);
```

## Subdivision

```cpp
#include "bromesh/manipulation/subdivide.h"

auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
auto smooth = bromesh::subdivideLoop(mesh, 2);          // 2 iterations of Loop subdivision
auto cc     = bromesh::subdivideCatmullClark(mesh, 1);   // Catmull-Clark
auto simple = bromesh::subdivideMidpoint(mesh, 1);       // No smoothing, pure refinement
```

## Smoothing and remeshing

```cpp
#include "bromesh/manipulation/smooth.h"
#include "bromesh/manipulation/remesh.h"

bromesh::smoothTaubin(mesh, 0.5f, -0.53f, 10);  // Volume-preserving smoothing
auto remeshed = bromesh::remeshIsotropic(mesh);  // Uniform triangle quality
```

## Sweep and extrusion

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

## Polygon triangulation

```cpp
#include "bromesh/manipulation/polygon.h"

// 2D outline (flat x,y,...) with optional holes, triangulated via manifold
auto fill = bromesh::triangulatePolygon2D(outerXY, /*holes=*/{innerXY});

// Planar 3D polygon projected to its plane and indexed back to 3D
float n[3] = {0, 1, 0};
auto cap = bromesh::triangulatePolygon3D(outerXYZ, holesXYZ, n);
```

## Half-edge edit mesh (PolyMesh)

```cpp
#include "bromesh/manipulation/poly_mesh.h"

// Lift a triangle soup into N-gon faces, edit topology, tessellate back
auto pm = bromesh::PolyMesh::fromMeshData(mesh.positions, mesh.indices, triToGroup);
pm.mergeFacesByGroup();                 // coplanar tris -> single N-gon faces
float up[3] = {0, 0.5f, 0};
pm.extrudeFace(faceIdx, up);            // SketchUp-style push/pull
auto tess = pm.tessellate();            // render-ready triangles + triToFace map
```

## Point cloud reconstruction

```cpp
#include "bromesh/reconstruction/reconstruct.h"

bromesh::ReconstructParams params;
params.gridResolution = 64;
auto mesh = bromesh::reconstructFromPointCloud(positions, normals, pointCount, params);
```

## Procedural foliage

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

## Transforms

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

## Skinning and morph targets

```cpp
#include "bromesh/manipulation/skin.h"

bromesh::normalizeWeights(skin);                          // Clean up weights
bromesh::applySkinning(mesh, skin, poseMatrices);         // Skeletal animation
bromesh::applyMorphTarget(mesh, morphTarget, 0.5f);      // 50% blend shape
```

## Auto-rigging

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

## Skin weight transfer

Project skin weights from a source mesh onto a target (e.g. swappable armor
that should ride the same skeleton as a base body).

```cpp
#include "bromesh/manipulation/skin_transfer.h"

auto armorSkin = bromesh::transferSkinWeights(armorMesh, bodyMesh, bodySkin);
```

## Pose evaluation and animation

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

// Socket world matrix (returns std::nullopt if the socket name is unknown)
if (auto socket = bromesh::socketWorldMatrix(skeleton, pose, "hand.R")) {
    // (*socket)[0..15] is the column-major world matrix
}
```

## Inverse kinematics

```cpp
#include "bromesh/animation/ik.h"

float target[3] = {0.3f, 1.2f, 0.4f};
bromesh::solveTwoBoneIK(skeleton, pose, shoulder, elbow, wrist, target);
bromesh::solveFABRIK(skeleton, pose, spineChain, target);
bromesh::solveLookAt(skeleton, pose, headBone, target);
```

## Retargeting and procedural locomotion

```cpp
#include "bromesh/animation/retarget.h"
#include "bromesh/animation/locomotion.h"

auto retargeted = bromesh::retargetAnimation(srcAnim, srcSkeleton, dstSkeleton);

bromesh::LocomotionParams params;
params.cycleDuration = 1.0f;
params.strideLength  = 0.30f;
auto walk = bromesh::generateLocomotionCycle(skeleton, spec, params);
```

## Shrinkwrap

Project one mesh's vertices onto another (armor to body, cloth to form, etc.).

```cpp
#include "bromesh/manipulation/shrinkwrap.h"

bromesh::shrinkwrap(armor, body, bromesh::ShrinkwrapMode::Nearest,
                    /*maxDistance=*/0.0f, /*offset=*/0.002f);
```

## Vertex color baking

```cpp
#include "bromesh/analysis/bake.h"

bromesh::bakeAmbientOcclusion(mesh, 64);   // AO into vertex colors
bromesh::bakeCurvature(mesh, 1.0f);        // Mean curvature visualization
bromesh::bakeThickness(mesh, 32);          // Thickness for SSS approximation
```

## Texture-space baking

```cpp
#include "bromesh/analysis/bake_texture.h"

auto aoMap  = bromesh::bakeAmbientOcclusionToTexture(mesh, 512, 512);
auto curvMap = bromesh::bakeCurvatureToTexture(mesh, 512, 512);
auto nrmMap  = bromesh::bakeNormalsToTexture(mesh, 1024, 1024);
auto posMap  = bromesh::bakePositionToTexture(mesh, 1024, 1024);
// Access pixels: aoMap.at(x, y) returns float*, aoMap.pixels for raw data
```

## High-poly to low-poly transfer

```cpp
#include "bromesh/analysis/bake_transfer.h"

auto nrmMap = bromesh::bakeNormalsFromReference(lowPoly, highPoly, 1024, 1024);
auto aoMap  = bromesh::bakeAOFromReference(lowPoly, highPoly, 512, 512, /*rays=*/64);
```

## UV quality metrics

```cpp
#include "bromesh/uv/uv_metrics.h"

auto metrics = bromesh::measureUVQuality(mesh);
// metrics.avgStretch, metrics.maxStretch, metrics.avgAngleDistortion, metrics.uvSpaceUsage

auto perTri = bromesh::computeUVDistortion(mesh);  // Per-triangle stretch/area/angle
```

## Surface sampling

```cpp
#include "bromesh/analysis/sample.h"

auto points = bromesh::sampleSurface(mesh, 10000, 42);  // 10k uniform random points
float area = bromesh::computeSurfaceArea(mesh);
```

## Progressive mesh (continuous LOD)

```cpp
#include "bromesh/optimization/progressive.h"

auto pm = bromesh::buildProgressiveMesh(mesh);             // Build once
auto lod = bromesh::progressiveMeshAtRatio(pm, 0.5f);     // Extract at 50%
auto lod2 = bromesh::progressiveMeshAtTriangleCount(pm, 1000); // Or by triangle count

// Serialize for streaming (send coarse first, then refine)
auto data = bromesh::serializeProgressiveMesh(pm);
auto pm2 = bromesh::deserializeProgressiveMesh(data.data(), data.size());
```

## Raycasting and queries

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

## Self-intersection and mesh-mesh intersection

```cpp
#include "bromesh/analysis/intersect.h"

bool bad = bromesh::hasSelfIntersections(mesh);           // Fast boolean check
auto pairs = bromesh::findSelfIntersections(mesh);        // All intersecting triangle pairs
bool overlap = bromesh::meshesIntersect(meshA, meshB);    // Cross-mesh test
```

## Boolean/CSG operations

```cpp
#include "bromesh/csg/boolean.h"

auto result = bromesh::booleanDifference(cube, sphere);
auto [top, bottom] = bromesh::splitByPlane(mesh, 0, 1, 0, 0);
```

## Convex decomposition

```cpp
#include "bromesh/analysis/convex_decomposition.h"

bromesh::ConvexDecompParams params;
params.maxHulls = 8;
auto hulls = bromesh::convexDecomposition(mesh, params);
// Each hull is a MeshData suitable for physics collision shapes
```

## PLY and FBX I/O

```cpp
#include "bromesh/io/ply.h"
#include "bromesh/io/fbx.h"

auto mesh = bromesh::loadPLY("scan.ply");
bromesh::savePLY(mesh, "output.ply");

auto fbxMeshes = bromesh::loadFBX("model.fbx");  // Returns all meshes in the scene
```

## Gaussian splats

`GaussianSplatCloud` holds a 3D Gaussian Splat cloud in the same separate-stream
layout as `MeshData`, stored render-ready: linear scales, `[0,1]` opacities,
normalized `xyzw` rotations, and spherical-harmonic coefficients interleaved by
coefficient (`[r0 g0 b0 r1 g1 b1 ...]`, stride `shStride()`). The PLY loader and
saver convert to and from the log-scale / logit-opacity, channel-major storage
that standard 3DGS files (INRIA / PlayCanvas) use.

```cpp
#include "bromesh/io/splat_ply.h"
#include "bromesh/manipulation/splat_ops.h"

auto cloud = bromesh::loadSplatPLY("scene.ply");  // SH degree inferred from f_rest_*
// cloud.count(), cloud.positions, cloud.scales, cloud.rotations,
// cloud.opacities, cloud.sh, cloud.shDegree
bromesh::saveSplatPLY(cloud, "scene_out.ply");    // binary little-endian, round-trips
```

### Splat operations and mesh conversion

Apply rigid or affine transformations, filter/crop by opacity, scale, or spatial bounding boxes, merge heterogeneous SH-degree clouds, or sample a triangle mesh directly into surface-aligned Gaussian splat disks:

```cpp
#include "bromesh/manipulation/splat_ops.h"

// 1. Transform: translate, scale, or apply a column-major 4x4 affine matrix
bromesh::translateSplats(cloud, 0.0f, 1.5f, 0.0f);
bromesh::scaleSplats(cloud, 2.0f, 2.0f, 2.0f);
// float m[16] = { ... };
// bromesh::transformSplats(cloud, m);

// 2. Filter & crop: drop low opacity, outside AABB, or excessively large splats
bromath::AABB3 cropBox{{-5.0f, -5.0f, -5.0f}, {5.0f, 5.0f, 5.0f}};
bromesh::SplatFilterOptions filterOpts;
filterOpts.minOpacity = 0.05f;
filterOpts.cropBox = &cropBox;
filterOpts.maxScale = 1.0f;
bromesh::filterSplats(cloud, filterOpts);

// 3. Merge clouds: unifies SH degrees to the maximum (zero-padding lower degrees)
std::vector<bromesh::GaussianSplatCloud> clouds = {cloudA, cloudB};
auto unified = bromesh::mergeSplats(clouds);

// 4. Convert MeshData into render-ready Gaussian splats
bromesh::MeshToSplatsOptions splatOpts;
splatOpts.splatCount = 10000;
splatOpts.opacity = 0.95f;
auto splatCloud = bromesh::meshToSplats(mesh, splatOpts);
```

## Draco mesh compression

Draco compression dramatically reduces mesh storage and transmission size via lossy attribute quantization and edgebreaker connectivity compression. Use `encodeDraco` to compress a `MeshData` to memory (.drc bytes) and `decodeDraco` to decompress.

```cpp
#include "bromesh/io/draco.h"

// Compress a mesh with customizable quantization and speed knobs
bromesh::DracoEncodeOptions opts;
opts.positionBits = 14;
opts.normalBits = 10;
opts.uvBits = 12;
opts.colorBits = 8;
opts.speed = 7;

std::string error;
std::vector<uint8_t> compressed = bromesh::encodeDraco(mesh, opts, &error);
if (!error.empty()) {
    // Handle error
}

// Decompress from memory
bromesh::DracoDecoded decoded = bromesh::decodeDraco(compressed.data(), compressed.size());
if (decoded.ok()) {
    bromesh::MeshData restoredMesh = decoded.mesh;
    // decoded.attributes provides access to all raw decoded attributes if needed
}
```
