# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`bromesh` is a C++20 static library for mesh generation, manipulation, I/O, rigging, skeletal animation, and Gaussian-splat data. It is consumed by other projects as a CMake subdirectory (`target_link_libraries(... bromesh)`); there is no public app, only the library + a single test binary.

## Build and test

```bash
# The bromath sibling must sit at ../bromath (header-only; configure FATAL_ERRORs
# without it — overridable via -DBROMATH_DIR=...). Submodules under third_party/
# must be present — features linked against missing submodules silently disable
# (see "Optional-by-submodule" below).
git submodule update --init --recursive

cmake -B build
cmake --build build --config Release

# Tests
cd build && ctest --build-config Release --output-on-failure
# or run the binary directly for full output
./build/tests/Release/bromesh_test          # Windows/MSVC
./build/tests/bromesh_test                  # single-config generators
```

CMake 3.24+, C++20 (MSVC 2022, GCC 12+, Clang 15+). MSVC uses static CRT (`/MT[d]`) **only** when bromesh is the top-level project — as a subdirectory it inherits the parent's runtime.

### Running a single test

The test harness is custom (`tests/test_framework.h`) — not GoogleTest. Tests self-register through the `TEST(name) { ... }` macro at static-init time, and `main()` in `test_main.cpp` walks the registry afterwards (deferring execution out of static init avoids cross-TU init-order failures on Windows). Assertions use `ASSERT(cond, msg)`, which counts rather than aborts, so a run reports `passed/run` and exits non-zero if they differ.

There is no filter flag and no `RUN_TEST` list to edit. To run a single test, temporarily filter inside `main()`'s registry loop on `t.name`, or comment out the `TEST` bodies you don't want.

Test sources: `test_mesh.cpp`, `test_animation.cpp`, `test_rigging.cpp`, `test_procedural.cpp`, `test_splat.cpp`, `test_coverage.cpp`, `test_coverage2.cpp` (shared rig fixtures in `synthetic_rigs.h`).

`ctest` runs the whole `bromesh_test` executable as one test case; `--output-on-failure` surfaces the per-assertion failure prints.

`scripts/coverage.ps1` runs the same suite under OpenCppCoverage (Windows/MSVC only) and writes an HTML report to `build/coverage/index.html`. It requires a Release build with PDBs (`-DCMAKE_MSVC_DEBUG_INFORMATION_FORMAT=ProgramDatabase`) — Debug builds hang on meshoptimizer's modal assert dialog.

## Architecture

### The core contract: `MeshData`

Every algorithm in the library produces or consumes `bromesh::MeshData` (`include/bromesh/mesh_data.h`) — a flat struct of separate `std::vector<float>` streams (positions, normals, uvs, colors, tangents) plus a `std::vector<uint32_t>` index buffer. The design goal is **direct GPU upload and zero-copy TypedArray transfer to JS**; do not introduce AoS vertex structs or interleaved buffers. Companion structs `SkinData`, `Skeleton`, `Bone`, `Socket`, `Animation`, `AnimChannel`, `MorphTarget` live in the same header and are the canonical types across the animation/rigging stack. The AABB type comes from `bromath` (`bromath::AABB3`) — bromesh does not define its own.

`GaussianSplatCloud` (`include/bromesh/gaussian_splat.h`) is the one peer type that is not a mesh: the same SoA-of-`std::vector<float>` discipline applied to 3D Gaussian Splats (positions, scales, rotations, opacities, SH coefficients). It is stored **render-ready** — linear scales, `[0,1]` opacity, normalized `xyzw` quaternions, SH interleaved by coefficient — and `io/splat_ply.h` converts to/from the log-scale / logit-opacity, channel-major encoding that on-disk 3DGS `.ply` files use. Keep the activation and the layout transpose confined to the loader/saver; everything in memory is already activated.

Public headers live under `include/bromesh/<module>/...` and mirror `src/<module>/...`. The only headers at the root are the shared data contracts, `mesh_data.h` and `gaussian_splat.h`, which are pure struct definitions with no `.cpp`. Each subsystem is a self-contained module with narrow entry points — they are composed by the caller, not by the library.

### Module layout (big picture)

```
isosurface/     marching cubes, surface nets, dual contouring (QEF), transvoxel
voxel/          greedy meshing, voxel chunks
primitives/     analytic (box/sphere/...) + par_shapes (platonic, knot, rock, ...)
manipulation/   normals, simplify, subdivide, smooth, remesh, repair, weld, merge,
                split_components, transform, skin, skin_transfer, shrinkwrap,
                sweep + bezier_sweep (profile-along-path extrusion / tube),
                polygon (triangulate 2D/planar-3D outlines via manifold),
                poly_mesh (PolyMesh half-edge N-gon edit mesh)
analysis/       bbox (+ volume), manifold_check, bvh, raycast, intersect,
                sample, bake (vertex + texture + hi→lo transfer),
                convex_decomposition (+ convexHull)
animation/      pose, ik (two-bone / FABRIK / look-at), retarget, locomotion
rigging/        auto_rig driver, landmarks + landmark_detect, rig_spec, skeleton_fit,
                weighting (bone_heat, bbw, voxel_bind), weight_smooth, skin_validate,
                mesh_laplacian
uv/             projection, unwrap (xatlas), uv_metrics
optimization/   optimize (meshopt), meshlets, analyze, spatial, encode, strips,
                progressive (continuous LOD with serialization)
reconstruction/ point cloud → mesh
csg/            boolean (manifold)
procedural/     stateless plant-shaped mesh constructors. plants.h:
                leafCard / flower (4x4-atlas foliage + petal cards, LeafShape,
                cup/silhouette controls), bladeStrip + bladePath (swept blades),
                tree (attractor cloud -> colonize -> pipe-model thicken -> tube).
                space_colonization + branches (attractor-grown branch tubes),
                lsystem + lsystem_turtle (rule-based geometry), leaf_scatter
                (branch-driven, per-segment density weighting),
                obstacle_field (CapsuleField for obstacle-aware growth/scatter)
io/             obj, stl, ply, gltf, fbx (read-only), vox (read-only),
                splat_ply (3DGS GaussianSplatCloud read/write)
```

The `procedural/` module is the renderer surface that the sibling
`broflora` ecosystem-simulation library composes — broflora ticks plant
state forward, then calls into these constructors to emit `MeshData`.
Standalone consumers (game backgrounds, UI icons, distant LODs) use the
same calls directly without dragging in broflora.

### The rigging + animation pipeline (multi-file flow)

This is the one area where the separation of modules matters because it spans many files. The end-to-end authoring flow is:

1. **Landmarks** (`rigging/landmark_detect.h` *or* user-supplied) → `Landmarks` dict.
2. **`rigging/rig_spec.h`** declares the skeleton template + landmark schema. Bundled specs: `builtinHumanoidSpec()`, `builtinQuadrupedSpec()`, `builtinHexapodSpec()`, `builtinOctopodSpec()`, or by name via `builtinRigSpec()`. Specs round-trip as JSON (`parseRigSpecJSON` / `serializeRigSpecJSON` / `loadRigSpecFile`); the bundled ones are mirrored as data under `data/rig_specs/`. Bone positions are expressions (`"landmark:wrist_L"`, `"mid:A,B"`, `"lerp:A,B,t"`, `"offset:A,dx,dy,dz"`) resolved against the landmark dict.
3. **`rigging/auto_rig.h`** — `autoRig(mesh, spec, landmarks, wopts)` is the one-call driver. It runs `skeleton_fit` → weighting (dispatched by `WeightingOptions::method`) → `weight_smooth`.
4. **Weighting methods** (`rigging/weighting.h` dispatcher):
   - `BoneHeat` (`bone_heat.cpp`) — default for manifold meshes.
   - `BBW` (`bbw.cpp`) — opt-in, uses OSQP QP solver per bone.
   - `VoxelBind` (`voxel_bind.cpp`) — fallback for non-manifold / messy topology.
   - `Auto` — `autoSelectWeightingMethod(mesh)` picks BoneHeat vs VoxelBind based on manifoldness.
5. **Runtime animation** (`animation/pose.h`): `bindPose` → `evaluateAnimation` → `blendPoses` (two-way, in place) or `blendPosesN` (weighted N-way into an out-pose, for blend spaces) → `computeSkinningMatrices` → `manipulation/skin.h::applySkinning`. IK (`ik.h`) mutates `Pose` in place between evaluation and skinning.

The strategy doc at `docs/auto-rig-strategy.md` is the north star for this subsystem.

### Optional-by-submodule pattern

The **one hard dependency** is the `bromath` sibling (`../bromath`, header-only,
linked as `bromath::bromath`) — it backs `MeshData`'s AABB plus the Vec/Quat/Mat
and `SpatialHash3D` types used throughout sweep/procedural/analysis. Everything
else is optional.

Every third-party dependency (meshoptimizer, V-HACD, tinygltf, par_shapes, xatlas, manifold, OpenFBX, OSQP) is a git submodule under `third_party/`. The top-level `CMakeLists.txt` gates each by an `EXISTS` check and sets `BROMESH_HAS_<DEP>` accordingly, forwarded as a public compile definition. Note `manifold` backs both Boolean/CSG (`csg/boolean.cpp`) and polygon triangulation (`manipulation/polygon.cpp`); the latter no-ops to an empty `MeshData` when `BROMESH_HAS_MANIFOLD==0`.

**Implication for new code**: any `.cpp` that uses an optional dep must compile to a working no-op (usually returning an empty `MeshData` or `false`) when its `BROMESH_HAS_...` macro is undefined. Do not add a hard dependency on any submodule. The library must build and link with every submodule absent.

tinygltf is compiled with `TINYGLTF_NO_EXTERNAL_IMAGE` / `TINYGLTF_NO_STB_IMAGE_WRITE` — embedded images are decoded into `GltfScene::images` but external-file images are intentionally not loaded (sandbox reasons).

### Conventions specific to this codebase

- **No AoS vertex structs.** Flat parallel arrays, always.
- **Column-major 4×4 matrices**, stored as `float[16]`. Quaternions are `xyzw`. Poses use a 10-float-per-bone stride (T3 + R4 + S3).
- **Bones are topologically ordered** (parents precede children). Code walking the hierarchy relies on this — preserve it when constructing or editing a `Skeleton`.
- **Skeleton has a `rootTransform`** applied to root bones (parent == -1) at world-matrix composition time. Used to keep scene-tree ancestry (e.g. a glTF armature scale) out of the bones' local TRS so animation keyframes stay in the original space.
- **Rigging/animation APIs are consumed from language bindings (JS/TypedArray).** When touching the public surface, preserve stable string names (see `weightingMethodName` / `parseWeightingMethod`) and avoid adding types that don't flatten to POD arrays.

## Dependency policy

Prefer git submodules under `third_party/` for new deps — do not vendor copies of third-party source.
