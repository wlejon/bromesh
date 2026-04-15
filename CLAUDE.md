# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`bromesh` is a C++20 static library for mesh generation, manipulation, I/O, rigging, and skeletal animation. It is consumed by other projects as a CMake subdirectory (`target_link_libraries(... bromesh)`); there is no public app, only the library + a single test binary.

## Build and test

```bash
# Submodules must be present — features linked against missing submodules silently disable (see "Optional-by-submodule" below).
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

The test harness is custom (`tests/test_framework.h`) — not GoogleTest. Tests register via `TEST(name) { ... }` macros and run from `main()` in `test_main.cpp`. There is no built-in filter flag. To run one test, either:
- Comment out `RUN_TEST(...)` entries, or
- Temporarily edit `test_main.cpp` to call a specific registered function.

`ctest` runs the whole `bromesh_test` executable as one test case; `--output-on-failure` surfaces the per-TEST pass/fail prints.

## Architecture

### The core contract: `MeshData`

Every algorithm in the library produces or consumes `bromesh::MeshData` (`include/bromesh/mesh_data.h`) — a flat struct of separate `std::vector<float>` streams (positions, normals, uvs, colors) plus a `std::vector<uint32_t>` index buffer. The design goal is **direct GPU upload and zero-copy TypedArray transfer to JS**; do not introduce AoS vertex structs or interleaved buffers. Companion structs `SkinData`, `Skeleton`, `Bone`, `Socket`, `Animation`, `AnimChannel`, `MorphTarget` live in the same header and are the canonical types across the animation/rigging stack.

Public headers live under `include/bromesh/<module>/...` and mirror `src/<module>/...`. Each subsystem is a self-contained module with narrow entry points — they are composed by the caller, not by the library.

### Module layout (big picture)

```
isosurface/     marching cubes, surface nets, dual contouring (QEF), transvoxel
voxel/          greedy meshing, voxel chunks
primitives/     analytic (box/sphere/...) + par_shapes (platonic, knot, rock, ...)
manipulation/   normals, simplify, subdivide, smooth, remesh, repair, weld, merge,
                split_components, transform, skin, skin_transfer, shrinkwrap
analysis/       bbox, bvh, raycast, intersect, sample, bake (vertex + texture +
                hi→lo transfer), convex_decomposition
animation/      pose, ik (two-bone / FABRIK / look-at), retarget, locomotion
rigging/        auto_rig driver, landmarks + landmark_detect, rig_spec, skeleton_fit,
                weighting (bone_heat, bbw, voxel_bind), weight_smooth, skin_validate,
                mesh_laplacian
uv/             projection, unwrap (xatlas), uv_metrics
optimization/   optimize (meshopt), meshlets, analyze, spatial, encode, strips,
                progressive (continuous LOD with serialization)
reconstruction/ point cloud → mesh
csg/            boolean (manifold)
io/             obj, stl, ply, gltf, fbx (read-only), vox (read-only)
```

### The rigging + animation pipeline (multi-file flow)

This is the one area where the separation of modules matters because it spans many files. The end-to-end authoring flow is:

1. **Landmarks** (`rigging/landmark_detect.h` *or* user-supplied) → `Landmarks` dict.
2. **`rigging/rig_spec.h`** declares the skeleton template + landmark schema. Bundled specs: `builtinHumanoidSpec()`, `builtinQuadrupedSpec()`, `builtinHexapodSpec()`, `builtinOctopodSpec()`. Bone positions are expressions (`"landmark:wrist_L"`, `"mid:A,B"`, `"lerp:A,B,t"`, `"offset:A,dx,dy,dz"`) resolved against the landmark dict.
3. **`rigging/auto_rig.h`** — `autoRig(mesh, spec, landmarks, wopts)` is the one-call driver. It runs `skeleton_fit` → weighting (dispatched by `WeightingOptions::method`) → `weight_smooth`.
4. **Weighting methods** (`rigging/weighting.h` dispatcher):
   - `BoneHeat` (`bone_heat.cpp`) — default for manifold meshes.
   - `BBW` (`bbw.cpp`) — opt-in, uses OSQP QP solver per bone.
   - `VoxelBind` (`voxel_bind.cpp`) — fallback for non-manifold / messy topology.
   - `Auto` — `autoSelectWeightingMethod(mesh)` picks BoneHeat vs VoxelBind based on manifoldness.
5. **Runtime animation** (`animation/pose.h`): `bindPose` → `evaluateAnimation` → `blendPoses` → `computeSkinningMatrices` → `manipulation/skin.h::applySkinning`. IK (`ik.h`) mutates `Pose` in place between evaluation and skinning.

The strategy doc at `docs/auto-rig-strategy.md` is the north star for this subsystem.

### Optional-by-submodule pattern

Every third-party dependency (meshoptimizer, V-HACD, tinygltf, par_shapes, xatlas, manifold, OpenFBX, OSQP) is a git submodule under `third_party/`. The top-level `CMakeLists.txt` gates each by an `EXISTS` check and sets `BROMESH_HAS_<DEP>` accordingly, forwarded as a public compile definition.

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
