#pragma once

// Shared test harness for bromesh_test. Every split test file includes
// this; test_main.cpp owns the counter definitions and the entry point.

#include "bromesh/mesh_data.h"
#include "bromesh/isosurface/marching_cubes.h"
#include "bromesh/isosurface/surface_nets.h"
#include "bromesh/isosurface/dual_contouring.h"
#include "bromesh/isosurface/transvoxel.h"
#include "bromesh/voxel/greedy_mesh.h"
#include "bromesh/primitives/primitives.h"
#include "bromesh/primitives/par_primitives.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/simplify.h"
#include "bromesh/manipulation/split_components.h"
#include "bromesh/manipulation/weld.h"
#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/repair.h"
#include "bromesh/manipulation/subdivide.h"
#include "bromesh/manipulation/smooth.h"
#include "bromesh/manipulation/remesh.h"
#include "bromesh/manipulation/skin.h"
#include "bromesh/manipulation/skin_transfer.h"
#include "bromesh/manipulation/shrinkwrap.h"
#include "bromesh/manipulation/transform.h"
#include "bromesh/analysis/bbox.h"
#include "bromesh/analysis/convex_decomposition.h"
#include "bromesh/analysis/bake.h"
#include "bromesh/analysis/bake_texture.h"
#include "bromesh/analysis/bake_transfer.h"
#include "bromesh/analysis/sample.h"
#include "bromesh/analysis/raycast.h"
#include "bromesh/analysis/bvh.h"
#include "bromesh/analysis/intersect.h"
#include "bromesh/animation/pose.h"
#include "bromesh/animation/ik.h"
#include "bromesh/animation/retarget.h"
#include "bromesh/animation/locomotion.h"
#include "bromesh/uv/projection.h"
#include "bromesh/uv/unwrap.h"
#include "bromesh/uv/uv_metrics.h"
#include "bromesh/optimization/optimize.h"
#include "bromesh/optimization/meshlets.h"
#include "bromesh/optimization/analyze.h"
#include "bromesh/optimization/spatial.h"
#include "bromesh/optimization/encode.h"
#include "bromesh/optimization/strips.h"
#include "bromesh/optimization/progressive.h"
#include "bromesh/reconstruction/reconstruct.h"
#include "bromesh/csg/boolean.h"
#include "bromesh/io/obj.h"
#include "bromesh/io/vox.h"
#include "bromesh/io/stl.h"
#include "bromesh/io/gltf.h"
#include "bromesh/io/ply.h"
#include "bromesh/io/fbx.h"
#include "bromesh/rigging/rig_spec.h"
#include "bromesh/rigging/landmarks.h"
#include "bromesh/rigging/skeleton_fit.h"
#include "bromesh/rigging/voxel_bind.h"
#include "bromesh/rigging/auto_rig.h"
#include "bromesh/rigging/landmark_detect.h"
#include "bromesh/rigging/weight_smooth.h"
#include "bromesh/rigging/bone_heat.h"
#include "bromesh/rigging/bbw.h"
#include "bromesh/rigging/weighting.h"
#include "bromesh/rigging/mesh_laplacian.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>

extern int tests_run;
extern int tests_passed;

// Test registry. Tests register themselves at static-init time via the TEST
// macro, but execution is deferred to main() — running them during static
// init leads to cross-TU init-order fiascos on Windows (e.g. calling into
// meshoptimizer before its own file-scope lookup tables have been built).
using TestFn = void(*)();
struct TestEntry { const char* name; TestFn fn; };
std::vector<TestEntry>& testRegistry();

struct TestRegistrar {
    TestRegistrar(const char* name, TestFn fn) {
        testRegistry().push_back({name, fn});
    }
};

#define TEST(name) \
    static void test_##name(); \
    static TestRegistrar reg_##name(#name, &test_##name); \
    static void test_##name()

#define ASSERT(cond, msg) do { \
    tests_run++; \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
    } else { \
        tests_passed++; \
    } \
} while(0)
