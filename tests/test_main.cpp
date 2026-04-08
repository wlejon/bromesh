#include "bromesh/mesh_data.h"
#include "bromesh/isosurface/marching_cubes.h"
#include "bromesh/isosurface/surface_nets.h"
#include "bromesh/isosurface/dual_contouring.h"
#include "bromesh/voxel/greedy_mesh.h"
#include "bromesh/primitives/primitives.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/simplify.h"
#include "bromesh/manipulation/split_components.h"
#include "bromesh/manipulation/weld.h"
#include "bromesh/analysis/bbox.h"
#include "bromesh/analysis/convex_decomposition.h"
#include "bromesh/uv/projection.h"
#include "bromesh/optimization/optimize.h"
#include "bromesh/io/obj.h"
#include "bromesh/io/vox.h"
#include "bromesh/io/stl.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) \
    static void test_##name(); \
    struct test_reg_##name { test_reg_##name() { test_##name(); } } reg_##name; \
    static void test_##name()

#define ASSERT(cond, msg) do { \
    tests_run++; \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
    } else { \
        tests_passed++; \
    } \
} while(0)

// --- Smoke tests: verify all headers compile and stubs link ---

TEST(mesh_data_basics) {
    bromesh::MeshData m;
    ASSERT(m.empty(), "new mesh should be empty");
    ASSERT(m.vertexCount() == 0, "vertex count should be 0");
    ASSERT(m.triangleCount() == 0, "triangle count should be 0");

    m.positions = {0,0,0, 1,0,0, 0,1,0};
    m.indices = {0, 1, 2};
    ASSERT(m.vertexCount() == 3, "should have 3 vertices");
    ASSERT(m.triangleCount() == 1, "should have 1 triangle");
    ASSERT(!m.hasNormals(), "no normals yet");

    m.clear();
    ASSERT(m.empty(), "cleared mesh should be empty");
}

TEST(bbox_basics) {
    bromesh::BBox b;
    b.min[0] = -1; b.min[1] = -2; b.min[2] = -3;
    b.max[0] =  1; b.max[1] =  2; b.max[2] =  3;
    ASSERT(std::fabs(b.centerX()) < 0.001f, "center X");
    ASSERT(std::fabs(b.centerY()) < 0.001f, "center Y");
    ASSERT(std::fabs(b.centerZ()) < 0.001f, "center Z");
    ASSERT(std::fabs(b.extentX() - 1.0f) < 0.001f, "extent X");
}

TEST(stubs_link) {
    // Just verify all stub functions link without crashing
    float field[8] = {-1,-1,-1,-1, 1,1,1,1};
    auto mc = bromesh::marchingCubes(field, 2, 2, 2);
    ASSERT(mc.empty(), "stub marching cubes returns empty");

    auto sn = bromesh::surfaceNets(field, 2, 2, 2);
    ASSERT(sn.empty(), "stub surface nets returns empty");

    auto dc = bromesh::dualContour(field, 2, 2, 2);
    ASSERT(dc.empty(), "stub dual contour returns empty");

    uint8_t voxels[8] = {1,1,1,1, 0,0,0,0};
    auto gm = bromesh::greedyMesh(voxels, 2, 2, 2);
    ASSERT(gm.empty(), "stub greedy mesh returns empty");

    auto b = bromesh::box(1, 1, 1);
    ASSERT(b.empty(), "stub box returns empty");

    auto s = bromesh::sphere(1);
    ASSERT(s.empty(), "stub sphere returns empty");
}

int main() {
    std::printf("bromesh tests: %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
