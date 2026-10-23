#include "test_framework.h"
#include "bromesh/isosurface/jit/sdf_node.h"
#include "bromesh/isosurface/jit/sdf_compiler.h"
#include "bromesh/isosurface/jit/jit_mesher.h"
#include "bromesh/analysis/bbox.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/target/spirv/spirv_ir.hpp>
#endif

#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>

TEST(sdf_jit_graph_construction_and_hash) {
    bromesh::SdfGraph g1;
    int s1 = g1.sphere(1.5f);
    ASSERT(s1 == 0, "sphere node is id 0");
    ASSERT(g1.root() == 0, "root is node 0");
    ASSERT(g1.size() == 1, "graph size 1");

    uint64_t h1 = g1.computeHash();
    ASSERT(h1 != 0, "hash non-zero");

    // Same graph should give identical hash
    bromesh::SdfGraph g2;
    g2.sphere(1.5f);
    ASSERT(g2.computeHash() == h1, "identical graph gives same hash");

    // Different radius should change hash
    bromesh::SdfGraph g3;
    g3.sphere(1.6f);
    ASSERT(g3.computeHash() != h1, "different parameter gives different hash");

    // Adding more nodes
    bromesh::SdfGraph gCsg;
    int a = gCsg.sphere(1.0f);
    int b = gCsg.box({0.8f, 0.8f, 0.8f});
    int u = gCsg.opSmoothUnion(a, b, 0.2f);
    ASSERT(u == 2, "smooth union is node 2");
    ASSERT(gCsg.root() == 2, "root defaults to last node");
    ASSERT(gCsg.computeHash() != h1, "csg hash is unique");

    gCsg.setRoot(0);
    ASSERT(gCsg.root() == 0, "explicit root change works");
    ASSERT(gCsg.computeHash() != h1, "hash changes when root changes");
}

TEST(sdf_jit_volume_evaluation_primitives) {
    // 1. Sphere volume
    bromesh::SdfGraph gSphere;
    gSphere.sphere(1.0f);

    bromath::AABB3 bounds{{-2.0f, -2.0f, -2.0f}, {2.0f, 2.0f, 2.0f}};
    int dim = 33; // center voxel is exactly at index 16
    auto compiled = bromesh::JitSdfCompiler::instance().compile(gSphere);
    ASSERT(compiled.isJit && compiled.gridFn != nullptr && compiled.pointFn != nullptr, "JIT compilation must succeed and produce valid native function pointers");
    auto vol = bromesh::JitSdfCompiler::instance().evaluateVolume(gSphere, dim, dim, dim, bounds);

    ASSERT(vol.dimX == dim && vol.dimY == dim && vol.dimZ == dim, "volume dimensions match");
    ASSERT(vol.field.size() == static_cast<size_t>(dim * dim * dim), "field size matches");

    // Center (0, 0, 0)
    float valCenter = vol.value(16, 16, 16);
    ASSERT(std::fabs(valCenter - (-1.0f)) < 0.05f, "center of sphere is -radius");

    // On surface (1, 0, 0) -> index (16 + 8, 16, 16) = (24, 16, 16) since step = 4/32 = 0.125
    float valSurface = vol.value(24, 16, 16);
    ASSERT(std::fabs(valSurface) < 0.05f, "surface of sphere is approx 0");

    // Outside (2, 0, 0) -> index (32, 16, 16)
    float valOutside = vol.value(32, 16, 16);
    ASSERT(std::fabs(valOutside - 1.0f) < 0.05f, "outside distance is approx 1.0");

    // 2. Box volume
    bromesh::SdfGraph gBox;
    gBox.box({1.0f, 1.0f, 1.0f});
    auto volBox = bromesh::JitSdfCompiler::instance().evaluateVolume(gBox, dim, dim, dim, bounds);
    ASSERT(volBox.value(16, 16, 16) < -0.9f, "box center is negative");
    ASSERT(volBox.value(32, 16, 16) > 0.9f, "box outside is positive");

    // 3. Smooth union of sphere + cylinder
    bromesh::SdfGraph gCombo;
    int sp = gCombo.sphere(1.0f);
    int cy = gCombo.cylinder(0.4f, 1.5f);
    gCombo.opSmoothUnion(sp, cy, 0.25f);

    auto volCombo = bromesh::JitSdfCompiler::instance().evaluateVolume(gCombo, dim, dim, dim, bounds);
    ASSERT(volCombo.value(16, 16, 16) < -0.9f, "smooth union center is inside");
    ASSERT(volCombo.value(16, 26, 16) < 0.0f, "cylinder extension is inside");
}

TEST(sdf_jit_mesher_marching_cubes_sphere) {
    bromesh::SdfGraph g;
    float r = 1.0f;
    g.sphere(r);

    bromath::AABB3 bounds{{-1.5f, -1.5f, -1.5f}, {1.5f, 1.5f, 1.5f}};
    int dim = 48;
    bromesh::MeshData mesh = bromesh::marchingCubesFromSDF(g, dim, dim, dim, bounds, 0.0f, true, true);

    ASSERT(!mesh.empty(), "MC sphere mesh is not empty");
    ASSERT(!mesh.indices.empty(), "MC sphere has indices");
    ASSERT(mesh.positions.size() == mesh.normals.size(), "Positions and normals have matching count");

    float volume = bromesh::computeVolume(mesh);
    float expectedVolume = (4.0f / 3.0f) * 3.14159265f * (r * r * r);
    ASSERT(volume > 0.0f, "Sphere volume must be positive");
    ASSERT(std::fabs(volume - expectedVolume) / expectedVolume < 0.08f, "Sphere volume matches analytical expectation within 8%");

    // Verify analytical normals: for sphere centered at origin, outward normal should equal pos / |pos|
    size_t numVerts = mesh.positions.size() / 3;
    float minDot = 1.0f;
    for (size_t i = 0; i < numVerts; ++i) {
        float px = mesh.positions[3 * i + 0];
        float py = mesh.positions[3 * i + 1];
        float pz = mesh.positions[3 * i + 2];
        float plen = std::sqrt(px * px + py * py + pz * pz);
        ASSERT(plen > 0.0f, "Vertex not at origin");
        float ex = px / plen;
        float ey = py / plen;
        float ez = pz / plen;

        float nx = mesh.normals[3 * i + 0];
        float ny = mesh.normals[3 * i + 1];
        float nz = mesh.normals[3 * i + 2];
        float dot = nx * ex + ny * ey + nz * ez;
        if (dot < minDot) minDot = dot;
    }
    ASSERT(minDot > 0.99f, "Analytical normals on sphere point exactly outward");
}

TEST(sdf_jit_mesher_surface_nets_sphere) {
    bromesh::SdfGraph g;
    float r = 1.0f;
    g.sphere(r);

    bromath::AABB3 bounds{{-1.5f, -1.5f, -1.5f}, {1.5f, 1.5f, 1.5f}};
    int dim = 48;
    bromesh::MeshData mesh = bromesh::surfaceNetsFromSDF(g, dim, dim, dim, bounds, 0.0f);

    ASSERT(!mesh.empty(), "Surface Nets sphere mesh is not empty");
    ASSERT(!mesh.indices.empty(), "Surface Nets sphere has indices");
    ASSERT(mesh.positions.size() == mesh.normals.size(), "Positions and normals have matching count");

    float volume = bromesh::computeVolume(mesh);
    float expectedVolume = (4.0f / 3.0f) * 3.14159265f * (r * r * r);
    ASSERT(volume > 0.0f, "Surface Nets volume must be positive");
    ASSERT(std::fabs(volume - expectedVolume) / expectedVolume < 0.10f, "Surface Nets sphere volume matches analytical expectation within 10%");

    // Verify analytical normals
    size_t numVerts = mesh.positions.size() / 3;
    float minDot = 1.0f;
    for (size_t i = 0; i < numVerts; ++i) {
        float px = mesh.positions[3 * i + 0];
        float py = mesh.positions[3 * i + 1];
        float pz = mesh.positions[3 * i + 2];
        float plen = std::sqrt(px * px + py * py + pz * pz);
        float ex = px / plen;
        float ey = py / plen;
        float ez = pz / plen;

        float nx = mesh.normals[3 * i + 0];
        float ny = mesh.normals[3 * i + 1];
        float nz = mesh.normals[3 * i + 2];
        float dot = nx * ex + ny * ey + nz * ez;
        if (dot < minDot) minDot = dot;
    }
    ASSERT(minDot > 0.98f, "Surface Nets analytical normals point outward");
}

TEST(sdf_jit_mesher_smooth_union_csg) {
    bromesh::SdfGraph g;
    int s = g.sphere(0.8f);
    int c = g.cylinder(0.4f, 1.2f);
    g.opSmoothUnion(s, c, 0.2f);

    bromath::AABB3 bounds{{-1.6f, -1.6f, -1.6f}, {1.6f, 1.6f, 1.6f}};
    int dim = 48;

    bromesh::MeshData mcMesh = bromesh::marchingCubesFromSDF(g, dim, dim, dim, bounds, 0.0f, true, true);
    ASSERT(!mcMesh.empty(), "MC smooth union mesh generated");
    float mcVol = bromesh::computeVolume(mcMesh);
    ASSERT(mcVol > 0.0f, "MC smooth union has positive volume");

    bromesh::MeshData snMesh = bromesh::surfaceNetsFromSDF(g, dim, dim, dim, bounds, 0.0f);
    ASSERT(!snMesh.empty(), "SN smooth union mesh generated");
    float snVol = bromesh::computeVolume(snMesh);
    ASSERT(snVol > 0.0f, "SN smooth union has positive volume");
}

TEST(sdf_jit_benchmark_128_cube) {
    bromesh::SdfGraph g;
    int s = g.sphere(0.85f);
    int b = g.box({0.6f, 0.6f, 0.6f});
    int u = g.opSmoothUnion(s, b, 0.15f);
    int noise = g.noise3D({3.0f, 3.0f, 3.0f}, 0.05f);
    g.displace(u, noise);

    bromath::AABB3 bounds{{-1.5f, -1.5f, -1.5f}, {1.5f, 1.5f, 1.5f}};
    const int dim = 128;
    const size_t totalVoxels = static_cast<size_t>(dim) * dim * dim; // 2,097,152 samples

    // Warm up / compilation
    auto warmupVol = bromesh::JitSdfCompiler::instance().evaluateVolume(g, 16, 16, 16, bounds);
    ASSERT(!warmupVol.field.empty(), "Warmup succeeded");

    // Timed run
    auto t0 = std::chrono::high_resolution_clock::now();
    auto vol = bromesh::JitSdfCompiler::instance().evaluateVolume(g, dim, dim, dim, bounds);
    auto t1 = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double mvoxPerSec = (static_cast<double>(totalVoxels) / 1e6) / (ms / 1000.0);

    ASSERT(vol.field.size() == totalVoxels, "Evaluated exactly 128^3 voxels");
    std::cout << "\n  [BENCHMARK] JIT 128^3 volume evaluation ("
              << totalVoxels << " voxels): "
              << ms << " ms (" << mvoxPerSec << " MVoxels/sec)" << std::endl;

    // Sanity check that execution is very fast (e.g. < 500 ms for 2M samples in JIT)
    ASSERT(ms < 2000.0, "128^3 evaluation must complete under 2 seconds");
}

TEST(sdf_jit_parallel_vs_sequential_equivalence_97) {
    bromesh::SdfGraph g;
    int s = g.sphere(0.85f);
    int b = g.roundedBox({0.6f, 0.5f, 0.7f}, 0.1f);
    int u = g.opSmoothUnion(s, b, 0.15f);
    int t = g.torus(0.7f, 0.2f);
    int sub = g.opSmoothSubtraction(u, t, 0.1f);
    int noise = g.noise3D({2.5f, 2.5f, 2.5f}, 0.04f);
    g.displace(sub, noise);

    bromath::AABB3 bounds{{-1.5f, -1.5f, -1.5f}, {1.5f, 1.5f, 1.5f}};
    const int dim = 36;

    auto volParallel = bromesh::SdfCompiler::instance().evaluateVolume(g, dim, dim, dim, bounds, -1, true);
    auto volSequential = bromesh::SdfCompiler::instance().evaluateVolume(g, dim, dim, dim, bounds, -1, false);

    ASSERT(volParallel.field.size() == volSequential.field.size(), "Field sizes match");
    ASSERT(!volParallel.field.empty(), "Volume field is non-empty");

    for (size_t i = 0; i < volParallel.field.size(); ++i) {
        uint32_t bitsPar = 0;
        uint32_t bitsSeq = 0;
        std::memcpy(&bitsPar, &volParallel.field[i], sizeof(uint32_t));
        std::memcpy(&bitsSeq, &volSequential.field[i], sizeof(uint32_t));
        ASSERT(bitsPar == bitsSeq, "Parallel Z-slice evaluation must be bit-identical to sequential evaluation");
    }
}

TEST(sdf_jit_gpu_spirv_and_ptx_codegen_99) {
#if BROMESH_HAS_BRASS_JIT
    bromesh::SdfGraph g;
    int s = g.sphere(1.0f);
    int b = g.box({0.75f, 0.75f, 0.75f});
    int u = g.opSmoothUnion(s, b, 0.2f);
    int c = g.cylinder(0.35f, 1.2f);
    int diff = g.opSubtraction(u, c);
    int xform = g.translate(diff, {0.1f, -0.2f, 0.05f});
    int noise = g.noise3D({2.0f, 2.0f, 2.0f}, 0.03f);
    g.displace(xform, noise);

    std::vector<uint32_t> spv = bromesh::SdfCompiler::instance().compile_spirv(g);
    ASSERT(!spv.empty(), "SPIR-V binary must not be empty");
    ASSERT(spv[0] == 0x07230203u, "SPIR-V magic number must match 0x07230203");
    std::string spvErr;
    ASSERT(brass::spirv::verify_spirv_module(spv, &spvErr), "Generated SPIR-V compute shader must pass structural verification");

    std::string ptx = bromesh::SdfCompiler::instance().compile_ptx(g);
    ASSERT(!ptx.empty(), "PTX source must not be empty");
    ASSERT(ptx.find(".visible .entry eval_sdf_gpu") != std::string::npos, "PTX must contain .visible .entry eval_sdf_gpu");
    ASSERT(ptx.find("st.global") != std::string::npos, "PTX must contain global store to output volume grid");
#endif
}

TEST(sdf_jit_interval_macro_block_pruning_106) {
    bromesh::SdfGraph g;
    int s1 = g.translate(g.sphere(0.75f), {-0.6f, 0.0f, 0.0f});
    int b1 = g.translate(g.roundedBox({0.5f, 0.5f, 0.5f}, 0.08f), {0.6f, 0.0f, 0.0f});
    int u = g.opSmoothUnion(s1, b1, 0.15f);
    int cyl = g.cylinder(0.25f, 1.0f);
    g.opSubtraction(u, cyl);

    auto kernel = bromesh::SdfCompiler::instance().compile(g);
    ASSERT(kernel.blockBoundsFn != nullptr, "CompiledSdfKernel must provide JIT blockBoundsFn");

    // 1. Verify interval soundness over macro-blocks against dense samples
    bromath::AABB3 testBlock{{-1.5f, -1.5f, -1.5f}, {-1.1f, -1.1f, -1.1f}};
    float dMin = 0.0f, dMax = 0.0f;
    int domNode = -1;
    bromesh::evalBlockBounds(g, -1, testBlock, dMin, dMax, &domNode);

    float jitMin = 0.0f, jitMax = 0.0f;
    kernel.evalBlockBounds(testBlock, jitMin, jitMax);
    ASSERT(std::fabs(dMin - jitMin) < 1e-4f && std::fabs(dMax - jitMax) < 1e-4f,
           "JIT blockBoundsFn matches C++ interval bounds");
    ASSERT(dMin > 0.0f, "Far corner block must be provably outside (> 0)");

    for (int iz = 0; iz <= 4; ++iz) {
        for (int iy = 0; iy <= 4; ++iy) {
            for (int ix = 0; ix <= 4; ++ix) {
                float px = testBlock.min.x + (testBlock.max.x - testBlock.min.x) * (ix / 4.0f);
                float py = testBlock.min.y + (testBlock.max.y - testBlock.min.y) * (iy / 4.0f);
                float pz = testBlock.min.z + (testBlock.max.z - testBlock.min.z) * (iz / 4.0f);
                float d = kernel.pointFn(px, py, pz);
                ASSERT(d >= dMin - 1e-4f && d <= dMax + 1e-4f,
                       "Point sample must lie inside [dMin, dMax] interval");
            }
        }
    }

    // 2. Verify CSG branch pruning when one branch provably dominates another
    bromesh::SdfGraph gDisjoint;
    int leftSphere = gDisjoint.translate(gDisjoint.sphere(0.5f), {-2.0f, 0.0f, 0.0f});
    int rightSphere = gDisjoint.translate(gDisjoint.sphere(0.5f), {2.0f, 0.0f, 0.0f});
    int rootUnion = gDisjoint.opSmoothUnion(leftSphere, rightSphere, 0.2f);

    bromath::AABB3 leftBlock{{-2.5f, -0.5f, -0.5f}, {-1.5f, 0.5f, 0.5f}};
    int activeSubtree = -1;
    float lbMin = 0.0f, lbMax = 0.0f;
    bromesh::evalBlockBounds(gDisjoint, rootUnion, leftBlock, lbMin, lbMax, &activeSubtree);
    ASSERT(activeSubtree == leftSphere, "CSG pruning must select leftSphere as dominating branch in leftBlock");

    bromesh::SdfGraph prunedLeft = bromesh::pruneCsgGraphForBlock(gDisjoint, leftBlock);
    ASSERT(prunedLeft.size() < gDisjoint.size(), "Pruned graph has fewer nodes than full CSG graph");
    auto kFull = bromesh::SdfCompiler::instance().compile(gDisjoint);
    auto kPruned = bromesh::SdfCompiler::instance().compile(prunedLeft);
    ASSERT(std::fabs(kFull.pointFn(-2.0f, 0.1f, -0.1f) - kPruned.pointFn(-2.0f, 0.1f, -0.1f)) < 1e-6f,
           "Pruned CSG graph evaluates identically inside macro-block");

    // 3. Verify 100% identical mesh topology & geometry vs dense evaluation
    bromath::AABB3 bounds{{-2.0f, -2.0f, -2.0f}, {2.0f, 2.0f, 2.0f}};
    int dim = 64;
    bromesh::MeshData mcPruned = bromesh::marchingCubesFromSDF(g, dim, dim, dim, bounds, 0.0f, true, true, true);
    bromesh::MeshData mcDense = bromesh::marchingCubesFromSDF(g, dim, dim, dim, bounds, 0.0f, true, true, false);
    ASSERT(!mcPruned.empty(), "Pruned MC mesh is non-empty");
    ASSERT(mcPruned.positions.size() == mcDense.positions.size(), "MC pruned vs dense vertex count identical");
    ASSERT(mcPruned.indices.size() == mcDense.indices.size(), "MC pruned vs dense index count identical");
    for (size_t i = 0; i < mcPruned.positions.size(); ++i) {
        ASSERT(std::fabs(mcPruned.positions[i] - mcDense.positions[i]) < 1e-6f, "MC vertex positions match");
        ASSERT(std::fabs(mcPruned.normals[i] - mcDense.normals[i]) < 1e-6f, "MC vertex normals match");
    }
    for (size_t i = 0; i < mcPruned.indices.size(); ++i) {
        ASSERT(mcPruned.indices[i] == mcDense.indices[i], "MC triangle indices match");
    }

    bromesh::MeshData snPruned = bromesh::surfaceNetsFromSDF(g, dim, dim, dim, bounds, 0.0f, true);
    bromesh::MeshData snDense = bromesh::surfaceNetsFromSDF(g, dim, dim, dim, bounds, 0.0f, false);
    ASSERT(snPruned.positions.size() == snDense.positions.size(), "SN pruned vs dense vertex count identical");
    ASSERT(snPruned.indices.size() == snDense.indices.size(), "SN pruned vs dense index count identical");
    for (size_t i = 0; i < snPruned.positions.size(); ++i) {
        ASSERT(std::fabs(snPruned.positions[i] - snDense.positions[i]) < 1e-6f, "SN vertex positions match");
    }

    // 4. Measure pruning ratio and speedup on 128^3 and 256^3 volumes
    auto t0 = std::chrono::high_resolution_clock::now();
    auto vol128Pruned = bromesh::SdfCompiler::instance().evaluateVolume(g, 128, 128, 128, bounds, -1, true, 0.0f, true);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto vol128Dense = bromesh::SdfCompiler::instance().evaluateVolume(g, 128, 128, 128, bounds, -1, true, 0.0f, false);
    auto t2 = std::chrono::high_resolution_clock::now();

    double ms128Pruned = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms128Dense = std::chrono::duration<double, std::milli>(t2 - t1).count();
    ASSERT(vol128Pruned.activeBlockCount < vol128Pruned.totalBlockCount / 2,
           "More than 50% of 8x8x8 macro-blocks must be pruned on 128^3 volume");

    auto t3 = std::chrono::high_resolution_clock::now();
    auto vol256Pruned = bromesh::SdfCompiler::instance().evaluateVolume(g, 256, 256, 256, bounds, -1, true, 0.0f, true);
    auto t4 = std::chrono::high_resolution_clock::now();
    double ms256Pruned = std::chrono::duration<double, std::milli>(t4 - t3).count();
    ASSERT(vol256Pruned.activeBlockCount < vol256Pruned.totalBlockCount / 3,
           "More than 66% of 8x8x8 macro-blocks must be pruned on 256^3 volume");

    std::cout << "\n  [PRUNING #106] 128^3 active blocks: " << vol128Pruned.activeBlockCount
              << "/" << vol128Pruned.totalBlockCount << " (" << ms128Pruned << " ms pruned vs "
              << ms128Dense << " ms dense); 256^3 active blocks: " << vol256Pruned.activeBlockCount
              << "/" << vol256Pruned.totalBlockCount << " (" << ms256Pruned << " ms)" << std::endl;
}

TEST(sdf_jit_symbolic_autodiff_and_contact_107) {
    // 1. Verify exact analytical gradient on Sphere, Torus, and Transforms
    bromesh::SdfGraph gSphere;
    gSphere.sphere(1.2f);
    auto kSphere = bromesh::SdfCompiler::instance().compile(gSphere);
    ASSERT(kSphere.gradFn != nullptr, "CompiledSdfKernel must provide 1-pass symbolic gradFn");

    float gx = 0.0f, gy = 0.0f, gz = 0.0f;
    float d = kSphere.gradFn(0.6f, -0.8f, 1.0f, &gx, &gy, &gz);
    float r = std::sqrt(0.6f * 0.6f + 0.8f * 0.8f + 1.0f * 1.0f);
    ASSERT(std::fabs(d - (r - 1.2f)) < 1e-5f, "gradFn returns exact SDF value");
    ASSERT(std::fabs(gx - 0.6f / r) < 1e-5f, "gradFn exact d/dx on sphere");
    ASSERT(std::fabs(gy - (-0.8f / r)) < 1e-5f, "gradFn exact d/dy on sphere");
    ASSERT(std::fabs(gz - 1.0f / r) < 1e-5f, "gradFn exact d/dz on sphere");

    // 2. Verify O(1) contact query (distance, normal, closestSurfacePoint)
    bromesh::SdfContact contact = kSphere.queryContact(0.0f, 2.0f, 0.0f);
    ASSERT(std::fabs(contact.distance - 0.8f) < 1e-5f, "Contact distance is 2.0 - 1.2 = 0.8");
    ASSERT(std::fabs(contact.normal.x) < 1e-5f &&
           std::fabs(contact.normal.y - 1.0f) < 1e-5f &&
           std::fabs(contact.normal.z) < 1e-5f, "Contact normal is (0, 1, 0)");
    ASSERT(std::fabs(contact.closestSurfacePoint.x) < 1e-5f &&
           std::fabs(contact.closestSurfacePoint.y - 1.2f) < 1e-5f &&
           std::fabs(contact.closestSurfacePoint.z) < 1e-5f,
           "Closest surface point lies on sphere of radius 1.2");

    // 3. Sharp CSG crease preservation: two intersecting planes p1 = x, p2 = y
    // Near the crease (e.g. x = 0.001, y = 0.0), Union(p1, p2) = min(x, y) = y (since 0.0 < 0.001).
    // Symbolic auto-diff gives exact normal (0, 1, 0) right up to the crease, whereas finite
    // differences with eps > 0.001 straddle x = y and blur the normal!
    bromesh::SdfGraph gCrease;
    int p1 = gCrease.plane({1.0f, 0.0f, 0.0f}, 0.0f);
    int p2 = gCrease.plane({0.0f, 1.0f, 0.0f}, 0.0f);
    gCrease.opUnion(p1, p2);
    auto kCrease = bromesh::SdfCompiler::instance().compile(gCrease);
    auto creaseContact = kCrease.queryContact(0.0005f, 0.0f, 0.0f);
    ASSERT(std::fabs(creaseContact.normal.x - 0.0f) < 1e-5f &&
           std::fabs(creaseContact.normal.y - 1.0f) < 1e-5f,
           "Symbolic auto-diff preserves sharp CSG crease normal without finite-difference blurring");
}

TEST(sdf_jit_custom_expression_gyroid_and_gpu_108) {
    bromesh::SdfGraph g;
    // Gyroid TPMS lattice intersected with a bounding sphere
    int gyroid = g.customExpr(
        "sin(x * 4.0) * cos(y * 4.0) + sin(y * 4.0) * cos(z * 4.0) + sin(z * 4.0) * cos(x * 4.0) - 0.2",
        4.0f
    );
    int boundSphere = g.sphere(1.3f);
    g.opIntersection(gyroid, boundSphere);

    auto kernel = bromesh::SdfCompiler::instance().compile(g);
    ASSERT(kernel.isJit && kernel.pointFn && kernel.gradFn && kernel.blockBoundsFn,
           "Custom expression SDF compiles to native CPU point, grid, grad, and interval kernels");

    // Verify point evaluation and symbolic gradient against analytical formula for the gyroid alone
    auto kGyroid = bromesh::SdfCompiler::instance().compile(g, gyroid);
    float x = 0.3f, y = -0.4f, z = 0.25f;
    float gx = 0.0f, gy = 0.0f, gz = 0.0f;
    float d = kGyroid.gradFn(x, y, z, &gx, &gy, &gz);

    float expectedD = std::sin(x * 4.0f) * std::cos(y * 4.0f) +
                      std::sin(y * 4.0f) * std::cos(z * 4.0f) +
                      std::sin(z * 4.0f) * std::cos(x * 4.0f) - 0.2f;
    float expectedGx = 4.0f * std::cos(x * 4.0f) * std::cos(y * 4.0f) -
                       4.0f * std::sin(z * 4.0f) * std::sin(x * 4.0f);
    float expectedGy = -4.0f * std::sin(x * 4.0f) * std::sin(y * 4.0f) +
                       4.0f * std::cos(y * 4.0f) * std::cos(z * 4.0f);
    float expectedGz = -4.0f * std::sin(y * 4.0f) * std::sin(z * 4.0f) +
                       4.0f * std::cos(z * 4.0f) * std::cos(x * 4.0f);

    ASSERT(std::fabs(d - expectedD) < 1e-4f, "Gyroid customExpr value matches analytical formula");
    ASSERT(std::fabs(gx - expectedGx) < 1e-3f, "Gyroid customExpr d/dx matches analytical derivative");
    ASSERT(std::fabs(gy - expectedGy) < 1e-3f, "Gyroid customExpr d/dy matches analytical derivative");
    ASSERT(std::fabs(gz - expectedGz) < 1e-3f, "Gyroid customExpr d/dz matches analytical derivative");

    // Verify JS arrow function syntax with Math.* intrinsics
    bromesh::SdfGraph gJs;
    gJs.customExpr("(x, y, z) => Math.sqrt(x*x + y*y + z*z) - 1.0");
    auto kJs = bromesh::SdfCompiler::instance().compile(gJs);
    ASSERT(std::fabs(kJs.pointFn(2.0f, 0.0f, 0.0f) - 1.0f) < 1e-5f, "JS arrow expression evaluates accurately");

    // Verify Marching Cubes meshing of the bounded Gyroid
    bromath::AABB3 bounds{{-1.6f, -1.6f, -1.6f}, {1.6f, 1.6f, 1.6f}};
    bromesh::MeshData mesh = bromesh::marchingCubesFromSDF(g, 40, 40, 40, bounds, 0.0f, true, true);
    ASSERT(!mesh.empty() && mesh.positions.size() > 300, "Bounded Gyroid custom SDF produces valid triangle mesh");

#if BROMESH_HAS_BRASS_JIT
    // Verify SPIR-V and PTX GPU codegen for custom expressions
    std::vector<uint32_t> spv = bromesh::SdfCompiler::instance().compile_spirv(g);
    ASSERT(!spv.empty() && spv[0] == 0x07230203u, "CustomExpr SPIR-V binary generated");
    std::string spvErr;
    ASSERT(brass::spirv::verify_spirv_module(spv, &spvErr), "CustomExpr SPIR-V passes verification");

    std::string ptx = bromesh::SdfCompiler::instance().compile_ptx(g);
    ASSERT(ptx.find(".visible .entry eval_sdf_gpu") != std::string::npos, "CustomExpr PTX entry generated");
    ASSERT(ptx.find("sin.approx.f32") != std::string::npos || ptx.find("cos.approx.f32") != std::string::npos,
           "CustomExpr PTX emits native hardware sin/cos instructions");
#endif
}
