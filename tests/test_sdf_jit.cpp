#include "test_framework.h"
#include "bromesh/isosurface/jit/sdf_node.h"
#include "bromesh/isosurface/jit/sdf_compiler.h"
#include "bromesh/isosurface/jit/jit_mesher.h"
#include "bromesh/analysis/bbox.h"

#include <chrono>
#include <cmath>
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
