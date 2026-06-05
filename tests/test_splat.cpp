// Gaussian-splat I/O: round-trip a synthesized cloud through saveSplatPLY /
// loadSplatPLY, verify SH-degree inference + activation handling, the ASCII
// reader path, and that a non-splat PLY is rejected.

#include "test_framework.h"
#include "bromesh/gaussian_splat.h"
#include "bromesh/io/splat_ply.h"
#include "bromesh/io/ply.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

// Build a small deterministic cloud at the given SH degree. Values are
// render-ready (linear scale, [0,1] opacity, unit quats).
bromesh::GaussianSplatCloud makeCloud(int degree, size_t n) {
    bromesh::GaussianSplatCloud c;
    c.shDegree = degree;
    c.reserve(n);
    const int stride = c.shStride();
    for (size_t i = 0; i < n; ++i) {
        float fi = static_cast<float>(i);
        c.positions.push_back(fi);
        c.positions.push_back(fi * 2.0f);
        c.positions.push_back(-fi);
        c.scales.push_back(0.10f + 0.01f * fi);
        c.scales.push_back(0.20f + 0.01f * fi);
        c.scales.push_back(0.30f + 0.01f * fi);
        // A non-trivial unit quaternion: rotate about Z by a small angle.
        float ang = 0.1f * fi;
        c.rotations.push_back(0.0f);
        c.rotations.push_back(0.0f);
        c.rotations.push_back(std::sin(ang * 0.5f));
        c.rotations.push_back(std::cos(ang * 0.5f));
        c.opacities.push_back(0.25f + 0.5f * (fi / static_cast<float>(n)));
        for (int k = 0; k < stride; ++k)
            c.sh.push_back(0.01f * static_cast<float>(k) + 0.001f * fi);
    }
    return c;
}

bool nearf(float a, float b, float tol = 1e-4f) { return std::fabs(a - b) <= tol; }

} // namespace

TEST(splat_ply_binary_roundtrip_deg3) {
    auto cloud = makeCloud(3, 5);
    const char* path = "test_splat_deg3.ply";
    ASSERT(bromesh::saveSplatPLY(cloud, path), "saveSplatPLY deg3");

    auto loaded = bromesh::loadSplatPLY(path);
    ASSERT(loaded.count() == cloud.count(), "splat: vertex count round-trips");
    ASSERT(loaded.shDegree == 3, "splat: degree inferred as 3");
    ASSERT(loaded.shStride() == cloud.shStride(), "splat: sh stride matches");

    bool posOk = true, scaleOk = true, opOk = true, rotOk = true, shOk = true;
    for (size_t v = 0; v < cloud.count(); ++v) {
        for (int a = 0; a < 3; ++a) {
            posOk &= nearf(loaded.positions[v * 3 + a], cloud.positions[v * 3 + a]);
            scaleOk &= nearf(loaded.scales[v * 3 + a], cloud.scales[v * 3 + a]);
        }
        opOk &= nearf(loaded.opacities[v], cloud.opacities[v]);
        for (int a = 0; a < 4; ++a)
            rotOk &= nearf(loaded.rotations[v * 4 + a], cloud.rotations[v * 4 + a]);
        for (int k = 0; k < cloud.shStride(); ++k)
            shOk &= nearf(loaded.sh[v * cloud.shStride() + k],
                          cloud.sh[v * cloud.shStride() + k]);
    }
    ASSERT(posOk, "splat: positions round-trip");
    ASSERT(scaleOk, "splat: scales round-trip (exp/log)");
    ASSERT(opOk, "splat: opacities round-trip (sigmoid/logit)");
    ASSERT(rotOk, "splat: rotations round-trip (w-first <-> xyzw)");
    ASSERT(shOk, "splat: SH coeffs round-trip (interleave <-> channel-major)");

    std::remove(path);
}

TEST(splat_ply_roundtrip_deg0) {
    auto cloud = makeCloud(0, 3);
    const char* path = "test_splat_deg0.ply";
    ASSERT(bromesh::saveSplatPLY(cloud, path), "saveSplatPLY deg0");
    auto loaded = bromesh::loadSplatPLY(path);
    ASSERT(loaded.shDegree == 0, "splat: degree 0 (no f_rest)");
    ASSERT(loaded.shStride() == 3, "splat: deg0 stride is 3 (DC only)");
    ASSERT(loaded.count() == 3, "splat: deg0 count");
    std::remove(path);
}

TEST(splat_bounds) {
    auto cloud = makeCloud(0, 4); // positions (i, 2i, -i) for i in 0..3
    auto b = cloud.bounds();
    ASSERT(nearf(b.min.x, 0.0f) && nearf(b.max.x, 3.0f), "splat bounds x");
    ASSERT(nearf(b.min.y, 0.0f) && nearf(b.max.y, 6.0f), "splat bounds y");
    ASSERT(nearf(b.min.z, -3.0f) && nearf(b.max.z, 0.0f), "splat bounds z");
}

TEST(splat_ply_ascii_load) {
    // Minimal degree-0 ASCII 3DGS ply: one splat. rot is w-first (rot_0=w=1).
    // opacity 0 -> sigmoid 0.5; scale_* 0 -> exp 1.0.
    const char* path = "test_splat_ascii.ply";
    FILE* f = std::fopen(path, "w");
    ASSERT(f != nullptr, "open ascii splat ply");
    if (!f) return;
    std::fputs(
        "ply\n"
        "format ascii 1.0\n"
        "element vertex 1\n"
        "property float x\nproperty float y\nproperty float z\n"
        "property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n"
        "property float opacity\n"
        "property float scale_0\nproperty float scale_1\nproperty float scale_2\n"
        "property float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\n"
        "end_header\n"
        "1 2 3 0.5 0.6 0.7 0 0 0 0 1 0 0 0\n", f);
    std::fclose(f);

    auto c = bromesh::loadSplatPLY(path);
    ASSERT(c.count() == 1, "ascii splat: 1 splat");
    ASSERT(c.shDegree == 0, "ascii splat: degree 0");
    ASSERT(nearf(c.positions[0], 1.0f) && nearf(c.positions[1], 2.0f) &&
           nearf(c.positions[2], 3.0f), "ascii splat: position");
    ASSERT(nearf(c.sh[0], 0.5f) && nearf(c.sh[1], 0.6f) && nearf(c.sh[2], 0.7f),
           "ascii splat: DC color");
    ASSERT(nearf(c.opacities[0], 0.5f), "ascii splat: opacity sigmoid(0)=0.5");
    ASSERT(nearf(c.scales[0], 1.0f), "ascii splat: scale exp(0)=1");
    ASSERT(nearf(c.rotations[3], 1.0f), "ascii splat: quat w=1 (identity)");
    std::remove(path);
}

TEST(splat_ply_rejects_plain_mesh) {
    // A plain triangle mesh PLY has no splat properties -> loadSplatPLY empty.
    const char* path = "test_plain_mesh.ply";
    FILE* f = std::fopen(path, "w");
    ASSERT(f != nullptr, "open plain mesh ply");
    if (!f) return;
    std::fputs(
        "ply\n"
        "format ascii 1.0\n"
        "element vertex 3\n"
        "property float x\nproperty float y\nproperty float z\n"
        "element face 1\n"
        "property list uchar uint vertex_indices\n"
        "end_header\n"
        "0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n", f);
    std::fclose(f);

    auto c = bromesh::loadSplatPLY(path);
    ASSERT(c.empty(), "plain mesh rejected by loadSplatPLY");
    // And it still loads fine as a mesh.
    auto m = bromesh::loadPLY(path);
    ASSERT(m.vertexCount() == 3, "plain mesh still loads via loadPLY");
    std::remove(path);
}
