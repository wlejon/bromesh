// Gaussian-splat I/O: round-trip a synthesized cloud through saveSplatPLY /
// loadSplatPLY, verify SH-degree inference + activation handling, the ASCII
// reader path, and that a non-splat PLY is rejected.

#include "test_framework.h"
#include "bromesh/gaussian_splat.h"
#include "bromesh/io/splat_ply.h"
#include "bromesh/io/ply.h"
#include "bromesh/manipulation/splat_ops.h"
#include "bromesh/primitives/primitives.h"

#include <algorithm>
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

TEST(splat_ply_degree1_channel_major_sh_and_activations) {
    // Degree-1 ASCII 3DGS ply with 1 vertex.
    // Tests:
    // 1. Scale log -> linear activation (exp)
    // 2. Opacity logit -> linear activation (sigmoid)
    // 3. SH rest channel-major to interleaved transpose
    // 4. Save and reload roundtrip
    const char* path = "test_splat_deg1_analytic.ply";
    FILE* f = std::fopen(path, "w");
    ASSERT(f != nullptr, "open ascii deg1 splat ply");
    if (!f) return;

    std::fputs(
        "ply\n"
        "format ascii 1.0\n"
        "element vertex 1\n"
        "property float x\nproperty float y\nproperty float z\n"
        "property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n"
        "property float f_rest_0\nproperty float f_rest_1\nproperty float f_rest_2\n"
        "property float f_rest_3\nproperty float f_rest_4\nproperty float f_rest_5\n"
        "property float f_rest_6\nproperty float f_rest_7\nproperty float f_rest_8\n"
        "property float opacity\n"
        "property float scale_0\nproperty float scale_1\nproperty float scale_2\n"
        "property float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\n"
        "end_header\n"
        "1.0 2.0 3.0 "                         // x y z
        "0.1 0.2 0.3 "                         // f_dc_0..2
        "1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 " // f_rest_0..8 (R1 R2 R3 G1 G2 G3 B1 B2 B3)
        "1.386294 "                            // opacity logit(0.8) ~ 1.386294
        "0.693147 -0.693147 0.0 "              // scale log(2.0), log(0.5), log(1.0)
        "1.0 0.0 0.0 0.0\n",                   // rot w=1, x=0, y=0, z=0
        f);
    std::fclose(f);

    auto c = bromesh::loadSplatPLY(path);
    ASSERT(c.count() == 1, "splat_deg1: 1 splat");
    ASSERT(c.shDegree == 1, "splat_deg1: shDegree is 1");
    ASSERT(c.shStride() == 12, "splat_deg1: shStride is 12");

    // Position verification
    ASSERT(nearf(c.positions[0], 1.0f, 1e-4f), "splat_deg1: pos x");
    ASSERT(nearf(c.positions[1], 2.0f, 1e-4f), "splat_deg1: pos y");
    ASSERT(nearf(c.positions[2], 3.0f, 1e-4f), "splat_deg1: pos z");

    // Scale verification: exp(ln(2.0))=2.0, exp(ln(0.5))=0.5, exp(0)=1.0
    ASSERT(nearf(c.scales[0], 2.0f, 1e-3f), "splat_deg1: scale_0 activated to 2.0");
    ASSERT(nearf(c.scales[1], 0.5f, 1e-3f), "splat_deg1: scale_1 activated to 0.5");
    ASSERT(nearf(c.scales[2], 1.0f, 1e-3f), "splat_deg1: scale_2 activated to 1.0");

    // Opacity verification: sigmoid(1.386294) = 0.80
    ASSERT(nearf(c.opacities[0], 0.80f, 1e-3f), "splat_deg1: opacity sigmoid(logit(0.8)) = 0.8");

    // Rotation verification: w-first in PLY (rot_0=w=1, rot_1=x=0, rot_2=y=0, rot_3=z=0)
    // -> xyzw in GaussianSplatCloud
    ASSERT(nearf(c.rotations[0], 0.0f, 1e-4f), "splat_deg1: rot x=0");
    ASSERT(nearf(c.rotations[1], 0.0f, 1e-4f), "splat_deg1: rot y=0");
    ASSERT(nearf(c.rotations[2], 0.0f, 1e-4f), "splat_deg1: rot z=0");
    ASSERT(nearf(c.rotations[3], 1.0f, 1e-4f), "splat_deg1: rot w=1");

    // SH DC verification
    ASSERT(nearf(c.sh[0], 0.1f, 1e-4f), "splat_deg1: sh DC R = 0.1");
    ASSERT(nearf(c.sh[1], 0.2f, 1e-4f), "splat_deg1: sh DC G = 0.2");
    ASSERT(nearf(c.sh[2], 0.3f, 1e-4f), "splat_deg1: sh DC B = 0.3");

    // SH rest channel-major to interleaved transpose verification:
    // INRIA channel-major:
    //   f_rest_0..2 = (R1, R2, R3) = (1, 2, 3)
    //   f_rest_3..5 = (G1, G2, G3) = (4, 5, 6)
    //   f_rest_6..8 = (B1, B2, B3) = (7, 8, 9)
    // Interleaved layout c.sh[3..11]: (R1, G1, B1, R2, G2, B2, R3, G3, B3)
    ASSERT(nearf(c.sh[3], 1.0f, 1e-4f), "splat_deg1: R1 = 1.0");
    ASSERT(nearf(c.sh[4], 4.0f, 1e-4f), "splat_deg1: G1 = 4.0");
    ASSERT(nearf(c.sh[5], 7.0f, 1e-4f), "splat_deg1: B1 = 7.0");
    ASSERT(nearf(c.sh[6], 2.0f, 1e-4f), "splat_deg1: R2 = 2.0");
    ASSERT(nearf(c.sh[7], 5.0f, 1e-4f), "splat_deg1: G2 = 5.0");
    ASSERT(nearf(c.sh[8], 8.0f, 1e-4f), "splat_deg1: B2 = 8.0");
    ASSERT(nearf(c.sh[9], 3.0f, 1e-4f), "splat_deg1: R3 = 3.0");
    ASSERT(nearf(c.sh[10], 6.0f, 1e-4f), "splat_deg1: G3 = 6.0");
    ASSERT(nearf(c.sh[11], 9.0f, 1e-4f), "splat_deg1: B3 = 9.0");

    // Binary roundtrip save & reload
    const char* binPath = "test_splat_deg1_roundtrip.ply";
    ASSERT(bromesh::saveSplatPLY(c, binPath), "saveSplatPLY binary");
    auto rt = bromesh::loadSplatPLY(binPath);
    ASSERT(rt.count() == 1, "rt count");
    ASSERT(rt.shDegree == 1, "rt shDegree");
    ASSERT(rt.shStride() == 12, "rt shStride");
    for (int a = 0; a < 3; ++a) {
        ASSERT(nearf(rt.positions[a], c.positions[a], 1e-4f), "rt pos");
        ASSERT(nearf(rt.scales[a], c.scales[a], 1e-4f), "rt scale");
    }
    ASSERT(nearf(rt.opacities[0], c.opacities[0], 1e-4f), "rt opacity");
    for (int a = 0; a < 4; ++a) {
        ASSERT(nearf(rt.rotations[a], c.rotations[a], 1e-4f), "rt rot");
    }
    for (int k = 0; k < 12; ++k) {
        ASSERT(nearf(rt.sh[k], c.sh[k], 1e-4f), "rt sh");
    }

    std::remove(path);
    std::remove(binPath);
}

TEST(splat_transforms) {
    auto cloud = makeCloud(0, 5);
    auto origPositions = cloud.positions;
    auto origScales = cloud.scales;
    auto b0 = cloud.bounds();

    // 1. translateSplats
    bromesh::translateSplats(cloud, 10.0f, -5.0f, 2.0f);
    auto b1 = cloud.bounds();
    ASSERT(nearf(b1.min.x, b0.min.x + 10.0f), "splat translate bounds min.x");
    ASSERT(nearf(b1.max.x, b0.max.x + 10.0f), "splat translate bounds max.x");
    ASSERT(nearf(b1.min.y, b0.min.y - 5.0f), "splat translate bounds min.y");
    ASSERT(nearf(b1.max.y, b0.max.y - 5.0f), "splat translate bounds max.y");
    ASSERT(nearf(b1.min.z, b0.min.z + 2.0f), "splat translate bounds min.z");
    ASSERT(nearf(b1.max.z, b0.max.z + 2.0f), "splat translate bounds max.z");

    for (size_t i = 0; i < cloud.count(); ++i) {
        ASSERT(nearf(cloud.positions[i * 3 + 0], origPositions[i * 3 + 0] + 10.0f), "splat trans x");
        ASSERT(nearf(cloud.positions[i * 3 + 1], origPositions[i * 3 + 1] - 5.0f), "splat trans y");
        ASSERT(nearf(cloud.positions[i * 3 + 2], origPositions[i * 3 + 2] + 2.0f), "splat trans z");
        ASSERT(nearf(cloud.scales[i * 3 + 0], origScales[i * 3 + 0]), "splat trans scales unmod 0");
        ASSERT(nearf(cloud.scales[i * 3 + 1], origScales[i * 3 + 1]), "splat trans scales unmod 1");
        ASSERT(nearf(cloud.scales[i * 3 + 2], origScales[i * 3 + 2]), "splat trans scales unmod 2");
    }

    // 2. scaleSplats
    bromesh::scaleSplats(cloud, 2.0f, 3.0f, 0.5f);
    for (size_t i = 0; i < cloud.count(); ++i) {
        float expectedX = (origPositions[i * 3 + 0] + 10.0f) * 2.0f;
        float expectedY = (origPositions[i * 3 + 1] - 5.0f) * 3.0f;
        float expectedZ = (origPositions[i * 3 + 2] + 2.0f) * 0.5f;
        ASSERT(nearf(cloud.positions[i * 3 + 0], expectedX), "splat scale pos x");
        ASSERT(nearf(cloud.positions[i * 3 + 1], expectedY), "splat scale pos y");
        ASSERT(nearf(cloud.positions[i * 3 + 2], expectedZ), "splat scale pos z");

        ASSERT(nearf(cloud.scales[i * 3 + 0], origScales[i * 3 + 0] * 2.0f), "splat scale stream 0");
        ASSERT(nearf(cloud.scales[i * 3 + 1], origScales[i * 3 + 1] * 3.0f), "splat scale stream 1");
        ASSERT(nearf(cloud.scales[i * 3 + 2], origScales[i * 3 + 2] * 0.5f), "splat scale stream 2");
    }

    // 3. transformSplats with 4x4 matrix (rot 90 deg about Z, scale 2x, translation 1,2,3)
    auto c2 = makeCloud(0, 3);
    auto c2Pos = c2.positions;
    auto c2Scale = c2.scales;
    float m[16] = {
        0.0f,  2.0f, 0.0f, 0.0f,
       -2.0f,  0.0f, 0.0f, 0.0f,
        0.0f,  0.0f, 2.0f, 0.0f,
        1.0f,  2.0f, 3.0f, 1.0f
    };
    bromesh::transformSplats(c2, m);
    ASSERT(c2.validate(), "transformed splats validate");
    for (size_t i = 0; i < c2.count(); ++i) {
        float x = c2Pos[i * 3 + 0];
        float y = c2Pos[i * 3 + 1];
        float z = c2Pos[i * 3 + 2];
        float tx = -2.0f * y + 1.0f;
        float ty =  2.0f * x + 2.0f;
        float tz =  2.0f * z + 3.0f;
        ASSERT(nearf(c2.positions[i * 3 + 0], tx), "affine pos x");
        ASSERT(nearf(c2.positions[i * 3 + 1], ty), "affine pos y");
        ASSERT(nearf(c2.positions[i * 3 + 2], tz), "affine pos z");

        ASSERT(nearf(c2.scales[i * 3 + 0], c2Scale[i * 3 + 0] * 2.0f), "affine scale 0");
        ASSERT(nearf(c2.scales[i * 3 + 1], c2Scale[i * 3 + 1] * 2.0f), "affine scale 1");
        ASSERT(nearf(c2.scales[i * 3 + 2], c2Scale[i * 3 + 2] * 2.0f), "affine scale 2");

        float qx = c2.rotations[i * 4 + 0];
        float qy = c2.rotations[i * 4 + 1];
        float qz = c2.rotations[i * 4 + 2];
        float qw = c2.rotations[i * 4 + 3];
        ASSERT(nearf(qx * qx + qy * qy + qz * qz + qw * qw, 1.0f), "affine quat unit len");
    }
}

TEST(splat_filter) {
    bromesh::GaussianSplatCloud cloud;
    cloud.shDegree = 0;
    cloud.positions = {
        0.0f, 0.0f, 0.0f,    // 0: low opacity (0.001) -> dropped
        0.5f, 0.5f, 0.5f,    // 1: kept
        5.0f, 5.0f, 5.0f,    // 2: outside cropBox -> dropped
        0.2f, 0.2f, 0.2f,    // 3: large scale (5.0) -> dropped
        0.8f, 0.8f, 0.8f,    // 4: kept
        0.1f, 0.1f, 0.1f     // 5: low opacity (0.002) -> dropped
    };
    cloud.scales = {
        0.1f, 0.1f, 0.1f,
        0.2f, 0.2f, 0.2f,
        0.1f, 0.1f, 0.1f,
        5.0f, 0.1f, 0.1f,
        0.3f, 0.3f, 0.3f,
        0.1f, 0.1f, 0.1f
    };
    cloud.rotations = {
        0, 0, 0, 1,
        0, 0, 0, 1,
        0, 0, 0, 1,
        0, 0, 0, 1,
        0, 0, 0, 1,
        0, 0, 0, 1
    };
    cloud.opacities = {
        0.001f,
        0.5f,
        0.8f,
        0.8f,
        0.9f,
        0.002f
    };
    cloud.sh = {
        1, 1, 1,
        2, 2, 2,
        3, 3, 3,
        4, 4, 4,
        5, 5, 5,
        6, 6, 6
    };

    ASSERT(cloud.validate(), "pre-filter cloud validates");

    bromath::AABB3 cropBox{{-1.0f, -1.0f, -1.0f}, {2.0f, 2.0f, 2.0f}};
    bromesh::SplatFilterOptions opts;
    opts.minOpacity = 0.01f;
    opts.cropBox = &cropBox;
    opts.maxScale = 2.0f;

    bromesh::filterSplats(cloud, opts);

    ASSERT(cloud.validate(), "post-filter cloud validates");
    ASSERT(cloud.count() == 2, "filter: exactly 2 splats survive");

    ASSERT(nearf(cloud.positions[0], 0.5f) && nearf(cloud.positions[1], 0.5f) && nearf(cloud.positions[2], 0.5f), "splat 1 pos");
    ASSERT(nearf(cloud.positions[3], 0.8f) && nearf(cloud.positions[4], 0.8f) && nearf(cloud.positions[5], 0.8f), "splat 4 pos");

    ASSERT(nearf(cloud.opacities[0], 0.5f), "splat 1 opacity");
    ASSERT(nearf(cloud.opacities[1], 0.9f), "splat 4 opacity");

    ASSERT(nearf(cloud.scales[0], 0.2f), "splat 1 scale");
    ASSERT(nearf(cloud.scales[3], 0.3f), "splat 4 scale");

    ASSERT(nearf(cloud.sh[0], 2.0f), "splat 1 sh");
    ASSERT(nearf(cloud.sh[3], 5.0f), "splat 4 sh");

    auto b = cloud.bounds();
    ASSERT(nearf(b.min.x, 0.5f) && nearf(b.max.x, 0.8f), "filtered bounds x");
    ASSERT(nearf(b.min.y, 0.5f) && nearf(b.max.y, 0.8f), "filtered bounds y");
    ASSERT(nearf(b.min.z, 0.5f) && nearf(b.max.z, 0.8f), "filtered bounds z");
}

TEST(splat_merge) {
    auto cloudA = makeCloud(0, 3);
    auto cloudB = makeCloud(2, 4);

    auto merged = bromesh::mergeSplats({cloudA, cloudB});

    ASSERT(merged.validate(), "merged splats validate");
    ASSERT(merged.count() == 7, "merged count is 3+4=7");
    ASSERT(merged.shDegree == 2, "merged degree is max(0, 2) = 2");
    ASSERT(merged.shStride() == 27, "merged stride is 3*(2+1)^2 = 27");
    ASSERT(merged.sh.size() == 7 * 27, "merged sh total size");
    ASSERT(merged.positions.size() == 7 * 3, "merged positions size");
    ASSERT(merged.scales.size() == 7 * 3, "merged scales size");
    ASSERT(merged.rotations.size() == 7 * 4, "merged rotations size");
    ASSERT(merged.opacities.size() == 7, "merged opacities size");

    for (size_t v = 0; v < 3; ++v) {
        for (int c = 0; c < 3; ++c) {
            ASSERT(nearf(merged.positions[v * 3 + c], cloudA.positions[v * 3 + c]), "merged pos A");
            ASSERT(nearf(merged.scales[v * 3 + c], cloudA.scales[v * 3 + c]), "merged scale A");
            ASSERT(nearf(merged.sh[v * 27 + c], cloudA.sh[v * 3 + c]), "merged sh DC A");
        }
        for (int k = 3; k < 27; ++k) {
            ASSERT(nearf(merged.sh[v * 27 + k], 0.0f), "merged higher SH padded 0");
        }
    }

    for (size_t v = 0; v < 4; ++v) {
        size_t mv = 3 + v;
        for (int c = 0; c < 3; ++c) {
            ASSERT(nearf(merged.positions[mv * 3 + c], cloudB.positions[v * 3 + c]), "merged pos B");
            ASSERT(nearf(merged.scales[mv * 3 + c], cloudB.scales[v * 3 + c]), "merged scale B");
        }
        for (int k = 0; k < 27; ++k) {
            ASSERT(nearf(merged.sh[mv * 27 + k], cloudB.sh[v * 27 + k]), "merged sh full B");
        }
    }

    const bromesh::GaussianSplatCloud arr[2] = {cloudA, cloudB};
    auto mergedArr = bromesh::mergeSplats(arr, 2);
    ASSERT(mergedArr.count() == 7, "mergeSplats array overload");
    ASSERT(mergedArr.validate(), "mergeSplats array validate");
}

TEST(mesh_to_splats) {
    auto sphereMesh = bromesh::sphere(1.0f, 32, 24);
    ASSERT(!sphereMesh.empty(), "sphere primitive valid");

    bromesh::MeshToSplatsOptions opts;
    opts.splatCount = 1000;
    opts.opacity = 0.95f;
    opts.seed = 42;

    auto splats = bromesh::meshToSplats(sphereMesh, opts);
    ASSERT(splats.validate(), "meshToSplats cloud validates");
    ASSERT(!splats.empty(), "meshToSplats non-empty");
    ASSERT(splats.count() == 1000, "splat count matches requested");
    ASSERT(splats.shDegree == 0, "meshToSplats degree 0");
    ASSERT(splats.shStride() == 3, "meshToSplats stride 3");
    ASSERT(splats.scales.size() == 3000, "scales size 3000");
    ASSERT(splats.rotations.size() == 4000, "rotations size 4000");
    ASSERT(splats.opacities.size() == 1000, "opacities size 1000");
    ASSERT(splats.sh.size() == 3000, "sh size 3000");

    auto b = splats.bounds();
    ASSERT(b.min.x >= -1.1f && b.min.x <= -0.8f, "sphere splats bound min x");
    ASSERT(b.max.x <= 1.1f && b.max.x >= 0.8f, "sphere splats bound max x");
    ASSERT(b.min.y >= -1.1f && b.min.y <= -0.8f, "sphere splats bound min y");
    ASSERT(b.max.y <= 1.1f && b.max.y >= 0.8f, "sphere splats bound max y");
    ASSERT(b.min.z >= -1.1f && b.min.z <= -0.8f, "sphere splats bound min z");
    ASSERT(b.max.z <= 1.1f && b.max.z >= 0.8f, "sphere splats bound max z");

    for (size_t i = 0; i < splats.count(); ++i) {
        ASSERT(nearf(splats.opacities[i], 0.95f), "splat opacity");
        float r = splats.scales[i * 3 + 0];
        ASSERT(r > 0.0f, "splat radius positive");
        ASSERT(nearf(splats.scales[i * 3 + 1], r), "splat scale y == radius");
        ASSERT(nearf(splats.scales[i * 3 + 2], r * 0.1f), "splat scale z == radius * 0.1");

        float qx = splats.rotations[i * 4 + 0];
        float qy = splats.rotations[i * 4 + 1];
        float qz = splats.rotations[i * 4 + 2];
        float qw = splats.rotations[i * 4 + 3];
        ASSERT(nearf(qx * qx + qy * qy + qz * qz + qw * qw, 1.0f, 1e-3f), "splat unit quat");
    }

    const char* plyPath = "test_mesh_to_splats_roundtrip.ply";
    ASSERT(bromesh::saveSplatPLY(splats, plyPath), "saveSplatPLY for meshToSplats");

    auto loaded = bromesh::loadSplatPLY(plyPath);
    ASSERT(loaded.validate(), "reloaded mesh splats validate");
    ASSERT(loaded.count() == splats.count(), "reloaded count matches");
    ASSERT(loaded.shDegree == 0, "reloaded degree 0");

    for (size_t i = 0; i < splats.count(); ++i) {
        for (int c = 0; c < 3; ++c) {
            ASSERT(nearf(loaded.positions[i * 3 + c], splats.positions[i * 3 + c], 1e-4f), "rt pos");
            ASSERT(nearf(loaded.scales[i * 3 + c], splats.scales[i * 3 + c], 1e-3f), "rt scale");
            ASSERT(nearf(loaded.sh[i * 3 + c], splats.sh[i * 3 + c], 1e-4f), "rt sh");
        }
        ASSERT(nearf(loaded.opacities[i], splats.opacities[i], 1e-3f), "rt opacity");
    }

    std::remove(plyPath);
}

