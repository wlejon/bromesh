#include "test_framework.h"

#include "bromesh/procedural/leaf_scatter.h"
#include "bromesh/procedural/obstacle_field.h"
#include "bromesh/procedural/plants.h"
#include "bromesh/procedural/space_colonization.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace bromesh;
using namespace bromath;

namespace {

void validateMesh(const MeshData& m, const char* label, bool expectColors = true) {
    ASSERT(!m.empty(), label);
    ASSERT(m.vertexCount() > 0, label);
    ASSERT(m.triangleCount() > 0, label);
    ASSERT(m.hasNormals(), label);
    ASSERT(m.hasUVs(), label);
    if (expectColors) {
        ASSERT(m.hasColors(), label);
    }

    for (size_t v = 0; v < m.vertexCount(); ++v) {
        float nx = m.normals[v * 3 + 0];
        float ny = m.normals[v * 3 + 1];
        float nz = m.normals[v * 3 + 2];
        float len2 = nx * nx + ny * ny + nz * nz;
        ASSERT(std::isfinite(len2) && len2 > 0.5f && len2 < 1.5f, "normal length ~ 1.0");

        float u = m.uvs[v * 2 + 0];
        float v_coord = m.uvs[v * 2 + 1];
        ASSERT(std::isfinite(u) && std::isfinite(v_coord), "uvs are finite");

        if (m.hasColors()) {
            float r = m.colors[v * 4 + 0];
            float g = m.colors[v * 4 + 1];
            float b = m.colors[v * 4 + 2];
            float a = m.colors[v * 4 + 3];
            ASSERT(r >= -1e-4f && r <= 1.0001f, "wind bend R channel in [0, 1]");
            ASSERT(std::isfinite(g) && std::isfinite(b) && std::isfinite(a), "colors are finite");
        }
    }

    for (size_t t = 0; t < m.triangleCount(); ++t) {
        uint32_t a = m.indices[t * 3 + 0];
        uint32_t b = m.indices[t * 3 + 1];
        uint32_t c = m.indices[t * 3 + 2];
        ASSERT(a != b && b != c && a != c, "no degenerate triangle indices");
        ASSERT(a < m.vertexCount() && b < m.vertexCount() && c < m.vertexCount(), "indices in bounds");
    }
}

std::vector<BranchSegment> createBranchSkeleton() {
    std::vector<BranchSegment> segs;
    BranchSegment trunk;
    trunk.from = {0.0f, 0.0f, 0.0f};
    trunk.to = {0.0f, 2.0f, 0.0f};
    trunk.radius = 0.08f;
    trunk.depth = 0;
    trunk.parent = -1;
    segs.push_back(trunk);

    BranchSegment b1;
    b1.from = {0.0f, 2.0f, 0.0f};
    b1.to = {1.5f, 3.0f, 0.0f};
    b1.radius = 0.03f;
    b1.depth = 1;
    b1.parent = 0;
    segs.push_back(b1);

    BranchSegment b2;
    b2.from = {0.0f, 2.0f, 0.0f};
    b2.to = {-1.5f, 3.0f, 0.0f};
    b2.radius = 0.03f;
    b2.depth = 1;
    b2.parent = 0;
    segs.push_back(b2);

    return segs;
}

} // namespace

TEST(leaf_scatter_placements_layout) {
    auto segs = createBranchSkeleton();
    LeafPlacementOptions opts;
    opts.minDepth = 1;
    opts.maxRadius = 0.05f;
    opts.perUnitLength = 15.0f;
    opts.seed = 42;

    LeafPlacements placements = placeLeavesOnBranches(segs, opts);
    ASSERT(placements.count() > 0, "placed leaves on branches");
    ASSERT(placements.transforms.size() == placements.count() * 16, "transform buffer 16 floats per instance");
    ASSERT(placements.branchRadius.size() == placements.count(), "branchRadius matches count");
    ASSERT(placements.branchDepth.size() == placements.count(), "branchDepth matches count");

    for (size_t i = 0; i < placements.count(); ++i) {
        const float* M = &placements.transforms[i * 16];
        ASSERT(M[12] == 1.0f && M[13] == 1.0f && M[14] == 1.0f && M[15] == 1.0f, "valid tint color");
        float px = M[3], py = M[7], pz = M[11];
        ASSERT(py >= 1.9f && py <= 3.2f, "leaf origin within branch bounds");
        ASSERT(std::isfinite(px) && std::isfinite(py) && std::isfinite(pz), "finite origin");
        ASSERT(placements.branchRadius[i] <= 0.05f + 1e-6f, "radius passes filter");
        ASSERT(placements.branchDepth[i] >= 1, "depth passes filter");
    }
}

// perUnitLength is caller data (a JS option in the bindings): a NaN or
// negative density places nothing, and a huge or infinite one is capped per
// segment instead of going through an int conversion it cannot survive.
TEST(leaf_scatter_density_is_clamped) {
    std::vector<BranchSegment> segs(1);
    segs[0].from = {0.0f, 1.0f, 0.0f};
    segs[0].to = {0.0f, 2.0f, 0.0f};
    segs[0].radius = 0.01f;
    segs[0].depth = 1;
    segs[0].parent = -1;

    LeafPlacementOptions opts;
    for (float bad : {std::nanf(""), -5.0f, 0.0f}) {
        opts.perUnitLength = bad;
        ASSERT(placeLeavesOnBranches(segs, opts).count() == 0, "no leaves for a NaN / non-positive density");
    }
    for (float huge : {1e30f, std::numeric_limits<float>::infinity()}) {
        opts.perUnitLength = huge;
        LeafPlacements pl = placeLeavesOnBranches(segs, opts);
        ASSERT(pl.count() == static_cast<size_t>(kMaxLeavesPerSegment), "a huge density is capped per segment");
    }
}

// The instance layout the header documents: row-major 3x4 affine plus a
// white tint, the translation at 3 / 7 / 11 -- a leaf's origin lies on its
// segment, which a column-major reading (translation at 12..14 = 1, 1, 1)
// would not give.
TEST(leaf_scatter_instance_layout_is_row_major) {
    std::vector<BranchSegment> segs(1);
    segs[0].from = {4.0f, 1.0f, -3.0f};
    segs[0].to = {4.0f, 2.0f, -3.0f};
    segs[0].radius = 0.01f;
    segs[0].depth = 1;
    segs[0].parent = -1;
    LeafPlacementOptions opts;
    opts.perUnitLength = 10.0f;
    LeafPlacements pl = placeLeavesOnBranches(segs, opts);
    ASSERT(pl.count() > 0, "placed leaves");
    for (size_t i = 0; i < pl.count(); ++i) {
        const float* M = &pl.transforms[i * 16];
        ASSERT(std::fabs(M[3] - 4.0f) < 1e-4f && std::fabs(M[11] + 3.0f) < 1e-4f, "translation at 3 / 11");
        ASSERT(M[7] >= 1.0f - 1e-4f && M[7] <= 2.0f + 1e-4f, "translation at 7 on the segment");
        ASSERT(M[12] == 1.0f && M[13] == 1.0f && M[14] == 1.0f && M[15] == 1.0f, "white tint at 12..15");
        // The basis columns (local X, Y, Z) are orthogonal and equally scaled.
        const float cx[3] = {M[0], M[4], M[8]}, cy[3] = {M[1], M[5], M[9]}, cz[3] = {M[2], M[6], M[10]};
        auto dot = [](const float* a, const float* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
        const float s2 = dot(cx, cx);
        ASSERT(s2 > 0.0f, "non-degenerate basis");
        ASSERT(std::fabs(dot(cy, cy) - s2) < 1e-4f * s2 && std::fabs(dot(cz, cz) - s2) < 1e-4f * s2, "uniform scale");
        ASSERT(std::fabs(dot(cx, cy)) < 1e-4f * s2 && std::fabs(dot(cx, cz)) < 1e-4f * s2 &&
               std::fabs(dot(cy, cz)) < 1e-4f * s2, "orthogonal basis");
    }
}

TEST(leaf_scatter_stamp_mesh) {
    auto segs = createBranchSkeleton();
    MeshData leaf = leafCard(LeafShape::Pointed);

    LeafPlacementOptions opts;
    opts.minDepth = 1;
    opts.maxRadius = 0.05f;
    opts.perUnitLength = 10.0f;
    opts.seed = 99;

    MeshData scattered = scatterLeaves(segs, leaf, opts);
    validateMesh(scattered, "scatterLeaves result");
}

TEST(leaf_scatter_obstacles_and_density_weights) {
    std::vector<BranchSegment> segs;
    BranchSegment s1;
    s1.from = {0.0f, 0.0f, 0.0f};
    s1.to = {0.0f, 3.0f, 0.0f};
    s1.radius = 0.02f;
    s1.depth = 1;
    s1.parent = -1;
    segs.push_back(s1);

    Capsule obs;
    obs.a = {-1.0f, 1.5f, 0.0f};
    obs.b = {1.0f, 1.5f, 0.0f};
    obs.radius = 0.5f;
    obs.tag = -1;
    CapsuleField field({obs});

    LeafPlacementOptions opts;
    opts.minDepth = 1;
    opts.perUnitLength = 30.0f;
    opts.avoid = &field;
    opts.obstacleClearance = 0.1f;
    opts.dedupRadius = 0.05f;

    LeafPlacements pl = placeLeavesOnBranches(segs, opts);
    ASSERT(pl.count() > 0, "placed leaves avoiding obstacle");

    for (size_t i = 0; i < pl.count(); ++i) {
        const float* M = &pl.transforms[i * 16];
        Vec3 pos{M[3], M[7], M[11]};
        ASSERT(!field.tooClose(pos, opts.obstacleClearance, -1), "leaf origin clears obstacle");
    }

    opts.densityWeight = {0.0f};
    LeafPlacements plZero = placeLeavesOnBranches(segs, opts);
    ASSERT(plZero.count() == 0, "densityWeight 0 produces zero placements");

    std::vector<BranchSegment> emptySegs;
    MeshData emptyScatter = scatterLeaves(emptySegs, leafCard(LeafShape::Oval));
    ASSERT(emptyScatter.empty(), "scattering on empty segments is empty");

    MeshData emptySrcScatter = scatterLeaves(segs, MeshData{});
    ASSERT(emptySrcScatter.empty(), "scattering empty source mesh is empty");
}

TEST(leaf_scatter_pack_anchors_functionality) {
    std::vector<Vec3> cand;
    for (int i = 0; i < 30; ++i) {
        cand.push_back({static_cast<float>(i) * 0.1f, 0.0f, 0.0f});
    }

    AnchorPackOptions opts;
    opts.minSpacing = 0.25f;
    opts.maxCount = 4;
    opts.seed = 1234;

    auto accepted = packAnchors(cand, nullptr, {}, opts);
    ASSERT(accepted.size() == 4, "maxCount respected");

    for (size_t i = 0; i < accepted.size(); ++i) {
        for (size_t j = i + 1; j < accepted.size(); ++j) {
            float dist = vdist(cand[accepted[i]], cand[accepted[j]]);
            ASSERT(dist >= opts.minSpacing - 1e-5f, "minSpacing respected");
        }
    }
}
