#include "test_framework.h"

#include "bromesh/procedural/leaf_cluster.h"
#include "bromesh/procedural/leaf_scatter.h"
#include "bromesh/procedural/plants.h"
#include "bromesh/procedural/space_colonization.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace bromesh;
using namespace bromath;

namespace {

// Helper to validate basic integrity of generated mesh
void validateMesh(const MeshData& m, const char* label, bool expectColors = true) {
    ASSERT(!m.empty(), label);
    ASSERT(m.vertexCount() > 0, label);
    ASSERT(m.triangleCount() > 0, label);
    ASSERT(m.hasNormals(), label);
    ASSERT(m.hasUVs(), label);
    if (expectColors) {
        ASSERT(m.hasColors(), label);
    }

    // Verify normals are normalized and finite
    for (size_t v = 0; v < m.vertexCount(); ++v) {
        float nx = m.normals[v * 3 + 0];
        float ny = m.normals[v * 3 + 1];
        float nz = m.normals[v * 3 + 2];
        float len2 = nx * nx + ny * ny + nz * nz;
        ASSERT(std::isfinite(len2) && len2 > 0.5f && len2 < 1.5f, "normal length ~ 1.0");

        // Verify UVs are finite
        float u = m.uvs[v * 2 + 0];
        float v_coord = m.uvs[v * 2 + 1];
        ASSERT(std::isfinite(u) && std::isfinite(v_coord), "uvs are finite");

        // Verify wind bend in R channel is in [0, 1]
        if (m.hasColors()) {
            float r = m.colors[v * 4 + 0];
            float g = m.colors[v * 4 + 1];
            float b = m.colors[v * 4 + 2];
            float a = m.colors[v * 4 + 3];
            ASSERT(r >= -1e-4f && r <= 1.0001f, "wind bend R channel in [0, 1]");
            ASSERT(std::isfinite(g) && std::isfinite(b) && std::isfinite(a), "colors are finite");
        }
    }

    // Verify no degenerate triangles
    for (size_t t = 0; t < m.triangleCount(); ++t) {
        uint32_t a = m.indices[t * 3 + 0];
        uint32_t b = m.indices[t * 3 + 1];
        uint32_t c = m.indices[t * 3 + 2];
        ASSERT(a != b && b != c && a != c, "no degenerate triangle indices");
        ASSERT(a < m.vertexCount() && b < m.vertexCount() && c < m.vertexCount(), "indices in bounds");
    }
}

} // namespace

// =============================================================================
// Test individual Phyllotaxy patterns
// =============================================================================

TEST(leaf_cluster_alternate) {
    LeafClusterOptions opts;
    opts.count = 6;
    opts.twigLength = 0.25f;
    opts.twigRadius = 0.005f;
    opts.petioleLength = 0.04f;
    opts.includeTwigMesh = true;

    MeshData m = leafCluster(Phyllotaxy::Alternate, opts);
    validateMesh(m, "leafCluster Alternate with twig");
    ASSERT(m.vertexCount() > 100, "expected substantial vertex count");

    // Test with includeTwigMesh = false
    opts.includeTwigMesh = false;
    MeshData mNoTwig = leafCluster(Phyllotaxy::Alternate, opts);
    validateMesh(mNoTwig, "leafCluster Alternate without twig");
    ASSERT(mNoTwig.vertexCount() < m.vertexCount(), "no-twig has fewer vertices");
}

TEST(leaf_cluster_opposite) {
    LeafClusterOptions opts;
    opts.count = 6; // 3 pairs of decussate leaves
    opts.twigLength = 0.30f;
    opts.petioleLength = 0.03f;
    opts.includeTwigMesh = true;

    MeshData m = leafCluster(Phyllotaxy::Opposite, opts);
    validateMesh(m, "leafCluster Opposite (even count)");

    // Test odd count
    opts.count = 5;
    MeshData mOdd = leafCluster(Phyllotaxy::Opposite, opts);
    validateMesh(mOdd, "leafCluster Opposite (odd count)");
}

TEST(leaf_cluster_spiral) {
    LeafClusterOptions opts;
    opts.count = 7;
    opts.twigLength = 0.35f;
    opts.shape = LeafShape::Pointed;
    opts.shapedSilhouette = true;

    MeshData m = leafCluster(Phyllotaxy::Spiral, opts);
    validateMesh(m, "leafCluster Spiral");
    ASSERT(m.triangleCount() > 50, "triangle count > 50");
}

TEST(leaf_cluster_fascicle) {
    LeafClusterOptions opts;
    opts.count = 5; // 5-needle pine fascicle
    opts.shape = LeafShape::Needle;
    opts.leafLength = 0.25f;
    opts.leafWidth = 0.015f;
    opts.petioleLength = 0.005f;
    opts.spread = 0.35f;

    MeshData m = leafCluster(Phyllotaxy::Fascicle, opts);
    validateMesh(m, "leafCluster Fascicle 5 needles");
}

TEST(leaf_cluster_compound_pinnate) {
    LeafClusterOptions opts;
    opts.count = 9; // 1 terminal + 4 pairs
    opts.twigLength = 0.40f;
    opts.petioleLength = 0.02f;
    opts.shape = LeafShape::Oval;

    MeshData m = leafCluster(Phyllotaxy::CompoundPinnate, opts);
    validateMesh(m, "leafCluster CompoundPinnate (odd count)");

    // Test even count
    opts.count = 6;
    MeshData mEven = leafCluster(Phyllotaxy::CompoundPinnate, opts);
    validateMesh(mEven, "leafCluster CompoundPinnate (even count)");
}

// =============================================================================
// Edge Cases & Options Coverage
// =============================================================================

TEST(leaf_cluster_edge_cases) {
    LeafClusterOptions opts;

    // Zero count with twig mesh: should only have twig geometry
    opts.count = 0;
    opts.includeTwigMesh = true;
    MeshData mTwigOnly = leafCluster(Phyllotaxy::Alternate, opts);
    ASSERT(!mTwigOnly.empty(), "count 0 with includeTwigMesh produces twig");
    validateMesh(mTwigOnly, "twig only");

    // Zero count without twig: empty
    opts.includeTwigMesh = false;
    MeshData mEmpty = leafCluster(Phyllotaxy::Alternate, opts);
    ASSERT(mEmpty.empty(), "count 0 without twig is empty");

    // Single leaf (count = 1) for all phyllotaxies
    opts.count = 1;
    opts.includeTwigMesh = true;
    for (int p = 0; p <= 4; ++p) {
        MeshData m1 = leafCluster(static_cast<Phyllotaxy>(p), opts);
        validateMesh(m1, "single leaf cluster");
    }

    // Zero petiole length
    opts.count = 4;
    opts.petioleLength = 0.0f;
    MeshData mNoPetiole = leafCluster(Phyllotaxy::Spiral, opts);
    validateMesh(mNoPetiole, "zero petiole length");

    // fullUV mode
    opts.fullUV = true;
    opts.shapedSilhouette = false;
    MeshData mFullUV = leafCluster(Phyllotaxy::Alternate, opts);
    validateMesh(mFullUV, "fullUV mode");
}

TEST(leaf_cluster_shapes_and_deformations) {
    LeafClusterOptions opts;
    opts.count = 4;
    opts.leafBend = 0.5f;
    opts.leafCurl = 0.3f;
    opts.leafCup = 0.4f;
    opts.droop = 0.5f;
    opts.upBias = 0.8f;
    opts.spread = 0.9f;

    LeafShape shapes[] = {
        LeafShape::Oval,
        LeafShape::Pointed,
        LeafShape::Lobed,
        LeafShape::Needle,
        LeafShape::Frond,
        LeafShape::Petal,
    };

    for (LeafShape shape : shapes) {
        opts.shape = shape;
        opts.shapedSilhouette = true;
        MeshData m = leafCluster(Phyllotaxy::Spiral, opts);
        validateMesh(m, "leafCluster deformed shape silhouette");
    }
}

// =============================================================================
// Scattering on Branch Skeletons
// =============================================================================

TEST(leaf_cluster_scatter_on_branches) {
    // Construct synthetic branch skeleton: 1 trunk + 2 child branches
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

    LeafPlacementOptions placeOpts;
    placeOpts.minDepth = 1; // only on branches, not trunk
    placeOpts.maxRadius = 0.05f;
    placeOpts.perUnitLength = 15.0f;
    placeOpts.seed = 42;

    LeafPlacements placements = placeLeafClustersOnBranches(segs, placeOpts);
    ASSERT(placements.count() > 0, "placed clusters on branches");
    ASSERT(placements.transforms.size() == placements.count() * 16, "transform buffer 16 floats per instance");

    // Verify transforms: rows 0..2 are affine, row 3 is tint (1,1,1,1)
    for (size_t i = 0; i < placements.count(); ++i) {
        const float* M = &placements.transforms[i * 16];
        ASSERT(M[12] == 1.0f && M[13] == 1.0f && M[14] == 1.0f && M[15] == 1.0f, "valid tint color");
        // Translation point (px, py, pz)
        float px = M[3], py = M[7], pz = M[11];
        ASSERT(py >= 1.9f && py <= 3.2f, "cluster origin within branch bounds");
        ASSERT(std::isfinite(px) && std::isfinite(py) && std::isfinite(pz), "finite origin");
    }

    // Scatter full cluster mesh across all phyllotaxies
    LeafClusterOptions clusterOpts;
    clusterOpts.count = 5;
    clusterOpts.twigLength = 0.15f;

    Phyllotaxy phyllotaxies[] = {
        Phyllotaxy::Alternate,
        Phyllotaxy::Opposite,
        Phyllotaxy::Spiral,
        Phyllotaxy::Fascicle,
        Phyllotaxy::CompoundPinnate,
    };

    for (Phyllotaxy p : phyllotaxies) {
        MeshData foliage = scatterLeafClusters(segs, p, clusterOpts, placeOpts);
        validateMesh(foliage, "scatterLeafClusters result");
        ASSERT(foliage.vertexCount() >= placements.count() * 50, "adequate vertex count from merged instances");
    }
}

TEST(leaf_cluster_scatter_with_obstacles_and_weights) {
    std::vector<BranchSegment> segs;
    BranchSegment s1;
    s1.from = {0.0f, 0.0f, 0.0f};
    s1.to = {0.0f, 3.0f, 0.0f};
    s1.radius = 0.02f;
    s1.depth = 1;
    s1.parent = -1;
    segs.push_back(s1);

    // Obstacle capsule blocking the middle of the branch
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

    LeafPlacements pl = placeLeafClustersOnBranches(segs, opts);
    ASSERT(pl.count() > 0, "placed clusters avoiding obstacle");

    // Verify none of the origins lie inside the obstacle capsule
    for (size_t i = 0; i < pl.count(); ++i) {
        const float* M = &pl.transforms[i * 16];
        Vec3 pos{M[3], M[7], M[11]};
        ASSERT(!field.tooClose(pos, opts.obstacleClearance, -1), "cluster origin clears obstacle");
    }

    // Density weight: weight 0 on segment produces 0 placements
    opts.densityWeight = {0.0f};
    LeafPlacements plZero = placeLeafClustersOnBranches(segs, opts);
    ASSERT(plZero.count() == 0, "densityWeight 0 produces zero placements");

    // Empty segments or empty mesh returns empty MeshData
    std::vector<BranchSegment> emptySegs;
    MeshData emptyScatter = scatterLeafClusters(emptySegs, Phyllotaxy::Alternate);
    ASSERT(emptyScatter.empty(), "scattering on empty segments is empty");
}
