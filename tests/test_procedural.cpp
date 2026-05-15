#include "test_framework.h"

#include "bromesh/optimization/spatial_hash.h"
#include "bromesh/manipulation/sweep.h"
#include "bromesh/manipulation/bezier_sweep.h"
#include "bromesh/procedural/lsystem.h"
#include "bromesh/procedural/lsystem_turtle.h"
#include "bromesh/procedural/space_colonization.h"
#include "bromesh/procedural/branches.h"
#include "bromesh/procedural/leaf_scatter.h"
#include "bromesh/procedural/obstacle_field.h"
#include "bromesh/procedural/plants.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <set>

using namespace bromesh;
using namespace bromath;

TEST(spatial_hash_radius_query) {
    SpatialHash3D hash(1.0f);
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> uni(-10.0f, 10.0f);
    std::vector<Vec3> pts;
    pts.reserve(1000);
    for (int i = 0; i < 1000; ++i) {
        Vec3 p{uni(rng), uni(rng), uni(rng)};
        pts.push_back(p);
        hash.insert(p, i);
    }
    ASSERT(hash.size() == 1000, "size after insert");

    Vec3 center{0, 0, 0};
    float radius = 3.0f;
    std::set<int32_t> brute;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (vdist(pts[i], center) <= radius) brute.insert(static_cast<int32_t>(i));
    }
    std::vector<int32_t> got;
    hash.radiusQuery(center, radius, got);
    std::set<int32_t> gotSet(got.begin(), got.end());
    ASSERT(brute == gotSet, "radius query matches brute force");

    // Brute-force nearest.
    int32_t bestId = -1;
    float bestD2 = 1e30f;
    for (size_t i = 0; i < pts.size(); ++i) {
        float d2 = vdist2(pts[i], center);
        if (d2 < bestD2) { bestD2 = d2; bestId = static_cast<int32_t>(i); }
    }
    int32_t got2 = hash.nearest(center, 100.0f);
    ASSERT(got2 == bestId, "nearest matches brute force");
}

TEST(spatial_hash_remove_clear) {
    SpatialHash3D h(0.5f);
    h.insert({0,0,0}, 1);
    h.insert({0.1f, 0, 0}, 2);
    h.insert({5, 5, 5}, 3);
    ASSERT(h.size() == 3, "three inserted");
    h.remove(2);
    ASSERT(h.size() == 2, "size after remove");
    std::vector<int32_t> ids;
    h.radiusQuery({0,0,0}, 1.0f, ids);
    ASSERT(ids.size() == 1 && ids[0] == 1, "removed id is gone");
    h.clear();
    ASSERT(h.size() == 0, "size after clear");
}

TEST(sweep_straight_triangle_profile) {
    std::vector<Vec2> profile = {{1, 0}, {-0.5f, 0.866f}, {-0.5f, -0.866f}};
    std::vector<Vec3> path = {{0,0,0}, {0,1,0}, {0,2,0}, {0,3,0}};
    MeshData m = sweep(profile, path);
    ASSERT(!m.empty(), "sweep produced output");
    // 3 verts/ring * 4 rings + 2 cap centroids
    ASSERT(m.vertexCount() == 3 * 4 + 2, "expected vertex count");
    ASSERT(m.hasNormals(), "has normals");
    // Sides: 3 quads * 3 strips = 9 quads = 18 tris; caps: 3 + 3 = 6 tris.
    ASSERT(m.triangleCount() == 18 + 6, "expected triangle count");
    // No degenerate triangles.
    bool degen = false;
    for (size_t t = 0; t < m.triangleCount(); ++t) {
        uint32_t a = m.indices[t*3], b = m.indices[t*3+1], c = m.indices[t*3+2];
        if (a == b || b == c || a == c) { degen = true; break; }
    }
    ASSERT(!degen, "no degenerate triangles");
}

TEST(sweep_circular_path_no_frame_flip) {
    // Quarter circle in XZ plane.
    std::vector<Vec3> path;
    const int N = 24;
    for (int i = 0; i <= N; ++i) {
        float a = 1.5707963f * static_cast<float>(i) / static_cast<float>(N);
        path.push_back({ std::cos(a), 0.0f, std::sin(a) });
    }
    std::vector<Vec2> profile = {{0.05f, 0}, {0, 0.05f}, {-0.05f, 0}, {0, -0.05f}};
    SweepOptions opts;
    opts.capStart = false;
    opts.capEnd = false;
    MeshData m = sweep(profile, path, opts);
    ASSERT(!m.empty(), "sweep produced output");
    // Check that consecutive ring vertex offsets from path points rotate
    // smoothly (no sudden flips). The "ring centroid -> first vertex" offset
    // should make a small angle (< 60 degrees) with the previous one.
    bool flipped = false;
    Vec3 prevOffset{0,0,0};
    for (int i = 0; i <= N; ++i) {
        Vec3 ringFirst{ m.positions[i*4*3 + 0], m.positions[i*4*3 + 1], m.positions[i*4*3 + 2] };
        Vec3 off = ringFirst - path[i];
        off = vnorm(off);
        if (i > 0) {
            float c = vdot(off, prevOffset);
            if (c < 0.5f) { flipped = true; break; }
        }
        prevOffset = off;
    }
    ASSERT(!flipped, "frames do not flip along circular path");
}

TEST(lsystem_algae_fibonacci) {
    LSystem ls;
    ls.setAxiom({{'A', {}}});
    ProductionRule a;
    a.predecessor = 'A';
    a.successor = [](const std::vector<float>&) {
        return std::vector<Module>{ {'A', {}}, {'B', {}} };
    };
    ProductionRule b;
    b.predecessor = 'B';
    b.successor = [](const std::vector<float>&) {
        return std::vector<Module>{ {'A', {}} };
    };
    ls.addRule(a);
    ls.addRule(b);
    // Lengths: 1, 2, 3, 5, 8, 13 ...
    auto r5 = ls.derive(5);
    ASSERT(r5.size() == 13, "fibonacci length after 5 iterations");
    auto r6 = ls.derive(6);
    ASSERT(r6.size() == 21, "fibonacci length after 6 iterations");
}

TEST(lsystem_stochastic_determinism) {
    LSystem ls;
    ls.setAxiom({{'X', {}}});
    ProductionRule r1; r1.predecessor='X'; r1.weight=1.0f;
    r1.successor = [](const std::vector<float>&){
        return std::vector<Module>{ {'X', {}}, {'A', {}} };
    };
    ProductionRule r2; r2.predecessor='X'; r2.weight=1.0f;
    r2.successor = [](const std::vector<float>&){
        return std::vector<Module>{ {'X', {}}, {'B', {}} };
    };
    ls.addRule(r1);
    ls.addRule(r2);
    auto a = ls.derive(8, 12345);
    auto b = ls.derive(8, 12345);
    ASSERT(a.size() == b.size(), "same seed: same length");
    bool same = true;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].symbol != b[i].symbol) { same = false; break; }
    }
    ASSERT(same, "same seed: same output");
    auto c = ls.derive(8, 999);
    bool different = (c.size() != a.size());
    if (!different) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (c[i].symbol != a[i].symbol) { different = true; break; }
        }
    }
    ASSERT(different, "different seeds eventually diverge");
}

TEST(lsystem_parse_modules) {
    auto m = parseModules("F(1.0)[+(25)F]F");
    ASSERT(m.size() == 6, "parsed 6 modules");
    ASSERT(m[0].symbol == 'F' && m[0].params.size() == 1, "F(1.0)");
    ASSERT(std::fabs(m[0].params[0] - 1.0f) < 1e-5f, "F param");
    ASSERT(m[1].symbol == '[', "bracket");
    ASSERT(m[2].symbol == '+' && m[2].params.size() == 1 && m[2].params[0] == 25.0f, "+(25)");
    ASSERT(m[3].symbol == 'F', "inner F");
    ASSERT(m[4].symbol == ']', "close bracket");
    ASSERT(m[5].symbol == 'F', "trailing F");
}

TEST(space_colonization_connectivity) {
    std::vector<Vec3> attractors;
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> uni(-1.0f, 1.0f);
    for (int i = 0; i < 200; ++i) {
        attractors.push_back({ uni(rng), 1.0f + uni(rng) * 0.5f, uni(rng) });
    }
    SpaceColonizationOptions opts;
    opts.attractionRadius = 1.5f;
    opts.killRadius = 0.15f;
    opts.segmentLength = 0.12f;
    opts.maxIterations = 80;
    opts.tropism = {0,1,0};
    opts.tropismWeight = 0.2f;
    auto segs = spaceColonize(attractors, {{0,0,0}}, {0,1,0}, opts);
    ASSERT(segs.size() > 1, "produced segments");
    // Every non-root segment's parent must be a valid earlier index.
    bool valid = true;
    int roots = 0;
    for (size_t i = 0; i < segs.size(); ++i) {
        if (segs[i].parent == -1) ++roots;
        else if (segs[i].parent < 0 || static_cast<size_t>(segs[i].parent) >= i) { valid = false; break; }
    }
    ASSERT(valid, "parent indices are well-formed");
    ASSERT(roots >= 1, "at least one root");

    thickenBranches(segs, 0.02f, 2.5f);
    // Every parent radius >= max child radius (monotonic non-decreasing).
    bool monotonic = true;
    std::vector<std::vector<int>> kids(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        if (segs[i].parent >= 0) kids[segs[i].parent].push_back(static_cast<int>(i));
    }
    for (size_t i = 0; i < segs.size(); ++i) {
        for (int c : kids[i]) {
            if (segs[i].radius + 1e-5f < segs[c].radius) { monotonic = false; break; }
        }
        if (!monotonic) break;
    }
    ASSERT(monotonic, "thicken: parent radius >= child radius");
}

TEST(mesh_branches_smoke) {
    // Three-segment Y: trunk + two children share a fork. meshBranches must
    // produce a non-empty merged mesh and not duplicate the parent ring at
    // the fork (chain-aware sweep).
    std::vector<BranchSegment> segs;
    BranchSegment root;
    root.parent = -1; root.from = {0, 0, 0}; root.to = {0, 0, 0};
    root.radius = 0.1f; root.depth = 0;
    segs.push_back(root);
    BranchSegment trunk;
    trunk.parent = 0; trunk.from = {0, 0, 0}; trunk.to = {0, 1, 0};
    trunk.radius = 0.08f; trunk.depth = 1;
    segs.push_back(trunk);
    BranchSegment left;
    left.parent = 1; left.from = {0, 1, 0}; left.to = {-0.5f, 1.7f, 0};
    left.radius = 0.04f; left.depth = 2;
    segs.push_back(left);
    BranchSegment right;
    right.parent = 1; right.from = {0, 1, 0}; right.to = {0.5f, 1.7f, 0};
    right.radius = 0.04f; right.depth = 2;
    segs.push_back(right);

    MeshData m = meshBranches(segs, 6);
    ASSERT(!m.empty(), "meshBranches Y mesh non-empty");
    ASSERT(m.vertexCount() > 0, "meshBranches has vertices");
}

static bool allFinite(const std::vector<float>& v) {
    for (float x : v) if (!std::isfinite(x)) return false;
    return true;
}

TEST(leaf_card_basic_grid) {
    LeafCardOptions o;
    o.widthSegments = 4;
    o.lengthSegments = 8;
    MeshData m = leafCard(LeafShape::Oval, o);
    ASSERT(!m.empty(), "leafCard non-empty");
    ASSERT(m.vertexCount() == 5u * 9u, "leafCard vertex count = (w+1)*(l+1)");
    ASSERT(m.triangleCount() == 4u * 8u * 2u, "leafCard tri count = w*l*2");
    ASSERT(m.hasNormals(), "leafCard has normals");
    ASSERT(m.hasUVs(), "leafCard has uvs");
    ASSERT(m.hasColors(), "leafCard has colors (windBend)");
    ASSERT(allFinite(m.positions), "leafCard positions finite");
    ASSERT(allFinite(m.normals), "leafCard normals finite");
    // Normals roughly unit length.
    bool unit = true;
    for (size_t i = 0; i < m.normals.size(); i += 3) {
        float L = std::sqrt(m.normals[i]*m.normals[i] + m.normals[i+1]*m.normals[i+1] + m.normals[i+2]*m.normals[i+2]);
        if (L < 0.9f || L > 1.1f) { unit = false; break; }
    }
    ASSERT(unit, "leafCard normals unit length");
    // Color R goes 0..1 along length.
    bool baseZero = m.colors[0] < 0.01f;
    size_t tipBase = (m.vertexCount() - 1u) * 4u;
    bool tipOne = m.colors[tipBase] > 0.99f;
    ASSERT(baseZero, "leafCard base windBend == 0");
    ASSERT(tipOne, "leafCard tip windBend == 1");
}

TEST(leaf_card_bend_curls_forward) {
    LeafCardOptions o;
    o.bend = 1.0f;
    MeshData straight = leafCard(LeafShape::Pointed, {});
    MeshData bent     = leafCard(LeafShape::Pointed, o);
    ASSERT(straight.vertexCount() == bent.vertexCount(), "bend preserves topology");
    // The bent mesh should have positive Y at the tip; the straight one stays at y=0.
    size_t tipIdx = (bent.vertexCount() - 1u) * 3u;
    ASSERT(bent.positions[tipIdx + 1] > 0.05f, "bent tip lifts in +Y");
    ASSERT(std::fabs(straight.positions[tipIdx + 1]) < 1e-5f, "straight tip at y=0");
}

TEST(flower_merges_petals_and_center) {
    FlowerOptions f;
    f.petalCount = 6;
    f.layers = 2;
    MeshData m = flower(f);
    ASSERT(!m.empty(), "flower non-empty");
    ASSERT(m.hasNormals(), "flower has normals");
    ASSERT(m.hasUVs(), "flower has uvs");
    ASSERT(m.hasColors(), "flower has colors");
    ASSERT(allFinite(m.positions), "flower positions finite");
    // 12 petals + 1 center; should be > 100 verts.
    ASSERT(m.vertexCount() > 100u, "flower has nontrivial vertex count");
    ASSERT(m.triangleCount() > 0u, "flower has triangles");
}

TEST(bezier_sweep_single_segment) {
    std::vector<Vec3> ctrl = {
        {0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}
    };
    std::vector<Vec2> profile = {{0.05f, 0}, {0, 0.05f}, {-0.05f, 0}, {0, -0.05f}};
    BezierSweepOptions o;
    o.samples = 16;
    MeshData m = bezierSweep(ctrl, profile, o);
    ASSERT(!m.empty(), "bezier sweep non-empty");
    ASSERT(m.hasNormals(), "bezier sweep normals");
    // 4 verts/ring * 16 rings + 2 cap centroids
    ASSERT(m.vertexCount() == 4u * 16u + 2u, "bezier sweep vertex count");
    ASSERT(allFinite(m.positions), "bezier sweep positions finite");
}

TEST(bezier_sweep_multi_segment) {
    // Two cubic segments sharing endpoint: 7 control points.
    std::vector<Vec3> ctrl = {
        {0, 0, 0}, {0.2f, 0.5f, 0}, {0.8f, 0.5f, 0}, {1.0f, 0.0f, 0},
        {1.2f, -0.5f, 0}, {1.8f, -0.5f, 0}, {2.0f, 0.0f, 0}
    };
    std::vector<Vec2> profile = {{0.04f, 0}, {0, 0.04f}, {-0.04f, 0}, {0, -0.04f}};
    BezierSweepOptions o;
    o.samples = 24;
    o.profileScale = {1.0f, 0.5f}; // taper
    MeshData m = bezierSweep(ctrl, profile, o);
    ASSERT(!m.empty(), "multi-segment bezier non-empty");
    ASSERT(m.vertexCount() == 4u * 24u + 2u, "multi-segment vertex count");
}

TEST(bezier_sweep_rejects_bad_input) {
    std::vector<Vec2> profile = {{0.1f, 0}, {0, 0.1f}, {-0.1f, 0}};
    // 3 control points -> not enough for a cubic.
    std::vector<Vec3> bad3 = {{0,0,0}, {0,1,0}, {0,2,0}};
    MeshData m1 = bezierSweep(bad3, profile);
    ASSERT(m1.empty(), "bezier rejects N<4");
    // 5 points -> (5-1)%3 != 0
    std::vector<Vec3> bad5 = {{0,0,0},{0,1,0},{0,2,0},{0,3,0},{0,4,0}};
    MeshData m2 = bezierSweep(bad5, profile);
    ASSERT(m2.empty(), "bezier rejects malformed segment count");
}

// --- Leaf scatter ----------------------------------------------------------

static std::vector<BranchSegment> yForkSegments() {
    // Same shape as mesh_branches_smoke: root marker, trunk, two children.
    std::vector<BranchSegment> segs;
    BranchSegment root;
    root.parent = -1; root.from = {0, 0, 0}; root.to = {0, 0, 0};
    root.radius = 0.10f; root.depth = 0;
    segs.push_back(root);
    BranchSegment trunk;
    trunk.parent = 0; trunk.from = {0, 0, 0}; trunk.to = {0, 1, 0};
    trunk.radius = 0.08f; trunk.depth = 1;
    segs.push_back(trunk);
    BranchSegment left;
    left.parent = 1; left.from = {0, 1, 0}; left.to = {-0.5f, 1.7f, 0};
    left.radius = 0.02f; left.depth = 2;
    segs.push_back(left);
    BranchSegment right;
    right.parent = 1; right.from = {0, 1, 0}; right.to = {0.5f, 1.7f, 0};
    right.radius = 0.02f; right.depth = 2;
    segs.push_back(right);
    return segs;
}

TEST(leaf_scatter_filters_by_radius_and_depth) {
    auto segs = yForkSegments();
    LeafPlacementOptions o;
    o.maxRadius = 0.05f;   // excludes trunk (0.08), root marker (0.10)
    o.minDepth  = 1;        // excludes root marker (depth 0)
    o.perUnitLength = 50.0f;
    o.seed = 1;
    LeafPlacements pl = placeLeavesOnBranches(segs, o);
    ASSERT(pl.count() > 0, "produced leaves on twigs");
    // Every accepted leaf must be on a thin enough branch.
    bool ok = true;
    for (size_t i = 0; i < pl.count(); ++i) {
        if (pl.branchRadius[i] > 0.05f + 1e-6f) { ok = false; break; }
        if (pl.branchDepth[i] < 1) { ok = false; break; }
    }
    ASSERT(ok, "all accepted leaves respect filters");

    // With a tighter maxRadius excluding the twigs too, count -> 0.
    o.maxRadius = 0.01f;
    LeafPlacements none = placeLeavesOnBranches(segs, o);
    ASSERT(none.count() == 0, "no leaves when nothing passes the filter");
}

TEST(leaf_scatter_density_scales_with_length) {
    // Single segment, depth 1, thin enough to pass default filter.
    auto makeSeg = [](float length) {
        std::vector<BranchSegment> segs;
        BranchSegment root;
        root.parent = -1; root.from = {0,0,0}; root.to = {0,0,0};
        root.radius = 0.0f; root.depth = 0;
        segs.push_back(root);
        BranchSegment s;
        s.parent = 0; s.from = {0,0,0}; s.to = {0, length, 0};
        s.radius = 0.02f; s.depth = 1;
        segs.push_back(s);
        return segs;
    };
    LeafPlacementOptions o;
    o.perUnitLength = 100.0f;
    o.seed = 7;
    o.dedupRadius = 0.0f;
    o.tiltJitter = 0.0f;
    o.rollJitter = 0.0f;
    o.scaleJitter = 0.0f;

    auto a = placeLeavesOnBranches(makeSeg(1.0f), o);
    auto b = placeLeavesOnBranches(makeSeg(2.0f), o);
    ASSERT(a.count() > 0 && b.count() > 0, "both produced leaves");
    // Allow ±20% slop from stochastic rounding.
    float ratio = static_cast<float>(b.count()) / static_cast<float>(a.count());
    ASSERT(ratio > 1.6f && ratio < 2.4f, "doubling length roughly doubles leaf count");
}

TEST(leaf_scatter_deterministic) {
    auto segs = yForkSegments();
    LeafPlacementOptions o;
    o.maxRadius = 0.05f;
    o.perUnitLength = 30.0f;
    o.seed = 123;
    LeafPlacements a = placeLeavesOnBranches(segs, o);
    LeafPlacements b = placeLeavesOnBranches(segs, o);
    ASSERT(a.count() == b.count(), "same seed: same leaf count");
    bool same = (a.transforms.size() == b.transforms.size());
    if (same) {
        for (size_t i = 0; i < a.transforms.size(); ++i) {
            if (a.transforms[i] != b.transforms[i]) { same = false; break; }
        }
    }
    ASSERT(same, "same seed: identical transforms");

    o.seed = 456;
    LeafPlacements c = placeLeavesOnBranches(segs, o);
    bool different = (c.transforms.size() != a.transforms.size());
    if (!different && !a.transforms.empty()) {
        for (size_t i = 0; i < a.transforms.size(); ++i) {
            if (std::fabs(c.transforms[i] - a.transforms[i]) > 1e-5f) {
                different = true; break;
            }
        }
    }
    ASSERT(different, "different seed: diverges");
}

TEST(leaf_scatter_dedup_enforces_spacing) {
    // One thin segment, packed densely, with dedup on. Origins must respect spacing.
    std::vector<BranchSegment> segs;
    BranchSegment root;
    root.parent = -1; root.from = {0,0,0}; root.to = {0,0,0};
    root.radius = 0.0f; root.depth = 0;
    segs.push_back(root);
    BranchSegment s;
    s.parent = 0; s.from = {0,0,0}; s.to = {0, 1.0f, 0};
    s.radius = 0.02f; s.depth = 1;
    segs.push_back(s);

    LeafPlacementOptions o;
    o.perUnitLength = 200.0f;
    o.dedupRadius = 0.05f;
    o.seed = 99;
    LeafPlacements pl = placeLeavesOnBranches(segs, o);
    ASSERT(pl.count() > 0, "dedup leaves something behind");
    // Brute-force check spacing.
    bool ok = true;
    for (size_t i = 0; i < pl.count() && ok; ++i) {
        Vec3 pi{pl.transforms[i*16+12], pl.transforms[i*16+13], pl.transforms[i*16+14]};
        for (size_t j = i + 1; j < pl.count(); ++j) {
            Vec3 pj{pl.transforms[j*16+12], pl.transforms[j*16+13], pl.transforms[j*16+14]};
            if (vdist(pi, pj) < o.dedupRadius - 1e-5f) { ok = false; break; }
        }
    }
    ASSERT(ok, "all leaf origins respect dedupRadius");
}

TEST(leaf_scatter_smoke_y_fork) {
    auto segs = yForkSegments();
    MeshData leaf = leafCard(LeafShape::Oval, {});
    LeafPlacementOptions o;
    o.maxRadius = 0.05f;
    o.perUnitLength = 20.0f;
    o.baseScale = 0.3f;
    o.seed = 42;
    MeshData m = scatterLeaves(segs, leaf, o);
    ASSERT(!m.empty(), "scatterLeaves non-empty");
    ASSERT(m.hasNormals(), "scattered leaves have normals");
    ASSERT(m.hasUVs(), "scattered leaves have UVs");
    ASSERT(m.hasColors(), "scattered leaves have colors (windBend)");
    ASSERT(allFinite(m.positions), "positions finite");
    ASSERT(allFinite(m.normals), "normals finite");
    // Every triangle is a unique-vertex triangle.
    bool degen = false;
    for (size_t t = 0; t < m.triangleCount(); ++t) {
        uint32_t a = m.indices[t*3], b = m.indices[t*3+1], c = m.indices[t*3+2];
        if (a == b || b == c || a == c) { degen = true; break; }
    }
    ASSERT(!degen, "no degenerate triangles in scattered output");
}

TEST(leaf_scatter_terminal_only) {
    auto segs = yForkSegments();
    LeafPlacementOptions o;
    o.maxRadius = 1.0f;     // pass everything by radius
    o.minDepth = 0;          // and by depth
    o.terminalOnly = true;   // only chain tips
    o.perUnitLength = 30.0f;
    o.seed = 5;
    LeafPlacements pl = placeLeavesOnBranches(segs, o);
    ASSERT(pl.count() > 0, "terminal-only still produces leaves on the two tips");
    bool ok = true;
    for (size_t i = 0; i < pl.count(); ++i) {
        // Tips are at depth 2 in our Y-fork.
        if (pl.branchDepth[i] != 2) { ok = false; break; }
    }
    ASSERT(ok, "terminalOnly selects only the Y-fork tips");
}

TEST(blade_strip_explicit_path) {
    // Straight upright path of 9 points (segs=8). Default opts: capStart
    // off, capEnd on, miterJoints on. Expect (segs+1)·4 ring vertices plus
    // a single capEnd centroid (since the diamond profile has 4 verts the
    // cap is a 4-tri fan from the centroid).
    std::vector<Vec3> path;
    const int segs = 8;
    for (int i = 0; i <= segs; ++i) {
        path.push_back({0.0f, 0.1f * i, 0.0f});
    }
    BladeStripOptions opts;
    opts.width = 0.05f;
    opts.thickness = 0.01f;
    MeshData m = bladeStrip(path, opts);
    ASSERT(!m.empty(), "bladeStrip produced output");
    ASSERT(m.vertexCount() == (size_t)((segs + 1) * 4 + 1),
           "ring verts + 1 cap centroid");
    ASSERT(m.hasNormals(), "bladeStrip has normals");
    ASSERT(allFinite(m.positions), "positions finite");
    ASSERT(allFinite(m.normals), "normals finite");
}

TEST(blade_strip_linear_taper_tip_radius) {
    // Linear taper from 1.0 at base to 0.1 at tip. The last ring's vertex
    // along profile-X should sit at width · 0.1.
    std::vector<Vec3> path;
    const int segs = 8;
    for (int i = 0; i <= segs; ++i) {
        path.push_back({0.0f, 0.1f * i, 0.0f});
    }
    BladeStripOptions opts;
    opts.width = 0.05f;
    opts.thickness = 0.0f;
    opts.profileScale.resize(path.size());
    for (size_t i = 0; i < path.size(); ++i) {
        float t = (float)i / (float)(path.size() - 1);
        opts.profileScale[i] = 1.0f - 0.9f * t;
    }
    MeshData m = bladeStrip(path, opts);
    ASSERT(!m.empty(), "tapered bladeStrip non-empty");
    // Last ring starts at index segs · 4. Profile vertex 0 is (+width, 0)
    // in profile-XY → first vertex of every ring lands on the lateral axis
    // at distance width · scaleAt(ring). For a vertical path, the lateral
    // axis lands on world +X (sweep's anyPerpendicular fallback).
    size_t lastRingFirst = (size_t)segs * 4;
    float x = m.positions[lastRingFirst * 3 + 0];
    float z = m.positions[lastRingFirst * 3 + 2];
    float r = std::sqrt(x*x + z*z);
    float expected = 0.05f * 0.1f;
    ASSERT(std::fabs(r - expected) < 1e-5f,
           "tip ring radius matches profileScale.back() · width");
}

TEST(blade_strip_zero_thickness) {
    std::vector<Vec3> path;
    const int segs = 6;
    for (int i = 0; i <= segs; ++i) {
        path.push_back({0.0f, 0.1f * i, 0.0f});
    }
    BladeStripOptions opts;
    opts.width = 0.05f;
    opts.thickness = 0.0f;
    MeshData m = bladeStrip(path, opts);
    ASSERT(!m.empty(), "zero-thickness bladeStrip non-empty");
    ASSERT(allFinite(m.positions), "zero-thickness positions finite");
    ASSERT(allFinite(m.normals), "zero-thickness normals finite");
}

TEST(blade_path_endpoints) {
    BladePathOptions opts;
    opts.base = {1.0f, 2.0f, 3.0f};
    opts.tipDir = {0.0f, 1.0f, 0.0f};
    opts.length = 5.0f;
    opts.bend = 0.0f;
    opts.lift = 0.0f;
    opts.segments = 8;
    auto pts = bladePath(opts);
    ASSERT(pts.size() == 9u, "segments+1 points");
    ASSERT(vdist(pts.front(), opts.base) < 1e-5f, "starts at base");
    Vec3 expectedTip = {1.0f, 7.0f, 3.0f};
    ASSERT(vdist(pts.back(), expectedTip) < 1e-5f, "ends at base + tipDir·length");
}

TEST(tree_archetype_smoke) {
    TreeOptions opts;
    opts.base = {0.0f, 0.0f, 0.0f};
    opts.canopyCenter = {0.0f, 4.0f, 0.0f};
    opts.canopyRadius = 2.0f;
    opts.attractorCount = 80;
    opts.sides = 6;
    opts.leafRadius = 0.04f;
    opts.pipeExp = 2.5f;
    opts.colonize.attractionRadius = 6.0f;   // base must see canopy attractors
    opts.colonize.killRadius = 0.3f;
    opts.colonize.segmentLength = 0.25f;
    opts.colonize.maxIterations = 300;
    opts.colonize.tropism = {0.0f, 1.0f, 0.0f};
    opts.colonize.tropismWeight = 0.4f;
    opts.seed = 42;

    TreeResult r = tree(opts);
    ASSERT(!r.segments.empty(), "tree produced segments");
    ASSERT(!r.branches.empty(), "tree produced a branch mesh");
    ASSERT(r.branches.hasNormals(), "branch mesh has normals");
    // First (root) segment should originate at base.
    ASSERT(vdist(r.segments.front().from, opts.base) < 1e-3f,
           "root segment starts at base");
    // Some branch should reach near the canopy.
    bool reachedCanopy = false;
    for (const auto& s : r.segments) {
        if (vdist(s.to, opts.canopyCenter) < opts.canopyRadius * 1.2f) {
            reachedCanopy = true;
            break;
        }
    }
    ASSERT(reachedCanopy, "at least one segment reaches the canopy region");
    ASSERT(allFinite(r.branches.positions), "branch positions finite");
    ASSERT(allFinite(r.branches.normals), "branch normals finite");
}

TEST(tree_archetype_deterministic) {
    TreeOptions opts;
    opts.canopyCenter = {0.0f, 4.0f, 0.0f};
    opts.canopyRadius = 2.0f;
    opts.attractorCount = 50;
    opts.seed = 7;
    auto a = tree(opts);
    auto b = tree(opts);
    ASSERT(a.segments.size() == b.segments.size(), "deterministic segment count");
    ASSERT(a.branches.vertexCount() == b.branches.vertexCount(),
           "deterministic mesh vertex count");
}

TEST(blade_path_bend_offsets_midpoint) {
    BladePathOptions opts;
    opts.base = {0.0f, 0.0f, 0.0f};
    opts.tipDir = {0.0f, 1.0f, 0.0f};
    opts.length = 1.0f;
    opts.bend = 0.5f;
    opts.lift = 0.0f;
    opts.segments = 8;
    auto pts = bladePath(opts);
    // Midpoint (index 4 of 9) should have positive x (lateral axis falls
    // on world +X for a +Y tipDir).
    ASSERT(pts[4].x > 0.1f, "bend pulls midpoint laterally");
}

// ── CapsuleField ─────────────────────────────────────────────────────────
// Occupancy primitive for obstacle-aware scatter / colonize.

TEST(capsule_field_distance_basic) {
    // Single capsule along Y from (0,0,0) to (0,1,0), radius 0.1.
    Capsule c; c.a = {0,0,0}; c.b = {0,1,0}; c.radius = 0.1f; c.tag = 0;
    CapsuleField f({c});

    // On the axis, midway: distance = -radius.
    float d_in = f.distance({0, 0.5f, 0});
    ASSERT(std::fabs(d_in + 0.1f) < 1e-4f, "inside axis distance == -radius");

    // 0.3 units perpendicular to the axis: distance = 0.3 - 0.1 = 0.2.
    float d_out = f.distance({0.3f, 0.5f, 0});
    ASSERT(std::fabs(d_out - 0.2f) < 1e-4f, "perpendicular distance correct");

    // Past the cap on +Y by 0.3: distance to nearest cap point is 0.3 - 0.1.
    float d_cap = f.distance({0, 1.3f, 0});
    ASSERT(std::fabs(d_cap - 0.2f) < 1e-4f, "endcap distance correct");
}

TEST(capsule_field_excludes_tag) {
    // Two overlapping capsules tagged 0 and 1, both containing the origin.
    Capsule c0; c0.a = {-1,0,0}; c0.b = {1,0,0}; c0.radius = 0.2f; c0.tag = 0;
    Capsule c1; c1.a = {0,-1,0}; c1.b = {0,1,0}; c1.radius = 0.2f; c1.tag = 1;
    CapsuleField f({c0, c1});

    // Without exclusion, origin is well inside.
    ASSERT(f.distance({0,0,0}, -1) < 0.0f, "no exclude: inside");
    // Exclude tag 0 → still inside via c1.
    ASSERT(f.distance({0,0,0}, 0) < 0.0f, "exclude 0: still inside via 1");
    // Exclude tag 1 → still inside via c0.
    ASSERT(f.distance({0,0,0}, 1) < 0.0f, "exclude 1: still inside via 0");

    // Point on the X axis only inside c0: excluding 0 makes it outside.
    Vec3 p{0.5f, 0, 0};   // inside c0 (radius 0.2 around X axis), outside c1
    ASSERT(f.distance(p, -1) < 0.0f, "p inside c0");
    ASSERT(f.distance(p, 0) > 0.0f, "exclude 0 → outside");
}

TEST(capsule_field_intersect_sphere) {
    Capsule c; c.a = {0,0,0}; c.b = {0,1,0}; c.radius = 0.1f;
    CapsuleField f({c});

    ASSERT(f.intersectsSphere({0, 0.5f, 0}, 0.05f),
           "sphere centered on axis intersects");
    ASSERT(!f.intersectsSphere({1.0f, 0.5f, 0}, 0.05f),
           "far sphere does not intersect");
    // Tangent: surface at x=0.1, sphere radius 0.05 at x=0.16 should not
    // intersect (gap = 0.06 - 0.05 = 0.01).
    ASSERT(!f.intersectsSphere({0.16f, 0.5f, 0}, 0.05f),
           "near-tangent miss");
    ASSERT(f.intersectsSphere({0.14f, 0.5f, 0}, 0.05f),
           "near-tangent hit");
}

TEST(capsule_field_from_segments) {
    // Two perpendicular branches forming an L, plus a synthetic root.
    std::vector<BranchSegment> segs;
    BranchSegment root; root.parent = -1; root.from = {0,0,0}; root.to = {0,0,0};
    root.radius = 0.0f; root.depth = 0;
    segs.push_back(root);
    BranchSegment a; a.parent = 0; a.from = {0,0,0}; a.to = {0,1,0};
    a.radius = 0.05f; a.depth = 1; segs.push_back(a);
    BranchSegment b; b.parent = 1; b.from = {0,1,0}; b.to = {1,1,0};
    b.radius = 0.04f; b.depth = 2; segs.push_back(b);

    auto caps = CapsuleField::capsulesFromSegments(segs);
    // Root's from==to and radius==0 → skipped.
    ASSERT(caps.size() == 2, "root segment skipped");
    ASSERT(caps[0].tag == 1 && caps[1].tag == 2, "tags equal segment indices");

    CapsuleField f(std::move(caps));
    // Point on segment a's axis sits inside that capsule.
    ASSERT(f.distance({0, 0.5f, 0}) < 0.0f, "on a-axis: inside");
    // Far away: positive.
    ASSERT(f.distance({5, 5, 5}) > 0.0f, "far point: outside");
    // Excluding the segment-a tag, that same point is outside.
    ASSERT(f.distance({0, 0.5f, 0}, 1) > 0.0f, "exclude segment a: outside");
}

// ── Obstacle-aware scatter ───────────────────────────────────────────────

TEST(leaf_scatter_avoids_capsule) {
    // Single straight branch along Y; obstacle capsule blocks the upper half.
    std::vector<BranchSegment> segs;
    BranchSegment root; root.parent = -1; root.from = {0,0,0}; root.to = {0,0,0};
    root.radius = 0.0f; root.depth = 0; segs.push_back(root);
    BranchSegment s; s.parent = 0; s.from = {0,0,0}; s.to = {0,1,0};
    s.radius = 0.02f; s.depth = 1; segs.push_back(s);

    Capsule blocker;
    // Cylinder centred at y=0.75 covering the upper quarter, well clear of
    // the branch's own segment (tag = 1).
    blocker.a = {-0.5f, 0.75f, 0};
    blocker.b = { 0.5f, 0.75f, 0};
    blocker.radius = 0.10f;
    blocker.tag = 999;   // not a real segment index
    CapsuleField field({blocker});

    LeafPlacementOptions o;
    o.maxRadius = 0.05f;
    o.perUnitLength = 200.0f;
    o.tiltJitter = 0; o.rollJitter = 0; o.scaleJitter = 0;
    o.seed = 7;

    LeafPlacements before = placeLeavesOnBranches(segs, o);
    ASSERT(before.count() > 0, "baseline produced leaves");

    // Without avoid, some origins land inside the blocker volume.
    int inside_before = 0;
    for (size_t i = 0; i < before.count(); ++i) {
        Vec3 p{before.transforms[i*16+12], before.transforms[i*16+13], before.transforms[i*16+14]};
        if (field.tooClose(p, 0.0f, 1)) ++inside_before;
    }
    ASSERT(inside_before > 0, "without avoid: collisions exist");

    o.avoid = &field;
    LeafPlacements after = placeLeavesOnBranches(segs, o);
    int inside_after = 0;
    for (size_t i = 0; i < after.count(); ++i) {
        Vec3 p{after.transforms[i*16+12], after.transforms[i*16+13], after.transforms[i*16+14]};
        if (field.tooClose(p, 0.0f, 1)) ++inside_after;
    }
    ASSERT(inside_after == 0, "with avoid: no collisions");
}

TEST(leaf_scatter_keepout_spheres) {
    std::vector<BranchSegment> segs;
    BranchSegment root; root.parent = -1; root.from = {0,0,0}; root.to = {0,0,0};
    root.radius = 0.0f; root.depth = 0; segs.push_back(root);
    BranchSegment s; s.parent = 0; s.from = {0,0,0}; s.to = {0,1,0};
    s.radius = 0.02f; s.depth = 1; segs.push_back(s);

    Sphere keep;
    keep.center = {0, 0.5f, 0};
    keep.radius = 0.15f;

    LeafPlacementOptions o;
    o.maxRadius = 0.05f;
    o.perUnitLength = 200.0f;
    o.seed = 7;
    o.keepOut.push_back(keep);

    LeafPlacements pl = placeLeavesOnBranches(segs, o);
    ASSERT(pl.count() > 0, "keepOut leaves something behind");
    bool ok = true;
    for (size_t i = 0; i < pl.count(); ++i) {
        Vec3 p{pl.transforms[i*16+12], pl.transforms[i*16+13], pl.transforms[i*16+14]};
        Vec3 d = p - keep.center;
        if (vdot(d, d) <= keep.radius * keep.radius) { ok = false; break; }
    }
    ASSERT(ok, "no leaf origin inside keep-out sphere");
}

TEST(leaf_scatter_pushout_recovers_some) {
    // Without pushout we drop everything inside the obstacle. With pushout,
    // those candidates are nudged outward and survive — so the count goes up
    // (or at least does not go down) and no surviving placement remains
    // inside.
    std::vector<BranchSegment> segs;
    BranchSegment root; root.parent = -1; root.from = {0,0,0}; root.to = {0,0,0};
    root.radius = 0.0f; root.depth = 0; segs.push_back(root);
    BranchSegment s; s.parent = 0; s.from = {0,0,0}; s.to = {0,1,0};
    s.radius = 0.02f; s.depth = 1; segs.push_back(s);

    Capsule blocker;
    blocker.a = {-0.5f, 0.5f, 0};
    blocker.b = { 0.5f, 0.5f, 0};
    blocker.radius = 0.05f;
    blocker.tag = 999;
    CapsuleField field({blocker});

    LeafPlacementOptions o;
    o.maxRadius = 0.05f;
    o.perUnitLength = 200.0f;
    o.tiltJitter = 0; o.rollJitter = 0; o.scaleJitter = 0;
    o.seed = 11;
    o.avoid = &field;

    LeafPlacements hard = placeLeavesOnBranches(segs, o);

    o.obstaclePushout = 0.05f;
    LeafPlacements push = placeLeavesOnBranches(segs, o);
    ASSERT(push.count() >= hard.count(), "pushout retains at least as many");

    bool ok = true;
    for (size_t i = 0; i < push.count(); ++i) {
        Vec3 p{push.transforms[i*16+12], push.transforms[i*16+13], push.transforms[i*16+14]};
        if (field.tooClose(p, 0.0f, 1)) { ok = false; break; }
    }
    ASSERT(ok, "after pushout no placement remains inside obstacle");
}

// ── packAnchors ──────────────────────────────────────────────────────────

TEST(pack_anchors_min_spacing) {
    // 100 candidates 0.01 apart along X. With minSpacing 0.05 we expect
    // roughly one accept per 5 candidates.
    std::vector<Vec3> cand;
    cand.reserve(100);
    for (int i = 0; i < 100; ++i) cand.push_back({i * 0.01f, 0, 0});

    AnchorPackOptions o;
    o.minSpacing = 0.05f;
    o.seed = 42;
    auto idx = packAnchors(cand, nullptr, {}, o);
    // 100 candidates over [0, 0.99] with minSpacing 0.05 gives a hard upper
    // bound around 20 and (with random visit order) a lower bound around 10
    // when an unlucky shuffle leaves gaps. Be generous on the lower side.
    ASSERT(idx.size() >= 8 && idx.size() <= 25, "spacing keeps count in expected band");

    // Verify pairwise spacing.
    bool ok = true;
    for (size_t i = 0; i < idx.size() && ok; ++i) {
        for (size_t j = i + 1; j < idx.size(); ++j) {
            float d = vdist(cand[idx[i]], cand[idx[j]]);
            if (d < o.minSpacing - 1e-5f) { ok = false; break; }
        }
    }
    ASSERT(ok, "all accepted anchors respect minSpacing");
}

TEST(pack_anchors_max_count) {
    std::vector<Vec3> cand;
    for (int i = 0; i < 50; ++i) cand.push_back({i * 1.0f, 0, 0});
    AnchorPackOptions o;
    o.maxCount = 5;
    o.seed = 1;
    auto idx = packAnchors(cand, nullptr, {}, o);
    ASSERT(idx.size() == 5, "maxCount caps acceptances");
}

TEST(pack_anchors_avoid_obstacle) {
    // Half of the candidates lie inside an obstacle capsule and should be
    // rejected.
    std::vector<Vec3> cand;
    for (int i = 0; i < 20; ++i) cand.push_back({i * 0.05f, 0, 0});
    Capsule c; c.a = {-1, 0, 0}; c.b = {0.5f, 0, 0}; c.radius = 0.05f; c.tag = 0;
    CapsuleField field({c});

    AnchorPackOptions o;
    o.minObstacleDistance = 0.0f;
    o.seed = 3;
    auto idx = packAnchors(cand, &field, {}, o);

    bool ok = true;
    for (int i : idx) {
        if (field.tooClose(cand[i], 0.0f, -1)) { ok = false; break; }
    }
    ASSERT(ok, "no accepted anchor lies inside obstacle");
    ASSERT(idx.size() > 0, "some anchors survive outside the obstacle");
}

// ── Obstacle-aware spaceColonize ──────────────────────────────────────────

TEST(space_colonize_obstacle_blocks_growth) {
    // Single attractor on +X axis; an obstacle capsule blocks the direct
    // path. With hard reject (no steer) the tree cannot reach the
    // attractor — no segment crosses the obstacle.
    std::vector<Vec3> attractors{{2.0f, 0.0f, 0.0f}};
    std::vector<Vec3> seeds{{0.0f, 0.0f, 0.0f}};

    Capsule blocker;
    blocker.a = {1.0f, -2.0f, 0.0f};
    blocker.b = {1.0f,  2.0f, 0.0f};
    blocker.radius = 0.3f;
    blocker.tag = -1;
    CapsuleField field({blocker});

    SpaceColonizationOptions opts;
    opts.attractionRadius = 5.0f;
    opts.killRadius       = 0.2f;
    opts.segmentLength    = 0.2f;
    opts.maxIterations    = 60;
    opts.obstacles        = &field;
    opts.obstacleClearance = 0.0f;
    opts.obstacleSteer     = 0.0f;   // hard reject

    auto segs = spaceColonize(attractors, seeds, {1, 0, 0}, opts);
    bool anyCross = false;
    for (const BranchSegment& s : segs) {
        if (vdist2(s.from, s.to) < 1e-12f) continue;
        Vec3 mid = (s.from + s.to) * 0.5f;
        if (field.tooClose(mid, 0.0f)) { anyCross = true; break; }
        if (field.tooClose(s.to, 0.0f))  { anyCross = true; break; }
    }
    ASSERT(!anyCross, "no segment lies inside the obstacle");
}

TEST(space_colonize_obstacle_steer_around) {
    // Same setup, but with steering. Tree should bend around the obstacle
    // and produce more segments than the hard-reject case.
    std::vector<Vec3> attractors{{2.0f, 0.0f, 0.0f}};
    std::vector<Vec3> seeds{{0.0f, 0.0f, 0.0f}};

    Capsule blocker;
    blocker.a = {1.0f, -2.0f, 0.0f};
    blocker.b = {1.0f,  2.0f, 0.0f};
    blocker.radius = 0.3f;
    CapsuleField field({blocker});

    SpaceColonizationOptions opts;
    opts.attractionRadius = 5.0f;
    opts.killRadius       = 0.2f;
    opts.segmentLength    = 0.2f;
    opts.maxIterations    = 100;
    opts.obstacles        = &field;
    opts.obstacleSteer    = 0.8f;   // ~46°

    auto segs = spaceColonize(attractors, seeds, {1, 0, 0}, opts);

    // No segment ends inside the obstacle.
    bool anyInside = false;
    for (const BranchSegment& s : segs) {
        if (vdist2(s.from, s.to) < 1e-12f) continue;
        if (field.tooClose(s.to, 0.0f)) { anyInside = true; break; }
    }
    ASSERT(!anyInside, "steering keeps every endpoint outside obstacle");

    // And the tree managed to grow past the obstacle (some endpoint with
    // x > 1.4 — past the blocker's outer surface).
    bool reached = false;
    for (const BranchSegment& s : segs) {
        if (s.to.x > 1.4f) { reached = true; break; }
    }
    ASSERT(reached, "steering let the tree grow past the obstacle");
}

TEST(space_colonize_obstacle_null_is_noop) {
    // Identical inputs, no obstacles: result must match the prior
    // behavior bit-exact (regression guard so the obstacle plumb-through
    // doesn't perturb the deterministic output).
    std::vector<Vec3> attractors{
        {1, 1, 0}, {-1, 1, 0}, {0, 1.5f, 1}, {0, 1.5f, -1}, {0, 2, 0}};
    std::vector<Vec3> seeds{{0, 0, 0}};
    SpaceColonizationOptions opts;
    opts.attractionRadius = 1.5f;
    opts.killRadius       = 0.2f;
    opts.segmentLength    = 0.2f;
    opts.maxIterations    = 50;

    auto a = spaceColonize(attractors, seeds, {0, 1, 0}, opts);

    opts.obstacles = nullptr;
    auto b = spaceColonize(attractors, seeds, {0, 1, 0}, opts);
    ASSERT(a.size() == b.size(), "null obstacles: same segment count");
    bool same = true;
    for (size_t i = 0; i < a.size() && same; ++i) {
        if (vdist2(a[i].from, b[i].from) > 1e-10f) same = false;
        if (vdist2(a[i].to,   b[i].to)   > 1e-10f) same = false;
        if (a[i].parent != b[i].parent) same = false;
    }
    ASSERT(same, "null obstacles: identical segments");
}

// ── L-system turtle ───────────────────────────────────────────────────────

TEST(lsystem_turtle_straight_line) {
    auto mods = parseModules("FFFF");
    TurtleOptions to;
    to.stepLength = 1.0f;
    auto segs = lsystemToBranches(mods, to);
    // Synthetic root + 4 forward segments.
    ASSERT(segs.size() == 5, "1 root + 4 F segments");
    ASSERT(segs[0].parent == -1, "root parent = -1");
    // Each F endpoint advances by 1 along +Y (default heading).
    for (int i = 1; i <= 4; ++i) {
        ASSERT(std::fabs(segs[i].to.y - static_cast<float>(i)) < 1e-4f, "F advances along +Y");
        ASSERT(segs[i].parent == i - 1, "chain parents");
    }
}

TEST(lsystem_turtle_branching) {
    // F[+F]F → main F, branch +F off the tip, then continue with F.
    auto mods = parseModules("F[+F]F");
    TurtleOptions to;
    to.stepLength = 1.0f;
    to.angle = 3.14159265358979323846f / 4.0f;   // 45°
    auto segs = lsystemToBranches(mods, to);

    // 1 root + 3 F segments.
    ASSERT(segs.size() == 4, "1 root + 3 F segments");

    // Find the bracketed +F. It must parent off segment index 1 (the first
    // F's endpoint) and not lie along +Y.
    int branchIdx = -1;
    for (size_t i = 1; i < segs.size(); ++i) {
        if (segs[i].parent == 1 && std::fabs(segs[i].to.x) > 0.1f) {
            branchIdx = static_cast<int>(i);
            break;
        }
    }
    ASSERT(branchIdx >= 0, "branch parents off first F's tip");

    // The post-bracket F continues from the first F's tip too, going
    // straight up (still along +Y after the [..] pop).
    int continuationIdx = -1;
    for (size_t i = 1; i < segs.size(); ++i) {
        if (static_cast<int>(i) == branchIdx) continue;
        if (segs[i].parent == 1) { continuationIdx = static_cast<int>(i); break; }
    }
    ASSERT(continuationIdx >= 0, "continuation parents off first F too");
    ASSERT(std::fabs(segs[continuationIdx].to.x) < 1e-4f, "continuation stays on Y axis");
    ASSERT(segs[continuationIdx].to.y > 1.5f, "continuation moves further up");
}

TEST(lsystem_turtle_param_overrides) {
    // F(0.5)+(45)F(0.5) — half-step, 45° yaw, half-step.
    auto mods = parseModules("F(0.5)+(45)F(0.5)");
    TurtleOptions to;
    to.stepLength = 1.0f;   // overridden by params
    auto segs = lsystemToBranches(mods, to);
    ASSERT(segs.size() == 3, "1 root + 2 F");
    ASSERT(std::fabs(segs[1].to.y - 0.5f) < 1e-4f, "first F is 0.5 long");
    // After F(0.5) heading is +Y; +(45) yaws +45° around `up` (=+Z by default).
    // Right-hand rule: +Y rotates toward -X. So heading becomes
    // (-sin45, cos45, 0). The second F lands at the first F's tip plus
    // 0.5·heading.
    float c45 = std::cos(3.14159265358979323846f / 4.0f);
    float s45 = std::sin(3.14159265358979323846f / 4.0f);
    Vec3 expected{-0.5f * s45, 0.5f + 0.5f * c45, 0.0f};
    ASSERT(vdist(segs[2].to, expected) < 1e-3f, "second F lands at 45° offset");
}

TEST(lsystem_turtle_radius_command) {
    auto mods = parseModules("!(0.05)F!(0.02)F");
    TurtleOptions to;
    to.stepLength = 1.0f;
    to.radius = 0.10f;
    auto segs = lsystemToBranches(mods, to);
    ASSERT(segs.size() == 3, "1 root + 2 F");
    ASSERT(std::fabs(segs[1].radius - 0.05f) < 1e-6f, "first F uses ! radius 0.05");
    ASSERT(std::fabs(segs[2].radius - 0.02f) < 1e-6f, "second F uses ! radius 0.02");
}

