#include "test_framework.h"

#include "bromesh/optimization/spatial_hash.h"
#include "bromesh/manipulation/sweep.h"
#include "bromesh/procedural/lsystem.h"
#include "bromesh/procedural/space_colonization.h"
#include "bromesh/procedural/plants/tree.h"
#include "bromesh/procedural/plants/conifer.h"
#include "bromesh/procedural/plants/shrub.h"
#include "bromesh/procedural/plants/grass_tuft.h"
#include "bromesh/procedural/plants/vine.h"
#include "bromesh/procedural/plants/fern.h"
#include "bromesh/procedural/plants/succulent.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <set>

using namespace bromesh;

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

TEST(plant_tree_smoke) {
    TreeParams p;
    p.seed = 1;
    p.height = 4.0f;
    p.canopyRadius = 1.8f;
    p.attractorCount = 200;
    PlantResult r = buildTree(p);
    ASSERT(!r.branchMesh.empty(), "tree mesh non-empty");
    // Canopy must actually grow. Trunk-only output (which produces exactly
    // one terminal leaf at the trunk tip) is a regression: it means
    // colonization couldn't reach the attractor cloud.
    ASSERT(r.leaves.size() >= 10, "canopy grew (>=10 leaves)");
    // AABB top must be near the requested height. Trunk-only output stops
    // around 0.45 * height where the pre-grown trunk ends.
    ASSERT(r.aabbMax.y >= p.height * 0.7f, "canopy reaches >=70% of height");
    ASSERT(r.aabbMin.x < r.aabbMax.x, "aabb x sane");
    ASSERT(r.aabbMin.y < r.aabbMax.y, "aabb y sane");
    ASSERT(r.aabbMin.z < r.aabbMax.z, "aabb z sane");
}

TEST(plant_tree_age_monotonic) {
    // Same seed at different ages must produce a structurally consistent
    // tree: more segments / taller AABB as age increases.
    TreeParams p; p.seed = 42; p.height = 6.0f; p.attractorCount = 400;
    p.age01 = 0.3f; PlantResult r03 = buildTree(p);
    p.age01 = 0.6f; PlantResult r06 = buildTree(p);
    p.age01 = 1.0f; PlantResult r10 = buildTree(p);
    ASSERT(!r10.branchMesh.empty(), "mature tree non-empty");
    ASSERT(r10.branchMesh.vertexCount() >= r06.branchMesh.vertexCount(),
           "mature tree has >= verts than mid-age");
    ASSERT(r06.branchMesh.vertexCount() >= r03.branchMesh.vertexCount(),
           "mid-age tree has >= verts than young");
    ASSERT(r10.aabbMax.y >= r03.aabbMax.y * 0.95f,
           "mature aabb top >= young aabb top");
}

TEST(plant_conifer_smoke) {
    ConiferParams p; p.seed = 2;
    PlantResult r = buildConifer(p);
    ASSERT(!r.branchMesh.empty(), "conifer mesh non-empty");
}

TEST(plant_shrub_smoke) {
    ShrubParams p; p.seed = 3;
    PlantResult r = buildShrub(p);
    ASSERT(!r.branchMesh.empty(), "shrub mesh non-empty");
}

TEST(plant_grass_smoke) {
    GrassTuftParams p; p.seed = 4;
    PlantResult r = buildGrassTuft(p);
    ASSERT(!r.branchMesh.empty(), "grass mesh non-empty");
}

TEST(plant_vine_smoke) {
    VineParams p; p.seed = 5;
    PlantResult r = buildVine(p);
    ASSERT(!r.branchMesh.empty(), "vine mesh non-empty");
    ASSERT(r.leaves.size() > 0, "vine leaves non-empty");
}

TEST(plant_fern_smoke) {
    FernParams p; p.seed = 6;
    PlantResult r = buildFern(p);
    ASSERT(!r.branchMesh.empty(), "fern mesh non-empty");
}

TEST(plant_succulent_smoke) {
    SucculentParams p; p.seed = 7;
    PlantResult r = buildSucculent(p);
    ASSERT(!r.branchMesh.empty(), "succulent mesh non-empty");
}
