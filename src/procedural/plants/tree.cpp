#include "bromesh/procedural/plants/tree.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

// Hash a 64-bit seed and a small integer index into a deterministic
// float in [0, 1). Used to jitter per-segment birth times so siblings
// don't all appear at the same instant of growth.
static inline float hash01(uint64_t seed, int idx) {
    uint64_t x = seed ^ (uint64_t)(uint32_t)idx * 0x9E3779B97F4A7C15ull;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return static_cast<float>((x >> 40) & 0xFFFFFFu) / static_cast<float>(0x1000000u);
}

PlantResult buildTree(const TreeParams& params) {
    PlantResult result;
    using namespace plant_internal;

    // Always build the mature structure. age01 trims an existing skeleton;
    // it does not regenerate one. Same seed -> same skeleton at every age.
    const float age01 = std::clamp(params.age01, 0.0f, 1.0f);
    const float H = params.height;
    const float CR = params.canopyRadius;
    const int attractorCount = std::max(8, params.attractorCount);

    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);
    const float canopyBase = H * 0.45f;
    const float canopyTop = H * 1.05f;
    std::vector<Vec3> attractors;
    attractors.reserve(attractorCount);
    for (int i = 0; i < attractorCount; ++i) {
        Vec3 p;
        do {
            p = {
                uni(rng) * 2.0f - 1.0f,
                uni(rng) * 2.0f - 1.0f,
                uni(rng) * 2.0f - 1.0f
            };
        } while (vdot(p, p) > 1.0f);
        Vec3 a = {
            p.x * CR,
            canopyBase + uni(rng) * (canopyTop - canopyBase),
            p.z * CR
        };
        attractors.push_back(a);
    }

    SpaceColonizationOptions opts;
    opts.attractionRadius = CR * 0.7f;
    opts.killRadius = std::max(params.trunkRadius * 1.5f, CR * 0.08f);
    opts.segmentLength = std::max(0.1f, H * 0.05f);
    opts.maxIterations = 220;
    opts.tropism = {0, 1, 0};
    opts.tropismWeight = 0.35f;

    std::vector<BranchSegment> segs;
    {
        const float trunkTop = canopyBase * 0.9f;
        const float step = opts.segmentLength;
        int trunkSegs = std::max(1, static_cast<int>(trunkTop / step));
        Vec3 prev{0, 0, 0};
        BranchSegment root;
        root.parent = -1;
        root.from = prev;
        root.to = prev;
        root.depth = 0;
        segs.push_back(root);
        for (int i = 1; i <= trunkSegs; ++i) {
            BranchSegment s;
            s.parent = static_cast<int>(segs.size()) - 1;
            s.from = prev;
            s.to = { 0.0f, static_cast<float>(i) * step, 0.0f };
            s.depth = i;
            segs.push_back(s);
            prev = s.to;
        }
    }

    {
        Vec3 trunkTopPos = segs.back().to;
        int trunkTopSeg = static_cast<int>(segs.size()) - 1;
        std::vector<BranchSegment> grown = spaceColonize(
            attractors, { trunkTopPos }, {0, 1, 0}, opts);
        const int offset = static_cast<int>(segs.size()) - 1;
        for (size_t i = 1; i < grown.size(); ++i) {
            BranchSegment s = grown[i];
            if (s.parent == 0) s.parent = trunkTopSeg;
            else if (s.parent > 0) s.parent = s.parent + offset;
            s.depth += static_cast<int>(segs.size()) - 1;
            segs.push_back(s);
        }
    }

    thickenBranches(segs, std::max(0.01f, params.trunkRadius * 0.12f), 2.5f);

    for (auto& s : segs) {
        if (s.parent == -1) s.radius = std::max(s.radius, params.trunkRadius);
    }

    // age01 trim: birth time per segment from depth (monotonic from root
    // outward) plus a small per-segment jitter so siblings don't pop in
    // simultaneously. Heuristic, not biologically rigorous; the contract
    // is determinism + monotonic-from-root order.
    int maxDepth = 1;
    for (const auto& s : segs) if (s.depth > maxDepth) maxDepth = s.depth;
    std::vector<bool> alive(segs.size(), false);
    if (age01 >= 1.0f) {
        std::fill(alive.begin(), alive.end(), true);
    } else {
        for (size_t i = 0; i < segs.size(); ++i) {
            float dt = static_cast<float>(segs[i].depth) / static_cast<float>(maxDepth);
            float jitter = (hash01(params.seed, static_cast<int>(i)) - 0.5f) * 0.05f;
            float t_birth = std::clamp(dt + jitter, 0.0f, 1.0f);
            alive[i] = t_birth <= age01;
        }
        // Drop descendants of any dead ancestor to keep parent indices valid.
        for (size_t i = 0; i < segs.size(); ++i) {
            int p = segs[i].parent;
            if (p >= 0 && !alive[p]) alive[i] = false;
        }
    }

    // Compact alive segments and remap parent indices.
    std::vector<int> remap(segs.size(), -1);
    std::vector<BranchSegment> kept;
    kept.reserve(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        if (!alive[i]) continue;
        BranchSegment s = segs[i];
        s.parent = (s.parent >= 0) ? remap[s.parent] : -1;
        remap[i] = static_cast<int>(kept.size());
        kept.push_back(s);
    }

    result.branchMesh = meshBranches(kept, 8);

    // Leaves: terminal-after-trim segments only.
    std::vector<int> childCount(kept.size(), 0);
    for (const auto& s : kept) if (s.parent >= 0) ++childCount[s.parent];
    // Leaf scale scaled with tree size — large enough that the alpha-cutout
    // atlas survives mipmap LOD averaging at typical viewing distances.
    const float leafScaleBase = std::max(0.15f, H * 0.10f);
    // Dedupe terminal positions: space-colonization can leave multiple
    // near-coincident terminal segments converging on the same attractor.
    // Drop near-duplicates so we don't stack 10 leaf quads on one point.
    // Light-touch dedup: only collapse leaves that are essentially
    // coincident (well under killRadius). Aggressive dedup empties the
    // canopy because terminal density is what makes a tree look full.
    const float dedupeR = std::max(0.001f, opts.killRadius * 0.15f);
    const float dedupeR2 = dedupeR * dedupeR;
    std::vector<Vec3> placedPositions;
    placedPositions.reserve(kept.size());
    for (size_t i = 0; i < kept.size(); ++i) {
        if (childCount[i] != 0) continue;
        if (kept[i].parent == -1) continue;
        // Skip degenerate (zero-length) terminals — they produce a leaf at
        // their parent's tip with no orientation information.
        if (vdist2(kept[i].from, kept[i].to) < 1e-8f) continue;
        bool dup = false;
        for (const Vec3& q : placedPositions) {
            if (vdist2(q, kept[i].to) < dedupeR2) { dup = true; break; }
        }
        if (dup) continue;
        placedPositions.push_back(kept[i].to);
        LeafInstance L;
        L.position = kept[i].to;
        Vec3 fwd = vnormOr(kept[i].to - kept[i].from, {0, 1, 0});
        L.orientation = quatLookDir(fwd);
        // Deterministic per-segment values from a stable hash, not the rng
        // (which is consumed elsewhere). Keeps a leaf's variant stable as
        // the tree grows: the same segment id always picks the same variant.
        float hScale = hash01(params.seed ^ 0xA1B2C3D4u, static_cast<int>(i));
        L.scale = leafScaleBase * (0.7f + hScale * 0.6f);
        uint64_t v = static_cast<uint64_t>(hash01(params.seed ^ 0x55AA55AAu,
                                                  static_cast<int>(i)) * 4.0f);
        L.variantIndex = static_cast<int>(v & 3);
        result.leaves.push_back(L);
    }

    if (!result.branchMesh.empty()) {
        plant_internal::aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    } else {
        result.aabbMin = {0, 0, 0};
        result.aabbMax = {0, H * std::max(0.05f, age01), 0};
    }
    for (const auto& L : result.leaves) {
        Vec3 r{L.scale, L.scale, L.scale};
        plant_internal::updateAabb(result.aabbMin, result.aabbMax, L.position - r);
        plant_internal::updateAabb(result.aabbMin, result.aabbMax, L.position + r);
    }
    return result;
}

} // namespace bromesh
