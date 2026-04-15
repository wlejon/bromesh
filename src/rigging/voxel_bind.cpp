#include "bromesh/rigging/voxel_bind.h"
#include "bromesh/manipulation/skin.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <queue>
#include <vector>

namespace bromesh {

namespace {

using Vec3 = std::array<float, 3>;

Vec3 operator+(Vec3 a, Vec3 b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
Vec3 operator-(Vec3 a, Vec3 b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 operator*(Vec3 a, float s) { return {a[0]*s, a[1]*s, a[2]*s}; }
float dot(Vec3 a, Vec3 b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
float lenSq(Vec3 a) { return dot(a, a); }

// ---- bone world position extraction ---------------------------------------
//
// We need each bone's world-space head + tail in the bind pose. The
// Skeleton only stores inverseBind matrices; world = inverse(inverseBind),
// and +Y axis is the bone direction. Tail = head + Y * 1 unit. We use a
// fixed unit length because bone-length info isn't guaranteed in
// Skeleton (fitSkeleton doesn't store it). For binding that's fine — the
// head position is what matters; we approximate the bone segment from head
// by a short step along +Y.
//
// For bromesh's fitter-produced Skeleton we can reconstruct the exact
// segment from the forward chain (parent head -> child head). We use that
// path since it's available without modifying SkinData/Bone shapes.

void matInverse4(const float* m, float* inv) {
    float a[16]; std::memcpy(a, m, sizeof(a));
    inv[0]  =  a[5]*a[10]*a[15] - a[5]*a[11]*a[14] - a[9]*a[6]*a[15] + a[9]*a[7]*a[14] + a[13]*a[6]*a[11] - a[13]*a[7]*a[10];
    inv[4]  = -a[4]*a[10]*a[15] + a[4]*a[11]*a[14] + a[8]*a[6]*a[15] - a[8]*a[7]*a[14] - a[12]*a[6]*a[11] + a[12]*a[7]*a[10];
    inv[8]  =  a[4]*a[9]*a[15]  - a[4]*a[11]*a[13] - a[8]*a[5]*a[15] + a[8]*a[7]*a[13] + a[12]*a[5]*a[11] - a[12]*a[7]*a[9];
    inv[12] = -a[4]*a[9]*a[14]  + a[4]*a[10]*a[13] + a[8]*a[5]*a[14] - a[8]*a[6]*a[13] - a[12]*a[5]*a[10] + a[12]*a[6]*a[9];
    inv[1]  = -a[1]*a[10]*a[15] + a[1]*a[11]*a[14] + a[9]*a[2]*a[15] - a[9]*a[3]*a[14] - a[13]*a[2]*a[11] + a[13]*a[3]*a[10];
    inv[5]  =  a[0]*a[10]*a[15] - a[0]*a[11]*a[14] - a[8]*a[2]*a[15] + a[8]*a[3]*a[14] + a[12]*a[2]*a[11] - a[12]*a[3]*a[10];
    inv[9]  = -a[0]*a[9]*a[15]  + a[0]*a[11]*a[13] + a[8]*a[1]*a[15] - a[8]*a[3]*a[13] - a[12]*a[1]*a[11] + a[12]*a[3]*a[9];
    inv[13] =  a[0]*a[9]*a[14]  - a[0]*a[10]*a[13] - a[8]*a[1]*a[14] + a[8]*a[2]*a[13] + a[12]*a[1]*a[10] - a[12]*a[2]*a[9];
    inv[2]  =  a[1]*a[6]*a[15]  - a[1]*a[7]*a[14]  - a[5]*a[2]*a[15] + a[5]*a[3]*a[14] + a[13]*a[2]*a[7]  - a[13]*a[3]*a[6];
    inv[6]  = -a[0]*a[6]*a[15]  + a[0]*a[7]*a[14]  + a[4]*a[2]*a[15] - a[4]*a[3]*a[14] - a[12]*a[2]*a[7]  + a[12]*a[3]*a[6];
    inv[10] =  a[0]*a[5]*a[15]  - a[0]*a[7]*a[13]  - a[4]*a[1]*a[15] + a[4]*a[3]*a[13] + a[12]*a[1]*a[7]  - a[12]*a[3]*a[5];
    inv[14] = -a[0]*a[5]*a[14]  + a[0]*a[6]*a[13]  + a[4]*a[1]*a[14] - a[4]*a[2]*a[13] - a[12]*a[1]*a[6]  + a[12]*a[2]*a[5];
    inv[3]  = -a[1]*a[6]*a[11]  + a[1]*a[7]*a[10]  + a[5]*a[2]*a[11] - a[5]*a[3]*a[10] - a[9]*a[2]*a[7]   + a[9]*a[3]*a[6];
    inv[7]  =  a[0]*a[6]*a[11]  - a[0]*a[7]*a[10]  - a[4]*a[2]*a[11] + a[4]*a[3]*a[10] + a[8]*a[2]*a[7]   - a[8]*a[3]*a[6];
    inv[11] = -a[0]*a[5]*a[11]  + a[0]*a[7]*a[9]   + a[4]*a[1]*a[11] - a[4]*a[3]*a[9]  - a[8]*a[1]*a[7]   + a[8]*a[3]*a[5];
    inv[15] =  a[0]*a[5]*a[10]  - a[0]*a[6]*a[9]   - a[4]*a[1]*a[10] + a[4]*a[2]*a[9]  + a[8]*a[1]*a[6]   - a[8]*a[2]*a[5];
    float det = a[0]*inv[0] + a[1]*inv[4] + a[2]*inv[8] + a[3]*inv[12];
    if (std::fabs(det) < 1e-20f) { for (int i=0;i<16;++i) inv[i] = (i%5==0)?1.0f:0.0f; return; }
    float id = 1.0f / det;
    for (int i = 0; i < 16; ++i) inv[i] *= id;
}

// Extract world head position + +Y axis from the bone's inverseBind.
// world = inverse(inverseBind). head = column 3. yAxis = column 1.
void extractHeadY(const float* invBind, Vec3& head, Vec3& yAxis) {
    float w[16];
    matInverse4(invBind, w);
    head  = { w[12], w[13], w[14] };
    yAxis = { w[4],  w[5],  w[6]  };
}

// ---- voxel grid ------------------------------------------------------------

struct Grid {
    int sx = 0, sy = 0, sz = 0;
    float cell = 0.0f;
    Vec3 origin = {0, 0, 0}; // world position of voxel (0,0,0) lower corner
    std::vector<uint8_t> solid; // 1 = inside/boundary, 0 = empty
    int idx(int x, int y, int z) const { return (z * sy + y) * sx + x; }
    bool in(int x, int y, int z) const {
        return x >= 0 && x < sx && y >= 0 && y < sy && z >= 0 && z < sz;
    }
    Vec3 voxelCenter(int x, int y, int z) const {
        return { origin[0] + (x + 0.5f) * cell,
                 origin[1] + (y + 0.5f) * cell,
                 origin[2] + (z + 0.5f) * cell };
    }
    void worldToVoxel(Vec3 p, int& x, int& y, int& z) const {
        x = (int)std::floor((p[0] - origin[0]) / cell);
        y = (int)std::floor((p[1] - origin[1]) / cell);
        z = (int)std::floor((p[2] - origin[2]) / cell);
    }
};

Grid makeGrid(const MeshData& mesh, int maxRes) {
    Grid g;
    if (mesh.empty()) return g;

    float mn[3] = { mesh.positions[0], mesh.positions[1], mesh.positions[2] };
    float mx[3] = { mn[0], mn[1], mn[2] };
    for (size_t i = 3; i < mesh.positions.size(); i += 3) {
        for (int k = 0; k < 3; ++k) {
            mn[k] = std::min(mn[k], mesh.positions[i + k]);
            mx[k] = std::max(mx[k], mesh.positions[i + k]);
        }
    }
    float ext[3] = { mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2] };
    float maxExt = std::max({ext[0], ext[1], ext[2]});
    if (maxExt < 1e-12f) return g;

    // Padding: 2 cells each side for flood-fill exterior guarantee.
    g.cell = maxExt / (float)std::max(8, maxRes - 4);
    int pad = 2;
    g.sx = (int)std::ceil(ext[0] / g.cell) + 2 * pad;
    g.sy = (int)std::ceil(ext[1] / g.cell) + 2 * pad;
    g.sz = (int)std::ceil(ext[2] / g.cell) + 2 * pad;
    g.origin = { mn[0] - pad * g.cell, mn[1] - pad * g.cell, mn[2] - pad * g.cell };
    g.solid.assign((size_t)g.sx * g.sy * g.sz, 0);
    return g;
}

// Surface rasterization by dense barycentric sampling. Simple & robust for
// any topology including non-manifold. Slight over-rasterization is fine —
// the flood-fill downstream handles it.
void rasterizeSurface(const MeshData& mesh, Grid& g) {
    const auto& pos = mesh.positions;
    const auto& idx = mesh.indices;
    for (size_t t = 0; t + 2 < idx.size(); t += 3) {
        Vec3 a = { pos[idx[t]*3], pos[idx[t]*3+1], pos[idx[t]*3+2] };
        Vec3 b = { pos[idx[t+1]*3], pos[idx[t+1]*3+1], pos[idx[t+1]*3+2] };
        Vec3 c = { pos[idx[t+2]*3], pos[idx[t+2]*3+1], pos[idx[t+2]*3+2] };
        float e1 = std::sqrt(lenSq(b - a));
        float e2 = std::sqrt(lenSq(c - a));
        float e3 = std::sqrt(lenSq(c - b));
        float eMax = std::max({e1, e2, e3});
        int N = std::max(2, (int)std::ceil(eMax / g.cell * 2.0f) + 1);
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N - i; ++j) {
                float u = (float)i / (float)N;
                float v = (float)j / (float)N;
                float w = 1.0f - u - v;
                Vec3 p = a * w + b * u + c * v;
                int x, y, z; g.worldToVoxel(p, x, y, z);
                if (g.in(x, y, z)) g.solid[g.idx(x, y, z)] = 1;
            }
        }
    }
}

// Flood-fill exterior from corner, then invert to get solid (interior +
// boundary). Boundary voxels stay 1; previously unreached air becomes 1
// only if it's enclosed (we flip exterior 0 from the fill).
void fillInterior(Grid& g) {
    std::vector<uint8_t> ext(g.solid.size(), 0);
    std::queue<int> q;
    auto push = [&](int x, int y, int z) {
        if (!g.in(x, y, z)) return;
        int i = g.idx(x, y, z);
        if (g.solid[i] || ext[i]) return;
        ext[i] = 1;
        q.push(i);
    };
    // Seed: the corner voxel. (All six faces would be safer but overkill
    // given the 2-cell padding.)
    push(0, 0, 0);

    const int dx[6] = { 1, -1, 0,  0, 0,  0 };
    const int dy[6] = { 0,  0, 1, -1, 0,  0 };
    const int dz[6] = { 0,  0, 0,  0, 1, -1 };
    while (!q.empty()) {
        int i = q.front(); q.pop();
        int z = i / (g.sx * g.sy);
        int rem = i - z * g.sx * g.sy;
        int y = rem / g.sx;
        int x = rem - y * g.sx;
        for (int k = 0; k < 6; ++k) push(x + dx[k], y + dy[k], z + dz[k]);
    }
    // Anything not exterior = solid.
    for (size_t i = 0; i < g.solid.size(); ++i) {
        if (!ext[i]) g.solid[i] = 1;
    }
}

// ---- per-bone geodesic BFS -------------------------------------------------

// Rasterize bone segment (head -> tail) into seed voxels on the grid.
// Uses Bresenham-like sampling along the segment at cell/2 spacing.
std::vector<int> seedVoxels(const Grid& g, Vec3 head, Vec3 tail) {
    std::vector<int> seeds;
    Vec3 d = tail - head;
    float L = std::sqrt(lenSq(d));
    int steps = std::max(2, (int)std::ceil(L / (g.cell * 0.5f)));
    for (int s = 0; s <= steps; ++s) {
        float t = (float)s / (float)steps;
        Vec3 p = head + d * t;
        int x, y, z; g.worldToVoxel(p, x, y, z);
        if (g.in(x, y, z) && g.solid[g.idx(x, y, z)]) {
            seeds.push_back(g.idx(x, y, z));
        }
    }
    // If no seed landed inside the solid volume (bone head outside mesh),
    // fall back to the nearest solid voxel to the head — ensures every
    // bone contributes something.
    if (seeds.empty()) {
        int bx, by, bz; g.worldToVoxel(head, bx, by, bz);
        int best = -1; int bestDistSq = std::numeric_limits<int>::max();
        for (int z = 0; z < g.sz; ++z) for (int y = 0; y < g.sy; ++y) for (int x = 0; x < g.sx; ++x) {
            if (!g.solid[g.idx(x, y, z)]) continue;
            int dxx = x - bx, dyy = y - by, dzz = z - bz;
            int d2 = dxx*dxx + dyy*dyy + dzz*dzz;
            if (d2 < bestDistSq) { bestDistSq = d2; best = g.idx(x, y, z); }
        }
        if (best >= 0) seeds.push_back(best);
    }
    return seeds;
}

// 6-neighbor BFS. Distances are in voxel hops (uint16_t; capped at 65535).
// Unreached voxels keep 65535 sentinel.
void bfsFromSeeds(const Grid& g,
                  const std::vector<int>& seeds,
                  std::vector<uint16_t>& dist) {
    constexpr uint16_t INF = std::numeric_limits<uint16_t>::max();
    dist.assign(g.solid.size(), INF);
    std::queue<int> q;
    for (int s : seeds) { dist[s] = 0; q.push(s); }

    const int dx[6] = { 1, -1, 0,  0, 0,  0 };
    const int dy[6] = { 0,  0, 1, -1, 0,  0 };
    const int dz[6] = { 0,  0, 0,  0, 1, -1 };
    while (!q.empty()) {
        int i = q.front(); q.pop();
        int z = i / (g.sx * g.sy);
        int rem = i - z * g.sx * g.sy;
        int y = rem / g.sx;
        int x = rem - y * g.sx;
        uint16_t nd = dist[i] + 1;
        if (nd == INF) continue;
        for (int k = 0; k < 6; ++k) {
            int nx = x + dx[k], ny = y + dy[k], nz = z + dz[k];
            if (!g.in(nx, ny, nz)) continue;
            int ni = g.idx(nx, ny, nz);
            if (!g.solid[ni]) continue;
            if (dist[ni] <= nd) continue;
            dist[ni] = nd;
            q.push(ni);
        }
    }
}

// ---- per-voxel top-K accumulator ------------------------------------------

struct TopK {
    // Parallel arrays, length maxInfluences, sorted descending by weight.
    std::vector<float>   w;
    std::vector<int32_t> b;
    int K;
    TopK(int k) : w(k, 0.0f), b(k, -1), K(k) {}
    void insert(float weight, int bone) {
        if (weight <= w[K - 1]) return; // too small to fit
        w[K - 1] = weight;
        b[K - 1] = bone;
        // Bubble up.
        for (int i = K - 1; i > 0 && w[i] > w[i - 1]; --i) {
            std::swap(w[i], w[i - 1]);
            std::swap(b[i], b[i - 1]);
        }
    }
};

} // namespace

SkinData voxelBindWeights(const MeshData& mesh,
                          const Skeleton& skeleton,
                          const VoxelBindOptions& opts) {
    SkinData out;
    if (mesh.empty() || skeleton.bones.empty()) return out;

    const int K = std::max(1, opts.maxInfluences);
    const size_t nVerts = mesh.vertexCount();
    const size_t nBones = skeleton.bones.size();

    // Build grid + rasterize + fill.
    Grid g = makeGrid(mesh, opts.maxResolution);
    if (g.solid.empty()) return out;
    rasterizeSurface(mesh, g);
    fillInterior(g);
    const size_t V = g.solid.size();

    // Extract per-bone world head. For segment we use head + 1-cell along +Y,
    // or parent->this if we have a parent (more accurate bone direction).
    std::vector<Vec3> headW(nBones), tailW(nBones);
    for (size_t i = 0; i < nBones; ++i) {
        Vec3 h, y;
        extractHeadY(skeleton.bones[i].inverseBind, h, y);
        headW[i] = h;
    }
    for (size_t i = 0; i < nBones; ++i) {
        // Tail: use first child's head if available; else head + yAxis * cell * 4.
        int firstChild = -1;
        for (size_t j = 0; j < nBones; ++j) {
            if (skeleton.bones[j].parent == (int)i) { firstChild = (int)j; break; }
        }
        if (firstChild >= 0) {
            tailW[i] = headW[firstChild];
        } else {
            Vec3 h, y;
            extractHeadY(skeleton.bones[i].inverseBind, h, y);
            tailW[i] = h + y * (g.cell * 4.0f);
        }
    }

    // Per-voxel top-K.
    std::vector<float> topW((size_t)V * K, 0.0f);
    std::vector<int32_t> topB((size_t)V * K, -1);

    std::vector<uint16_t> dist;
    constexpr uint16_t INF = std::numeric_limits<uint16_t>::max();
    const float eps = 1e-6f;
    const float pow2 = opts.falloffPower;

    for (size_t b = 0; b < nBones; ++b) {
        auto seeds = seedVoxels(g, headW[b], tailW[b]);
        if (seeds.empty()) continue;
        bfsFromSeeds(g, seeds, dist);
        for (size_t v = 0; v < V; ++v) {
            if (!g.solid[v]) continue;
            if (dist[v] == INF) continue;
            float d = (float)dist[v];
            float w = 1.0f / (std::pow(d + 1.0f, pow2) + eps);
            // Insert into voxel's top-K.
            float* vw = &topW[v * K];
            int32_t* vb = &topB[v * K];
            if (w <= vw[K - 1]) continue;
            vw[K - 1] = w;
            vb[K - 1] = (int32_t)b;
            for (int i = K - 1; i > 0 && vw[i] > vw[i - 1]; --i) {
                std::swap(vw[i], vw[i - 1]);
                std::swap(vb[i], vb[i - 1]);
            }
        }
    }

    // Normalize per-voxel weights so they sum to 1.
    for (size_t v = 0; v < V; ++v) {
        if (!g.solid[v]) continue;
        float sum = 0.0f;
        for (int k = 0; k < K; ++k) sum += topW[v * K + k];
        if (sum > 0.0f) {
            float inv = 1.0f / sum;
            for (int k = 0; k < K; ++k) topW[v * K + k] *= inv;
        }
    }

    // Per-vertex sampling. For each vertex, find nearest solid voxel
    // (start at its exact voxel; expand to 3x3x3 neighborhood if empty).
    // Then union the top-K contributions from the sampled voxels weighted
    // by trilinear distance, keep final top-K, normalize.
    out.boneWeights.assign(nVerts * (size_t)K, 0.0f);
    out.boneIndices.assign(nVerts * (size_t)K, 0);
    out.boneCount = nBones;
    out.inverseBindMatrices.assign(nBones * 16, 0.0f);
    for (size_t b = 0; b < nBones; ++b) {
        std::memcpy(&out.inverseBindMatrices[b * 16],
                    skeleton.bones[b].inverseBind, 16 * sizeof(float));
    }

    for (size_t v = 0; v < nVerts; ++v) {
        Vec3 p = { mesh.positions[v*3], mesh.positions[v*3+1], mesh.positions[v*3+2] };
        int vx, vy, vz; g.worldToVoxel(p, vx, vy, vz);

        // Collect contributions from the vertex voxel plus 26 neighbors.
        // Weight each by 1/(1+cellDist).
        std::vector<float> acc(nBones, 0.0f);
        for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            int nx = vx + dx, ny = vy + dy, nz = vz + dz;
            if (!g.in(nx, ny, nz)) continue;
            int vi = g.idx(nx, ny, nz);
            if (!g.solid[vi]) continue;
            float weight = 1.0f / (1.0f + (float)(dx*dx + dy*dy + dz*dz));
            for (int k = 0; k < K; ++k) {
                int32_t b = topB[vi * K + k];
                if (b < 0) continue;
                acc[b] += weight * topW[vi * K + k];
            }
        }

        // If we got nothing (vertex far outside any solid neighborhood),
        // fall back to nearest solid voxel via expanding search.
        float accSum = 0.0f;
        for (float a : acc) accSum += a;
        if (accSum <= 0.0f) {
            // Linear scan for the closest solid voxel. Slow path but rarely
            // taken and bounded to outlier vertices.
            int bestVI = -1; float bestD2 = std::numeric_limits<float>::max();
            for (int z = 0; z < g.sz; ++z) for (int y = 0; y < g.sy; ++y) for (int x = 0; x < g.sx; ++x) {
                int vi = g.idx(x, y, z);
                if (!g.solid[vi]) continue;
                Vec3 c = g.voxelCenter(x, y, z);
                float d2 = lenSq(c - p);
                if (d2 < bestD2) { bestD2 = d2; bestVI = vi; }
            }
            if (bestVI >= 0) {
                for (int k = 0; k < K; ++k) {
                    int32_t bi = topB[bestVI * K + k];
                    if (bi < 0) continue;
                    acc[bi] += topW[bestVI * K + k];
                }
            }
        }

        // Top-K select from acc, prune below minWeight, normalize.
        TopK tk(K);
        for (size_t b = 0; b < nBones; ++b) tk.insert(acc[b], (int)b);

        // Prune below minWeight.
        float sum = 0.0f;
        for (int k = 0; k < K; ++k) {
            if (tk.w[k] < opts.minWeight) { tk.w[k] = 0.0f; tk.b[k] = -1; }
            sum += tk.w[k];
        }
        if (sum <= 0.0f) {
            // Everything pruned: fall back to the single best bone at weight 1.
            // Pick the bone whose world head is closest to the vertex.
            int bestBone = 0; float bestD2 = std::numeric_limits<float>::max();
            for (size_t bi = 0; bi < nBones; ++bi) {
                float d2 = lenSq(headW[bi] - p);
                if (d2 < bestD2) { bestD2 = d2; bestBone = (int)bi; }
            }
            out.boneWeights[v * K + 0] = 1.0f;
            out.boneIndices[v * K + 0] = (uint32_t)bestBone;
            continue;
        }
        float inv = 1.0f / sum;
        for (int k = 0; k < K; ++k) {
            float w = tk.w[k] * inv;
            out.boneWeights[v * K + k] = w;
            out.boneIndices[v * K + k] = (tk.b[k] < 0) ? 0u : (uint32_t)tk.b[k];
        }
    }

    // Final belt-and-braces normalization from the existing helper.
    normalizeWeights(out);
    return out;
}

} // namespace bromesh
