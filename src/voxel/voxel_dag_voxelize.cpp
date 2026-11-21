#include "bromesh/voxel/voxel_dag.h"
#include "bromesh/analysis/bbox.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <deque>
#include <unordered_map>
#include <vector>

namespace bromesh {

namespace {

// 64-bit brick coordinate key (21 bits per axis in [0..2^20)).
inline uint64_t packBrickKey(int bx, int by, int bz) noexcept {
    return (static_cast<uint64_t>(bx & 0xFFFFF) << 40u) |
           (static_cast<uint64_t>(by & 0xFFFFF) << 20u) |
           static_cast<uint64_t>(bz & 0xFFFFF);
}

inline int roundDepthForResolution(int res) noexcept {
    int r = std::clamp(res, 8, 256);
    int depth = 1;
    while ((4 << depth) < r && depth < 6) {
        ++depth;
    }
    return depth;
}

struct AccumVoxel {
    uint8_t state = 0; // bit 0 = Solid, bit 1 = Carve
    float r = 0.0f, g = 0.0f, b = 0.0f, a = 0.0f;
    float nx = 0.0f, ny = 0.0f, nz = 0.0f;
    float roughness = 0.8f;
    float emissive = 0.0f;
    float weight = 0.0f;
};

struct SparseBrick {
    std::array<AccumVoxel, 64> v{};
};

struct PackedBrick {
    uint64_t key;
    uint32_t solidLo, solidHi, carveLo, carveHi;
    uint32_t attrFirst, attrCount;
    int32_t same;    // an earlier brick it is identical to, or -1
};

struct PackedBricks {
    int depth = 0;
    std::vector<PackedBrick> bricks;
    std::vector<DagVoxelAttr> attrs;
};

// The bricks a voxelization touches, found through a dense index over the
// brick grid ((res / 4)^3 slots, at most 64^3), not a hash map: a mesh adds a
// voxel several times over and each add is one array read.
class SparseBrickBuilder {
public:
    explicit SparseBrickBuilder(int depth)
        : depth_(depth), res_(4 << depth), n_(1 << depth),
          index_(static_cast<size_t>(n_) * n_ * n_, -1) {}

    void addVoxel(int gx, int gy, int gz, VoxelState state,
                  float r, float g, float b, float a,
                  float nx, float ny, float nz, float roughness,
                  float emissive = 0.0f) {
        if (gx < 0 || gy < 0 || gz < 0 || gx >= res_ || gy >= res_ || gz >= res_) return;
        if (state == VoxelState::Empty) return;

        SparseBrick& sb = brickAt(gx >> 2, gy >> 2, gz >> 2);
        const int bit = (gx & 3) | ((gy & 3) << 2) | ((gz & 3) << 4);
        AccumVoxel& av = sb.v[bit];
        if (state == VoxelState::Solid) {
            av.state |= 1u;
            av.r += r;
            av.g += g;
            av.b += b;
            av.a += a;
            av.nx += nx;
            av.ny += ny;
            av.nz += nz;
            av.roughness = (av.weight > 0.0f) ? ((av.roughness * av.weight + roughness) / (av.weight + 1.0f)) : roughness;
            av.emissive = std::max(av.emissive, emissive);
            av.weight += 1.0f;
        } else if (state == VoxelState::Carve) {
            av.state |= 2u;
        }
    }

    bool hasVoxelAt(int gx, int gy, int gz, uint8_t bitMask = 1u) const {
        if (gx < 0 || gy < 0 || gz < 0 || gx >= res_ || gy >= res_ || gz >= res_) return false;
        const int32_t slot = index_[slotOf(gx >> 2, gy >> 2, gz >> 2)];
        if (slot < 0) return false;
        const int bit = (gx & 3) | ((gy & 3) << 2) | ((gz & 3) << 4);
        return (bricks_[static_cast<size_t>(slot)].v[bit].state & bitMask) != 0u;
    }

    int resolution() const noexcept { return res_; }
    int depth() const noexcept { return depth_; }
    bool empty() const noexcept { return bricks_.empty(); }

    // Make solid every cell the surface encloses: those no 6-connected path of
    // empty cells joins to the grid's boundary. They take one colour, the
    // surface voxels' mean, so the full bricks inside are one brick.
    void fillInside(float roughness, float emissive) {
        const size_t R = static_cast<size_t>(res_);
        std::vector<uint8_t> cell(R * R * R, 0);    // 1: surface, 2: outside
        double sum[4] = {0.0, 0.0, 0.0, 0.0};
        double count = 0.0;
        for (size_t i = 0; i < bricks_.size(); ++i) {
            const uint64_t key = keys_[i];
            const size_t bx = static_cast<size_t>((key >> 40u) & 0xFFFFFu);
            const size_t by = static_cast<size_t>((key >> 20u) & 0xFFFFFu);
            const size_t bz = static_cast<size_t>(key & 0xFFFFFu);
            for (int bit = 0; bit < 64; ++bit) {
                const AccumVoxel& av = bricks_[i].v[bit];
                if (!(av.state & 1u)) continue;
                const size_t x = bx * 4 + (bit & 3), y = by * 4 + ((bit >> 2) & 3), z = bz * 4 + (bit >> 4);
                cell[(z * R + y) * R + x] = 1;
                const float w = std::max(av.weight, 1.0f);
                sum[0] += av.r / w; sum[1] += av.g / w; sum[2] += av.b / w; sum[3] += av.a / w;
                count += 1.0;
            }
        }
        if (count == 0.0) return;
        std::vector<uint32_t> stack;
        stack.reserve(R * R * 6);
        auto push = [&](size_t x, size_t y, size_t z) {
            const size_t c = (z * R + y) * R + x;
            if (cell[c] == 0) { cell[c] = 2; stack.push_back(static_cast<uint32_t>(c)); }
        };
        for (size_t a = 0; a < R; ++a) {
            for (size_t b = 0; b < R; ++b) {
                push(a, b, 0); push(a, b, R - 1);
                push(a, 0, b); push(a, R - 1, b);
                push(0, a, b); push(R - 1, a, b);
            }
        }
        while (!stack.empty()) {
            const size_t c = stack.back();
            stack.pop_back();
            const size_t x = c % R, y = (c / R) % R, z = c / (R * R);
            if (x > 0) push(x - 1, y, z);
            if (x + 1 < R) push(x + 1, y, z);
            if (y > 0) push(x, y - 1, z);
            if (y + 1 < R) push(x, y + 1, z);
            if (z > 0) push(x, y, z - 1);
            if (z + 1 < R) push(x, y, z + 1);
        }
        const float r = static_cast<float>(sum[0] / count), g = static_cast<float>(sum[1] / count);
        const float b = static_cast<float>(sum[2] / count), a = static_cast<float>(sum[3] / count);
        for (size_t z = 0; z < R; ++z) {
            for (size_t y = 0; y < R; ++y) {
                for (size_t x = 0; x < R; ++x) {
                    if (cell[(z * R + y) * R + x] != 0) continue;
                    addVoxel(static_cast<int>(x), static_cast<int>(y), static_cast<int>(z), VoxelState::Solid,
                             r, g, b, a, 0.0f, 0.0f, 0.0f, roughness, emissive);
                }
            }
        }
    }

    uint32_t commit(VoxelDagPool& pool) const { return commitPacked(pool, pack()); }

    // The bricks as the pool takes them: masks and packed attributes, no pool
    // touched, so any thread. A brick whose 64 voxels are solid and alike (a
    // filled inside) is packed once; the rest name it (`same`).
    PackedBricks pack() const {
        PackedBricks out;
        out.depth = depth_;
        out.bricks.reserve(bricks_.size());
        std::unordered_map<uint64_t, int32_t> uniform;
        for (size_t i = 0; i < bricks_.size(); ++i) {
            const SparseBrick& sb = bricks_[i];
            const AccumVoxel& v0 = sb.v[0];
            bool alike = v0.state == 1u;
            for (int bit = 1; alike && bit < 64; ++bit) {
                const AccumVoxel& av = sb.v[bit];
                alike = av.state == v0.state && av.r == v0.r && av.g == v0.g && av.b == v0.b && av.a == v0.a &&
                        av.nx == v0.nx && av.ny == v0.ny && av.nz == v0.nz && av.weight == v0.weight &&
                        av.roughness == v0.roughness && av.emissive == v0.emissive;
            }
            PackedBrick pb{};
            pb.key = keys_[i];
            pb.same = -1;
            uint64_t uniformKey = 0;
            if (alike) {
                const float w = std::max(v0.weight, 1.0f);
                float nx = v0.nx, ny = v0.ny, nz = v0.nz;
                if (nx * nx + ny * ny + nz * nz < 1e-6f) ny = 1.0f;
                uniformKey = (static_cast<uint64_t>(dag_pack::packRgba8(v0.r / w, v0.g / w, v0.b / w, std::clamp(v0.a / w, 0.0f, 1.0f))) << 32u) |
                             dag_pack::packNormalMat(nx, ny, nz, v0.roughness, v0.emissive);
                if (const auto it = uniform.find(uniformKey); it != uniform.end()) {
                    pb.same = it->second;
                    out.bricks.push_back(pb);
                    continue;
                }
                uniform.emplace(uniformKey, static_cast<int32_t>(out.bricks.size()));
            }
            uint64_t solidBits = 0ull, carveBits = 0ull;
            pb.attrFirst = static_cast<uint32_t>(out.attrs.size());
            for (int bit = 0; bit < 64; ++bit) {
                const AccumVoxel& av = sb.v[bit];
                if (av.state & 2u) carveBits |= (1ull << bit);
                if (av.state & 1u) {
                    solidBits |= (1ull << bit);
                    const float w = std::max(av.weight, 1.0f);
                    const float r = av.r / w;
                    const float g = av.g / w;
                    const float b = av.b / w;
                    const float a = std::clamp(av.a / w, 0.0f, 1.0f);
                    float nx = av.nx, ny = av.ny, nz = av.nz;
                    if (nx * nx + ny * ny + nz * nz < 1e-6f) {
                        ny = 1.0f;
                    }
                    out.attrs.push_back(DagVoxelAttr{
                        dag_pack::packRgba8(r, g, b, a),
                        dag_pack::packNormalMat(nx, ny, nz, av.roughness, av.emissive),
                    });
                }
            }
            pb.solidLo = static_cast<uint32_t>(solidBits);
            pb.solidHi = static_cast<uint32_t>(solidBits >> 32u);
            pb.carveLo = static_cast<uint32_t>(carveBits);
            pb.carveHi = static_cast<uint32_t>(carveBits >> 32u);
            pb.attrCount = static_cast<uint32_t>(out.attrs.size()) - pb.attrFirst;
            out.bricks.push_back(pb);
        }
        return out;
    }

    // Intern packed bricks (2^depth of them along each axis), then the nodes
    // above them level by level; returns the root.
    static uint32_t commitPacked(VoxelDagPool& pool, const PackedBricks& packed) {
        if (packed.bricks.empty()) return kDagNullNode;
        const int depth = packed.depth;
        std::unordered_map<uint64_t, uint32_t> curNodes;
        curNodes.reserve(packed.bricks.size());
        std::vector<uint32_t> internedAt(packed.bricks.size(), kDagNullNode);
        std::vector<DagVoxelAttr> solidAttrs;
        for (size_t i = 0; i < packed.bricks.size(); ++i) {
            const PackedBrick& pb = packed.bricks[i];
            uint32_t brickIdx = kDagNullNode;
            if (pb.same >= 0) {
                brickIdx = internedAt[static_cast<size_t>(pb.same)];
            } else {
                solidAttrs.assign(packed.attrs.begin() + pb.attrFirst, packed.attrs.begin() + pb.attrFirst + pb.attrCount);
                brickIdx = pool.internBrick(pb.solidLo, pb.solidHi, pb.carveLo, pb.carveHi, solidAttrs);
            }
            internedAt[i] = brickIdx;
            if (brickIdx != kDagNullNode) curNodes.emplace(pb.key, brickIdx);
        }

        if (curNodes.empty()) return kDagNullNode;

        // Bubble up bottom-up through `depth` internal levels (bricks are leaves at level == depth).
        for (int level = depth - 1; level >= 0; --level) {
            const bool childrenAreLeaves = (level == depth - 1);
            struct ParentBucket {
                uint32_t refs[8];
                ParentBucket() {
                    for (uint32_t& r : refs) r = kDagNullNode;
                }
            };
            std::unordered_map<uint64_t, ParentBucket> parents;
            for (const auto& [ckey, childIdx] : curNodes) {
                const int cx = static_cast<int>((ckey >> 40u) & 0xFFFFFu);
                const int cy = static_cast<int>((ckey >> 20u) & 0xFFFFFu);
                const int cz = static_cast<int>(ckey & 0xFFFFFu);
                const int px = cx >> 1;
                const int py = cy >> 1;
                const int pz = cz >> 1;
                const int oct = (cx & 1) | ((cy & 1) << 1) | ((cz & 1) << 2);
                parents[packBrickKey(px, py, pz)].refs[oct] = childIdx;
            }

            std::unordered_map<uint64_t, uint32_t> nextNodes;
            nextNodes.reserve(parents.size());
            for (const auto& [pkey, bucket] : parents) {
                uint8_t leaves[8] = {}, solids[8] = {}, carves[8] = {}, uniCarve[8] = {};
                for (int o = 0; o < 8; ++o) {
                    const uint32_t cref = bucket.refs[o];
                    if (cref == kDagNullNode) continue;
                    if (childrenAreLeaves) {
                        const DagBrick& b = pool.bricks()[cref];
                        leaves[o] = 1u;
                        solids[o] = ((b.solidLo | b.solidHi) != 0u) ? 1u : 0u;
                        carves[o] = ((b.carveLo | b.carveHi) != 0u) ? 1u : 0u;
                        uniCarve[o] = ((b.flags & kBrickFlagUniformCarve) != 0u) ? 1u : 0u;
                    } else {
                        const DagNode& n = pool.nodes()[cref];
                        leaves[o] = 0u;
                        solids[o] = (dag_pack::solidMaskOf(n) != 0u) ? 1u : 0u;
                        carves[o] = (dag_pack::carveMaskOf(n) != 0u) ? 1u : 0u;
                        uniCarve[o] = ((dag_pack::flagsOf(n) & kDagFlagUniformCarve) != 0u) ? 1u : 0u;
                    }
                }
                const uint32_t nodeIdx = pool.internNode(bucket.refs, leaves, solids, carves, uniCarve);
                if (nodeIdx != kDagNullNode) {
                    nextNodes.emplace(pkey, nodeIdx);
                }
            }
            curNodes = std::move(nextNodes);
        }

        auto rootIt = curNodes.find(packBrickKey(0, 0, 0));
        return (rootIt != curNodes.end()) ? rootIt->second : kDagNullNode;
    }

private:
    size_t slotOf(int bx, int by, int bz) const noexcept {
        return (static_cast<size_t>(bz) * n_ + by) * n_ + bx;
    }

    SparseBrick& brickAt(int bx, int by, int bz) {
        int32_t& slot = index_[slotOf(bx, by, bz)];
        if (slot < 0) {
            slot = static_cast<int32_t>(bricks_.size());
            bricks_.emplace_back();
            keys_.push_back(packBrickKey(bx, by, bz));
        }
        return bricks_[static_cast<size_t>(slot)];
    }

    int depth_;
    int res_;
    int n_;
    std::vector<int32_t> index_;
    std::deque<SparseBrick> bricks_;    // never moves an element as it grows
    std::vector<uint64_t> keys_;
};

// Conservative triangle vs axis-aligned voxel-cell overlap (Akenine-Möller SAT).
bool planeBoxOverlap(const float normal[3], const float vert[3], const float halfSize[3]) {
    float vmin[3], vmax[3];
    for (int q = 0; q < 3; ++q) {
        const float v = vert[q];
        if (normal[q] > 0.0f) {
            vmin[q] = -halfSize[q] - v;
            vmax[q] = halfSize[q] - v;
        } else {
            vmin[q] = halfSize[q] - v;
            vmax[q] = -halfSize[q] - v;
        }
    }
    const float dmin = normal[0] * vmin[0] + normal[1] * vmin[1] + normal[2] * vmin[2];
    if (dmin > 0.0f) return false;
    const float dmax = normal[0] * vmax[0] + normal[1] * vmax[1] + normal[2] * vmax[2];
    return dmax >= 0.0f;
}

bool triBoxOverlap(const float boxCenter[3], const float halfSize[3],
                   const float p0[3], const float p1[3], const float p2[3]) {
    const float v0[3] = {p0[0] - boxCenter[0], p0[1] - boxCenter[1], p0[2] - boxCenter[2]};
    const float v1[3] = {p1[0] - boxCenter[0], p1[1] - boxCenter[1], p1[2] - boxCenter[2]};
    const float v2[3] = {p2[0] - boxCenter[0], p2[1] - boxCenter[1], p2[2] - boxCenter[2]};

    const float e0[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    const float e1[3] = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
    const float e2[3] = {v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]};

    auto axisTest = [&](float a, float b, float fa, float fb, int c0, int c1) -> bool {
        const float d0 = a * v0[c0] - b * v0[c1];
        const float d1 = a * v1[c0] - b * v1[c1];
        const float d2 = a * v2[c0] - b * v2[c1];
        const float rad = fa * halfSize[c0] + fb * halfSize[c1];
        return std::min({d0, d1, d2}) <= rad && std::max({d0, d1, d2}) >= -rad;
    };

    const float* edges[3] = {e0, e1, e2};
    for (int ei = 0; ei < 3; ++ei) {
        const float* e = edges[ei];
        const float fex = std::abs(e[0]), fey = std::abs(e[1]), fez = std::abs(e[2]);
        // X-cross: (1, 0, 0) x e = (0, -ez, ey) -> ez*v_y - ey*v_z -> axisTest(ez, ey, fez, fey, 1, 2)
        if (!axisTest(e[2], e[1], fez, fey, 1, 2)) return false;
        // Y-cross: (0, 1, 0) x e = (ez, 0, -ex) -> ez*v_x - ex*v_z -> axisTest(ez, ex, fez, fex, 0, 2)
        if (!axisTest(e[2], e[0], fez, fex, 0, 2)) return false;
        // Z-cross: (0, 0, 1) x e = (-ey, ex, 0) -> ex*v_y - ey*v_x -> axisTest(ex, ey, fex, fey, 1, 0)
        if (!axisTest(e[0], e[1], fex, fey, 1, 0)) return false;
    }

    for (int q = 0; q < 3; ++q) {
        const float lo = std::min({v0[q], v1[q], v2[q]});
        const float hi = std::max({v0[q], v1[q], v2[q]});
        if (lo > halfSize[q] || hi < -halfSize[q]) return false;
    }

    const float normal[3] = {
        e0[1] * e1[2] - e0[2] * e1[1],
        e0[2] * e1[0] - e0[0] * e1[2],
        e0[0] * e1[1] - e0[1] * e1[0],
    };
    return planeBoxOverlap(normal, v0, halfSize);
}

void closestBarycentricOnTriangle(const float p[3], const float a[3], const float b[3], const float c[3],
                                 float& u, float& v, float& w) {
    const float ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const float ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const float ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
    const float d00 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    const float d01 = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];
    const float d11 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
    const float d20 = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2];
    const float d21 = ap[0] * ac[0] + ap[1] * ac[1] + ap[2] * ac[2];
    const float denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < 1e-12f) {
        u = 1.0f / 3.0f;
        v = 1.0f / 3.0f;
        w = 1.0f / 3.0f;
        return;
    }
    v = std::clamp((d11 * d20 - d01 * d21) / denom, 0.0f, 1.0f);
    w = std::clamp((d00 * d21 - d01 * d20) / denom, 0.0f, 1.0f);
    if (v + w > 1.0f) {
        const float s = 1.0f / (v + w);
        v *= s;
        w *= s;
    }
    u = 1.0f - v - w;
}

uint32_t readU32LE(const uint8_t* p) noexcept {
    return static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8u) |
           (static_cast<uint32_t>(p[2]) << 16u) | (static_cast<uint32_t>(p[3]) << 24u);
}

int32_t readI32LE(const uint8_t* p) noexcept {
    const uint32_t u = readU32LE(p);
    int32_t s = 0;
    std::memcpy(&s, &u, sizeof(s));
    return s;
}

} // namespace

uint32_t buildDagFromChunk(
    VoxelDagPool& pool,
    const VoxelChunk& chunk,
    const float* palette,
    int paletteCount,
    int* outDepth) {
    const int sx = chunk.sizeX(), sy = chunk.sizeY(), sz = chunk.sizeZ();
    const int maxDim = std::max({sx, sy, sz, 1});
    const int depth = roundDepthForResolution(maxDim);
    const int res = 4 << depth;
    const int offX = (res - sx) / 2;
    const int offY = (res - sy) / 2;
    const int offZ = (res - sz) / 2;

    SparseBrickBuilder builder(depth);
    auto isSolid = [&](int x, int y, int z) -> bool {
        return chunk.getVoxel(x, y, z) != 0;
    };

    for (int z = 0; z < sz; ++z) {
        for (int y = 0; y < sy; ++y) {
            for (int x = 0; x < sx; ++x) {
                const uint8_t mat = chunk.getVoxel(x, y, z);
                if (mat == 0) continue;

                float r = 0.72f, g = 0.72f, b = 0.72f, a = 1.0f;
                if (palette && mat < paletteCount) {
                    r = palette[mat * 4 + 0];
                    g = palette[mat * 4 + 1];
                    b = palette[mat * 4 + 2];
                    a = palette[mat * 4 + 3];
                } else {
                    r = 0.35f + 0.6f * static_cast<float>((mat * 73) % 255) / 255.0f;
                    g = 0.35f + 0.6f * static_cast<float>((mat * 151) % 255) / 255.0f;
                    b = 0.35f + 0.6f * static_cast<float>((mat * 211) % 255) / 255.0f;
                }

                const float nx = (isSolid(x - 1, y, z) ? 1.0f : 0.0f) - (isSolid(x + 1, y, z) ? 1.0f : 0.0f);
                const float ny = (isSolid(x, y - 1, z) ? 1.0f : 0.0f) - (isSolid(x, y + 1, z) ? 1.0f : 0.0f);
                const float nz = (isSolid(x, y, z - 1) ? 1.0f : 0.0f) - (isSolid(x, y, z + 1) ? 1.0f : 0.0f);
                builder.addVoxel(x + offX, y + offY, z + offZ, VoxelState::Solid, r, g, b, a, nx, ny, nz, 0.8f);
            }
        }
    }

    if (outDepth) *outDepth = depth;
    return builder.commit(pool);
}

struct DagBrickSet::Impl {
    PackedBricks packed;
};

DagBrickSet::DagBrickSet() = default;
DagBrickSet::DagBrickSet(std::unique_ptr<Impl> impl) noexcept : impl_(std::move(impl)) {}
DagBrickSet::~DagBrickSet() = default;
DagBrickSet::DagBrickSet(DagBrickSet&&) noexcept = default;
DagBrickSet& DagBrickSet::operator=(DagBrickSet&&) noexcept = default;

int DagBrickSet::depth() const noexcept { return impl_ ? impl_->packed.depth : 0; }
bool DagBrickSet::empty() const noexcept { return !impl_ || impl_->packed.bricks.empty(); }

uint32_t commitVoxelBricks(VoxelDagPool& pool, const DagBrickSet& set) {
    return set.impl() ? SparseBrickBuilder::commitPacked(pool, set.impl()->packed) : kDagNullNode;
}

uint32_t voxelizeMeshToDag(
    VoxelDagPool& pool,
    const MeshData& mesh,
    const DagVoxelizeOptions& opts,
    int* outDepth) {
    const DagBrickSet set = voxelizeMeshToBricks(mesh, opts);
    if (outDepth) *outDepth = roundDepthForResolution(opts.resolution);
    return commitVoxelBricks(pool, set);
}

DagBrickSet voxelizeMeshToBricks(const MeshData& mesh, const DagVoxelizeOptions& opts) {
    const int depth = roundDepthForResolution(opts.resolution);
    if (mesh.empty() || mesh.triangleCount() == 0) return DagBrickSet{};
    SparseBrickBuilder builder(depth);
    const int res = 4 << depth;
    float invRes = 1.0f, ox = 0.0f, oy = 0.0f, oz = 0.0f;
    if (opts.unitBounds) {
        invRes = 2.0f / static_cast<float>(res);
        ox = -1.0f;
        oy = -1.0f;
        oz = -1.0f;
    } else {
        const bromath::AABB3 box = computeBBox(mesh);
        const float extX = box.max.x - box.min.x;
        const float extY = box.max.y - box.min.y;
        const float extZ = box.max.z - box.min.z;
        const float maxExt = std::max({extX, extY, extZ, 1e-5f});
        invRes = maxExt / static_cast<float>(res - 2);
        ox = 0.5f * (box.min.x + box.max.x) - 0.5f * static_cast<float>(res) * invRes;
        oy = opts.fitToCube ? (box.min.y - invRes) : (0.5f * (box.min.y + box.max.y) - 0.5f * static_cast<float>(res) * invRes);
        oz = 0.5f * (box.min.z + box.max.z) - 0.5f * static_cast<float>(res) * invRes;
    }

    const bool hasNormals = mesh.hasNormals();
    const bool hasColors = mesh.hasColors();
    const float halfSize[3] = {0.5f, 0.5f, 0.5f};

    const size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        const uint32_t i0 = mesh.indices[t * 3 + 0];
        const uint32_t i1 = mesh.indices[t * 3 + 1];
        const uint32_t i2 = mesh.indices[t * 3 + 2];

        const float p0[3] = {
            (mesh.positions[i0 * 3 + 0] - ox) / invRes,
            (mesh.positions[i0 * 3 + 1] - oy) / invRes,
            (mesh.positions[i0 * 3 + 2] - oz) / invRes,
        };
        const float p1[3] = {
            (mesh.positions[i1 * 3 + 0] - ox) / invRes,
            (mesh.positions[i1 * 3 + 1] - oy) / invRes,
            (mesh.positions[i1 * 3 + 2] - oz) / invRes,
        };
        const float p2[3] = {
            (mesh.positions[i2 * 3 + 0] - ox) / invRes,
            (mesh.positions[i2 * 3 + 1] - oy) / invRes,
            (mesh.positions[i2 * 3 + 2] - oz) / invRes,
        };

        const float e1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const float e2[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        float fN[3] = {
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        };
        const float fnLen = std::sqrt(fN[0] * fN[0] + fN[1] * fN[1] + fN[2] * fN[2]);
        if (fnLen > 1e-6f) {
            fN[0] /= fnLen;
            fN[1] /= fnLen;
            fN[2] /= fnLen;
        } else {
            fN[1] = 1.0f;
        }

        const float minX = std::min({p0[0], p1[0], p2[0]}), maxX = std::max({p0[0], p1[0], p2[0]});
        const float minY = std::min({p0[1], p1[1], p2[1]}), maxY = std::max({p0[1], p1[1], p2[1]});
        const float minZ = std::min({p0[2], p1[2], p2[2]}), maxZ = std::max({p0[2], p1[2], p2[2]});
        if (maxX < 0.0f || minX >= static_cast<float>(res) ||
            maxY < 0.0f || minY >= static_cast<float>(res) ||
            maxZ < 0.0f || minZ >= static_cast<float>(res)) {
            continue;
        }
        const int gx0 = std::clamp(static_cast<int>(std::floor(minX)), 0, res - 1);
        const int gx1 = std::clamp(static_cast<int>(std::floor(maxX)), 0, res - 1);
        const int gy0 = std::clamp(static_cast<int>(std::floor(minY)), 0, res - 1);
        const int gy1 = std::clamp(static_cast<int>(std::floor(maxY)), 0, res - 1);
        const int gz0 = std::clamp(static_cast<int>(std::floor(minZ)), 0, res - 1);
        const int gz1 = std::clamp(static_cast<int>(std::floor(maxZ)), 0, res - 1);

        for (int gz = gz0; gz <= gz1; ++gz) {
            for (int gy = gy0; gy <= gy1; ++gy) {
                for (int gx = gx0; gx <= gx1; ++gx) {
                    const float c[3] = {gx + 0.5f, gy + 0.5f, gz + 0.5f};
                    if (!triBoxOverlap(c, halfSize, p0, p1, p2)) continue;

                    float bu, bv, bw;
                    closestBarycentricOnTriangle(c, p0, p1, p2, bu, bv, bw);

                    float nx = fN[0], ny = fN[1], nz = fN[2];
                    if (hasNormals) {
                        nx = bu * mesh.normals[i0 * 3 + 0] + bv * mesh.normals[i1 * 3 + 0] + bw * mesh.normals[i2 * 3 + 0];
                        ny = bu * mesh.normals[i0 * 3 + 1] + bv * mesh.normals[i1 * 3 + 1] + bw * mesh.normals[i2 * 3 + 1];
                        nz = bu * mesh.normals[i0 * 3 + 2] + bv * mesh.normals[i1 * 3 + 2] + bw * mesh.normals[i2 * 3 + 2];
                    }
                    float cr = opts.defaultColor[0];
                    float cg = opts.defaultColor[1];
                    float cb = opts.defaultColor[2];
                    float ca = opts.defaultColor[3];
                    if (hasColors) {
                        cr = bu * mesh.colors[i0 * 4 + 0] + bv * mesh.colors[i1 * 4 + 0] + bw * mesh.colors[i2 * 4 + 0];
                        cg = bu * mesh.colors[i0 * 4 + 1] + bv * mesh.colors[i1 * 4 + 1] + bw * mesh.colors[i2 * 4 + 1];
                        cb = bu * mesh.colors[i0 * 4 + 2] + bv * mesh.colors[i1 * 4 + 2] + bw * mesh.colors[i2 * 4 + 2];
                        ca = bu * mesh.colors[i0 * 4 + 3] + bv * mesh.colors[i1 * 4 + 3] + bw * mesh.colors[i2 * 4 + 3];
                    }
                    builder.addVoxel(gx, gy, gz, opts.state, cr, cg, cb, ca, nx, ny, nz, opts.roughness, opts.emissive);
                }
            }
        }
    }

    if (opts.fillInterior && opts.state == VoxelState::Solid) builder.fillInside(opts.roughness, opts.emissive);
    auto impl = std::make_unique<DagBrickSet::Impl>();
    impl->packed = builder.pack();
    return DagBrickSet(std::move(impl));
}

uint32_t voxelizeSdfToDag(
    VoxelDagPool& pool,
    const std::function<float(float, float, float)>& sdfFn,
    const bromath::AABB3& bounds,
    const DagVoxelizeOptions& opts,
    int* outDepth) {
    const int depth = roundDepthForResolution(opts.resolution);
    const int res = 4 << depth;
    if (!sdfFn) {
        if (outDepth) *outDepth = depth;
        return kDagNullNode;
    }

    const float dx = (bounds.max.x - bounds.min.x) / static_cast<float>(res);
    const float dy = (bounds.max.y - bounds.min.y) / static_cast<float>(res);
    const float dz = (bounds.max.z - bounds.min.z) / static_cast<float>(res);
    const float eps = std::max(std::min({dx, dy, dz}) * 0.5f, 1e-4f);
    const float shellBand = std::max({dx, dy, dz}) * 1.75f;

    SparseBrickBuilder builder(depth);
    for (int gz = 0; gz < res; ++gz) {
        const float wz = bounds.min.z + (static_cast<float>(gz) + 0.5f) * dz;
        for (int gy = 0; gy < res; ++gy) {
            const float wy = bounds.min.y + (static_cast<float>(gy) + 0.5f) * dy;
            for (int gx = 0; gx < res; ++gx) {
                const float wx = bounds.min.x + (static_cast<float>(gx) + 0.5f) * dx;
                const float d = sdfFn(wx, wy, wz);
                if (d > 0.0f) continue;

                float nx = 0.0f, ny = 1.0f, nz = 0.0f;
                if (opts.state == VoxelState::Solid && d > -shellBand) {
                    nx = sdfFn(wx + eps, wy, wz) - sdfFn(wx - eps, wy, wz);
                    ny = sdfFn(wx, wy + eps, wz) - sdfFn(wx, wy - eps, wz);
                    nz = sdfFn(wx, wy, wz + eps) - sdfFn(wx, wy, wz - eps);
                }
                builder.addVoxel(
                    gx, gy, gz, opts.state,
                    opts.defaultColor[0], opts.defaultColor[1], opts.defaultColor[2], opts.defaultColor[3],
                    nx, ny, nz, opts.roughness, opts.emissive);
            }
        }
    }

    if (outDepth) *outDepth = depth;
    return builder.commit(pool);
}

uint32_t loadVoxToDag(
    VoxelDagPool& pool,
    const uint8_t* data,
    size_t size,
    int* outDepth,
    float roughness) {
    if (!data || size < 20 || std::memcmp(data, "VOX ", 4) != 0) {
        if (outDepth) *outDepth = 1;
        return kDagNullNode;
    }

    size_t pos = 8; // skip "VOX " + version
    int sx = 0, sy = 0, sz = 0;
    struct VoxCell { uint8_t x, y, z, c; };
    std::vector<VoxCell> cells;
    float palette[256 * 4] = {};
    bool hasPalette = false;

    auto readChunk = [&](auto& self, size_t limit) -> void {
        while (pos + 12 <= limit) {
            const char* id = reinterpret_cast<const char*>(data + pos);
            const uint32_t contentSize = readU32LE(data + pos + 4);
            const uint32_t childSize = readU32LE(data + pos + 8);
            pos += 12;
            const size_t contentStart = pos;
            const size_t chunkEnd = std::min(limit, pos + static_cast<size_t>(contentSize) + static_cast<size_t>(childSize));

            if (std::memcmp(id, "MAIN", 4) == 0) {
                pos = std::min(chunkEnd, contentStart + contentSize);
                self(self, chunkEnd);
            } else if (std::memcmp(id, "SIZE", 4) == 0 && contentSize >= 12 && contentStart + 12 <= limit) {
                sx = readI32LE(data + contentStart + 0);
                sy = readI32LE(data + contentStart + 4);
                sz = readI32LE(data + contentStart + 8);
                pos = chunkEnd;
            } else if (std::memcmp(id, "XYZI", 4) == 0 && contentSize >= 4 && contentStart + 4 <= limit) {
                const uint32_t count = readU32LE(data + contentStart);
                const size_t maxCount = (limit - (contentStart + 4)) / 4;
                const size_t n = std::min(static_cast<size_t>(count), maxCount);
                cells.reserve(n);
                for (size_t i = 0; i < n; ++i) {
                    const uint8_t* p = data + contentStart + 4 + i * 4;
                    cells.push_back({p[0], p[1], p[2], p[3]});
                }
                pos = chunkEnd;
            } else if (std::memcmp(id, "RGBA", 4) == 0) {
                const size_t avail = std::min<size_t>({contentSize, 255 * 4, limit - contentStart});
                for (size_t i = 0; i * 4 + 3 < avail; ++i) {
                    const size_t pi = (i + 1) * 4;
                    palette[pi + 0] = static_cast<float>(data[contentStart + i * 4 + 0]) / 255.0f;
                    palette[pi + 1] = static_cast<float>(data[contentStart + i * 4 + 1]) / 255.0f;
                    palette[pi + 2] = static_cast<float>(data[contentStart + i * 4 + 2]) / 255.0f;
                    palette[pi + 3] = static_cast<float>(data[contentStart + i * 4 + 3]) / 255.0f;
                }
                hasPalette = true;
                pos = chunkEnd;
            } else {
                pos = chunkEnd;
            }
        }
    };

    readChunk(readChunk, size);
    if (sx <= 0 || sy <= 0 || sz <= 0 || cells.empty()) {
        if (outDepth) *outDepth = 1;
        return kDagNullNode;
    }

    if (!hasPalette) {
        for (int i = 1; i < 256; ++i) {
            palette[i * 4 + 0] = 0.25f + 0.7f * static_cast<float>((i * 67) % 255) / 255.0f;
            palette[i * 4 + 1] = 0.25f + 0.7f * static_cast<float>((i * 131) % 255) / 255.0f;
            palette[i * 4 + 2] = 0.25f + 0.7f * static_cast<float>((i * 197) % 255) / 255.0f;
            palette[i * 4 + 3] = 1.0f;
        }
    }

    // MagicaVoxel .vox is Z-up; convert to Y-up centered in the octree grid.
    const int maxDim = std::max({sx, sy, sz});
    const int depth = roundDepthForResolution(maxDim);
    const int res = 4 << depth;
    const int offX = (res - sx) / 2;
    const int offY = 0;
    const int offZ = (res - sy) / 2;

    std::vector<uint8_t> occ(static_cast<size_t>(sx) * sy * sz, 0);
    for (const VoxCell& v : cells) {
        if (v.x < sx && v.y < sy && v.z < sz && v.c != 0) {
            occ[(static_cast<size_t>(v.z) * sy + v.y) * sx + v.x] = v.c;
        }
    }
    auto solidAt = [&](int x, int y, int z) -> bool {
        if (x < 0 || y < 0 || z < 0 || x >= sx || y >= sy || z >= sz) return false;
        return occ[(static_cast<size_t>(z) * sy + y) * sx + x] != 0;
    };

    SparseBrickBuilder builder(depth);
    for (const VoxCell& v : cells) {
        if (v.c == 0 || v.x >= sx || v.y >= sy || v.z >= sz) continue;
        const int vx = v.x, vy = v.y, vz = v.z;
        const float nx = (solidAt(vx - 1, vy, vz) ? 1.0f : 0.0f) - (solidAt(vx + 1, vy, vz) ? 1.0f : 0.0f);
        const float ny = (solidAt(vx, vy, vz - 1) ? 1.0f : 0.0f) - (solidAt(vx, vy, vz + 1) ? 1.0f : 0.0f);
        const float nz = (solidAt(vx, vy + 1, vz) ? 1.0f : 0.0f) - (solidAt(vx, vy - 1, vz) ? 1.0f : 0.0f);
        const int gx = offX + vx;
        const int gy = offY + vz;
        const int gz = offZ + (sy - 1 - vy);
        const size_t pi = static_cast<size_t>(v.c) * 4u;
        builder.addVoxel(gx, gy, gz, VoxelState::Solid,
                         palette[pi + 0], palette[pi + 1], palette[pi + 2], palette[pi + 3],
                         nx, ny, nz, roughness);
    }

    if (outDepth) *outDepth = depth;
    return builder.commit(pool);
}

} // namespace bromesh
