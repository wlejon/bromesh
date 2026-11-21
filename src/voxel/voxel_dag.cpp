#include "bromesh/voxel/voxel_dag.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

namespace {

inline uint8_t quantizeU8(float v) noexcept {
    const float c = std::clamp(v, 0.0f, 1.0f);
    return static_cast<uint8_t>(std::lround(c * 255.0f));
}

inline uint32_t fnvMix(uint32_t h, uint32_t v) noexcept {
    h ^= v;
    h *= 16777619u;
    return h;
}

inline float octSignNotZero(float v) noexcept {
    return (v >= 0.0f) ? 1.0f : -1.0f;
}

} // namespace

namespace dag_pack {

uint32_t packRgba8(float r, float g, float b, float a) {
    const uint32_t ur = quantizeU8(r);
    const uint32_t ug = quantizeU8(g);
    const uint32_t ub = quantizeU8(b);
    const uint32_t ua = quantizeU8(a);
    return ur | (ug << 8u) | (ub << 16u) | (ua << 24u);
}

void unpackRgba8(uint32_t packed, float& r, float& g, float& b, float& a) {
    r = static_cast<float>(packed & 0xFFu) / 255.0f;
    g = static_cast<float>((packed >> 8u) & 0xFFu) / 255.0f;
    b = static_cast<float>((packed >> 16u) & 0xFFu) / 255.0f;
    a = static_cast<float>((packed >> 24u) & 0xFFu) / 255.0f;
}

uint16_t packOctNormal(float nx, float ny, float nz) {
    const float len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len < 1e-6f) {
        nx = 0.0f;
        ny = 1.0f;
        nz = 0.0f;
    } else {
        nx /= len;
        ny /= len;
        nz /= len;
    }
    const float invL1 = 1.0f / (std::abs(nx) + std::abs(ny) + std::abs(nz));
    float px = nx * invL1;
    float py = ny * invL1;
    if (nz < 0.0f) {
        const float ox = (1.0f - std::abs(py)) * octSignNotZero(px);
        const float oy = (1.0f - std::abs(px)) * octSignNotZero(py);
        px = ox;
        py = oy;
    }
    const uint32_t ux = quantizeU8(px * 0.5f + 0.5f);
    const uint32_t uy = quantizeU8(py * 0.5f + 0.5f);
    return static_cast<uint16_t>(ux | (uy << 8u));
}

void unpackOctNormal(uint16_t packed, float& nx, float& ny, float& nz) {
    const float fx = static_cast<float>(packed & 0xFFu) / 255.0f * 2.0f - 1.0f;
    const float fy = static_cast<float>((packed >> 8u) & 0xFFu) / 255.0f * 2.0f - 1.0f;
    float fz = 1.0f - std::abs(fx) - std::abs(fy);
    float rx = fx;
    float ry = fy;
    if (fz < 0.0f) {
        rx = (1.0f - std::abs(fy)) * octSignNotZero(fx);
        ry = (1.0f - std::abs(fx)) * octSignNotZero(fy);
    }
    const float len = std::sqrt(rx * rx + ry * ry + fz * fz);
    if (len > 1e-6f) {
        nx = rx / len;
        ny = ry / len;
        nz = fz / len;
    } else {
        nx = 0.0f;
        ny = 1.0f;
        nz = 0.0f;
    }
}

uint32_t packNormalMat(float nx, float ny, float nz, float roughness, float emissive) {
    const uint32_t oct = packOctNormal(nx, ny, nz);
    const uint32_t r8 = quantizeU8(roughness);
    const uint32_t e8 = quantizeU8(emissive);
    return oct | (r8 << 16u) | (e8 << 24u);
}

void unpackNormalMat(uint32_t packed, float& nx, float& ny, float& nz, float& roughness, float& emissive) {
    unpackOctNormal(static_cast<uint16_t>(packed & 0xFFFFu), nx, ny, nz);
    roughness = static_cast<float>((packed >> 16u) & 0xFFu) / 255.0f;
    emissive = static_cast<float>((packed >> 24u) & 0xFFu) / 255.0f;
}

} // namespace dag_pack

// ─── Keys ───────────────────────────────────────────────────────────────────

bool VoxelDagPool::BrickKey::operator==(const BrickKey& o) const noexcept {
    if (solidLo != o.solidLo || solidHi != o.solidHi || carveLo != o.carveLo || carveHi != o.carveHi ||
        packedColorOpacity != o.packedColorOpacity || packedNormalMat != o.packedNormalMat ||
        flags != o.flags || inlineAttrs.size() != o.inlineAttrs.size()) {
        return false;
    }
    for (size_t i = 0; i < inlineAttrs.size(); ++i) {
        if (inlineAttrs[i].colorOpacity != o.inlineAttrs[i].colorOpacity ||
            inlineAttrs[i].normalMat != o.inlineAttrs[i].normalMat) {
            return false;
        }
    }
    return true;
}

size_t VoxelDagPool::BrickKeyHash::operator()(const BrickKey& k) const noexcept {
    uint32_t h = 2166136261u;
    h = fnvMix(h, k.solidLo);
    h = fnvMix(h, k.solidHi);
    h = fnvMix(h, k.carveLo);
    h = fnvMix(h, k.carveHi);
    h = fnvMix(h, k.packedColorOpacity);
    h = fnvMix(h, k.packedNormalMat);
    h = fnvMix(h, k.flags);
    h = fnvMix(h, k.attrHash);
    uint64_t x = h;    // FNV's low bits see only the inputs' low bits: spread them
    x ^= x >> 16u;
    x *= 0x94D049BB133111EBull;
    x ^= x >> 31u;
    return static_cast<size_t>(x);
}

bool VoxelDagPool::NodeKey::operator==(const NodeKey& o) const noexcept {
    if (masks != o.masks || count != o.count) return false;
    for (uint32_t i = 0; i < count; ++i) {
        if (childRefs[i] != o.childRefs[i]) return false;
    }
    return true;
}

size_t VoxelDagPool::NodeKeyHash::operator()(const NodeKey& k) const noexcept {
    // Every input bit reaches every output bit: the table buckets by the low
    // bits, and child refs differ anywhere.
    uint64_t h = 0x9E3779B97F4A7C15ull ^ (static_cast<uint64_t>(k.masks) << 32u | k.count);
    for (uint32_t i = 0; i < k.count; ++i) {
        h = (h ^ k.childRefs[i]) * 0xBF58476D1CE4E5B9ull;
        h ^= h >> 31u;
    }
    h ^= h >> 30u;
    h *= 0x94D049BB133111EBull;
    h ^= h >> 31u;
    return static_cast<size_t>(h);
}

// ─── VoxelDagPool ───────────────────────────────────────────────────────────

void VoxelDagPool::clear() {
    nodes_.clear();
    children_.clear();
    bricks_.clear();
    attrs_.clear();
    animClips_.clear();
    animFrames_.clear();
    brickMap_.clear();
    nodeMap_.clear();
    uniformCarveNodes_.clear();
    freeNodes_.clear();
    freeBricks_.clear();
    freeChildren_.clear();
    freeAttrs_.clear();
    nodeRefs_.clear();
    brickRefs_.clear();
    heldNodes_.clear();
    heldBricks_.clear();
    freshNodes_.clear();
    freshBricks_.clear();
    for (auto& p : pageGen_) p.clear();
    ++generation_;
}

void VoxelDagPool::touch(Array a, uint32_t first, uint32_t count) {
    if (count == 0u) return;
    auto& gens = pageGen_[static_cast<uint32_t>(a)];
    const uint32_t p0 = first / kPageElems;
    const uint32_t p1 = (first + count - 1u) / kPageElems;
    if (gens.size() <= p1) gens.resize(p1 + 1u, 0u);
    for (uint32_t p = p0; p <= p1; ++p) gens[p] = generation_;
}

uint32_t VoxelDagPool::allocNode() {
    if (!freeNodes_.empty()) {
        const uint32_t i = freeNodes_.back();
        freeNodes_.pop_back();
        nodeRefs_[i] = 0u;
        return i;
    }
    nodes_.emplace_back();
    nodeRefs_.push_back(0u);
    return static_cast<uint32_t>(nodes_.size() - 1);
}

uint32_t VoxelDagPool::allocBrick() {
    if (!freeBricks_.empty()) {
        const uint32_t i = freeBricks_.back();
        freeBricks_.pop_back();
        brickRefs_[i] = 0u;
        return i;
    }
    bricks_.emplace_back();
    brickRefs_.push_back(0u);
    return static_cast<uint32_t>(bricks_.size() - 1);
}

uint32_t VoxelDagPool::allocChildren(uint32_t count) {
    if (count < freeChildren_.size() && !freeChildren_[count].empty()) {
        const uint32_t base = freeChildren_[count].back();
        freeChildren_[count].pop_back();
        return base;
    }
    const uint32_t base = static_cast<uint32_t>(children_.size());
    children_.resize(children_.size() + count, kDagNullNode);
    return base;
}

uint32_t VoxelDagPool::allocAttrs(uint32_t count) {
    if (count < freeAttrs_.size() && !freeAttrs_[count].empty()) {
        const uint32_t base = freeAttrs_[count].back();
        freeAttrs_[count].pop_back();
        return base;
    }
    const uint32_t base = static_cast<uint32_t>(attrs_.size());
    attrs_.resize(attrs_.size() + count);
    return base;
}

void VoxelDagPool::retain(uint32_t ref, bool isLeaf) {
    if (ref == kDagNullNode || (!isLeaf && (ref & kDagAnimBit) != 0u)) return;
    if (isLeaf) {
        if (ref >= bricks_.size()) return;
        ++brickRefs_[ref];
        ++heldBricks_[ref];
    } else {
        if (ref >= nodes_.size()) return;
        ++nodeRefs_[ref];
        ++heldNodes_[ref];
    }
}

void VoxelDagPool::release(uint32_t ref, bool isLeaf) {
    if (ref == kDagNullNode || (!isLeaf && (ref & kDagAnimBit) != 0u)) return;
    auto& held = isLeaf ? heldBricks_ : heldNodes_;
    auto it = held.find(ref);
    if (it == held.end()) return;    // never held: not ours to drop
    if (--it->second == 0u) held.erase(it);
    uint32_t& refs = isLeaf ? brickRefs_[ref] : nodeRefs_[ref];
    if (refs > 0u && --refs == 0u) freeUnreferenced(ref, isLeaf);
}

size_t VoxelDagPool::sweep() {
    size_t freed = 0;
    // A fresh slot may have been freed (and even reused) since it was made:
    // whatever is in it now, it is garbage if it is live and unreferenced.
    for (uint32_t b : freshBricks_) {
        if (b < bricks_.size() && brickRefs_[b] == 0u) {
            const DagBrick& br = bricks_[b];
            if ((br.solidLo | br.solidHi | br.carveLo | br.carveHi) != 0u) freed += freeUnreferenced(b, true);
        }
    }
    for (uint32_t n : freshNodes_) {
        if (n < nodes_.size() && nodeRefs_[n] == 0u && dag_pack::validMaskOf(nodes_[n]) != 0u) {
            freed += freeUnreferenced(n, false);
        }
    }
    freshBricks_.clear();
    freshNodes_.clear();
    return freed;
}

// `ref` has no references left: free it, and what that leaves unreferenced.
size_t VoxelDagPool::freeUnreferenced(uint32_t ref, bool isLeaf) {
    size_t freed = 0;
    std::vector<std::pair<uint32_t, bool>> stack{{ref, isLeaf}};
    while (!stack.empty()) {
        const auto [r, leaf] = stack.back();
        stack.pop_back();
        if (leaf) {
            freeBrickSlot(r);
            ++freed;
            continue;
        }
        const DagNode node = nodes_[r];
        const uint32_t vm = dag_pack::validMaskOf(node);
        const uint32_t lm = dag_pack::leafMaskOf(node);
        uint32_t slot = 0u;
        for (uint32_t o = 0u; o < 8u; ++o) {
            if ((vm & (1u << o)) == 0u) continue;
            const uint32_t c = children_[node.childBase + slot++];
            if (c == kDagNullNode || (c & kDagAnimBit) != 0u) continue;
            const bool cl = (lm & (1u << o)) != 0u;
            uint32_t& refs = cl ? brickRefs_[c] : nodeRefs_[c];
            if (refs > 0u && --refs == 0u) stack.emplace_back(c, cl);
        }
        freeNodeSlot(r);
        ++freed;
    }
    return freed;
}

void VoxelDagPool::freeNodeSlot(uint32_t i) {
    DagNode& node = nodes_[i];
    const uint32_t count = dag_pack::popcount8(dag_pack::validMaskOf(node));
    if (count == 0u) return;    // already free
    NodeKey key;
    key.masks = node.masks;
    key.count = count;
    for (uint32_t k = 0; k < count; ++k) key.childRefs[k] = children_[node.childBase + k];
    if (auto it = nodeMap_.find(key); it != nodeMap_.end() && it->second == i) nodeMap_.erase(it);
    if (freeChildren_.size() < 9) freeChildren_.resize(9);
    freeChildren_[count].push_back(node.childBase);
    node = DagNode{};
    freeNodes_.push_back(i);
}

void VoxelDagPool::freeBrickSlot(uint32_t i) {
    DagBrick& b = bricks_[i];
    if ((b.solidLo | b.solidHi | b.carveLo | b.carveHi) == 0u) return;    // already free
    const uint32_t solid = dag_pack::popcount32(b.solidLo) + dag_pack::popcount32(b.solidHi);
    const uint32_t attrCount = solid == 0u ? 0u : ((b.flags & kBrickFlagUniformAttr) ? 1u : solid);
    BrickKey key;
    key.solidLo = b.solidLo;
    key.solidHi = b.solidHi;
    key.carveLo = b.carveLo;
    key.carveHi = b.carveHi;
    key.packedColorOpacity = b.packedColorOpacity;
    key.packedNormalMat = b.packedNormalMat;
    key.flags = b.flags;
    key.inlineAttrs.assign(attrs_.begin() + b.attrOffset, attrs_.begin() + b.attrOffset + attrCount);
    uint32_t attrH = 0u;
    for (const DagVoxelAttr& a : key.inlineAttrs) {
        attrH = fnvMix(attrH, a.colorOpacity);
        attrH = fnvMix(attrH, a.normalMat);
    }
    key.attrHash = attrH;
    if (auto it = brickMap_.find(key); it != brickMap_.end() && it->second == i) brickMap_.erase(it);
    if (freeAttrs_.size() < 65) freeAttrs_.resize(65);
    if (attrCount > 0u) freeAttrs_[attrCount].push_back(b.attrOffset);
    b = DagBrick{};
    freeBricks_.push_back(i);
}

// Every live node's and brick's references, from the live nodes and the holds.
void VoxelDagPool::recountRefs() {
    std::fill(nodeRefs_.begin(), nodeRefs_.end(), 0u);
    std::fill(brickRefs_.begin(), brickRefs_.end(), 0u);
    for (const DagNode& node : nodes_) {
        const uint32_t vm = dag_pack::validMaskOf(node);
        if (vm == 0u) continue;
        const uint32_t lm = dag_pack::leafMaskOf(node);
        uint32_t slot = 0u;
        for (uint32_t o = 0u; o < 8u; ++o) {
            if ((vm & (1u << o)) == 0u) continue;
            const uint32_t c = children_[node.childBase + slot++];
            if (c == kDagNullNode || (c & kDagAnimBit) != 0u) continue;
            if (lm & (1u << o)) ++brickRefs_[c]; else ++nodeRefs_[c];
        }
    }
    for (const auto& [r, n] : heldNodes_) nodeRefs_[r] += n;
    for (const auto& [r, n] : heldBricks_) brickRefs_[r] += n;
}

size_t VoxelDagPool::collect(const std::vector<uint32_t>& rootNodes) {
    std::vector<uint8_t> liveNode(nodes_.size(), 0u);
    std::vector<uint8_t> liveBrick(bricks_.size(), 0u);
    std::vector<uint32_t> stack;
    const auto markBrick = [&](uint32_t b) {
        if (b < liveBrick.size()) liveBrick[b] = 1u;
    };
    const auto markNode = [&](uint32_t n) {
        if (n < liveNode.size() && !liveNode[n]) {
            liveNode[n] = 1u;
            stack.push_back(n);
        }
    };
    for (uint32_t r : rootNodes) markNode(r);
    for (const auto& [r, _] : heldNodes_) markNode(r);
    for (const auto& [r, _] : heldBricks_) markBrick(r);
    for (const auto& [_, n] : uniformCarveNodes_) markNode(n);
    // A clip's frames are kept whole, read either way (a frame of leaf
    // bricks or of nodes): keeping too much is safe, freeing a frame is not.
    for (uint32_t f : animFrames_) {
        markNode(f);
        markBrick(f);
    }
    while (!stack.empty()) {
        const uint32_t n = stack.back();
        stack.pop_back();
        const DagNode& node = nodes_[n];
        const uint32_t vm = dag_pack::validMaskOf(node);
        const uint32_t lm = dag_pack::leafMaskOf(node);
        uint32_t slot = 0u;
        for (uint32_t o = 0u; o < 8u; ++o) {
            if ((vm & (1u << o)) == 0u) continue;
            const uint32_t c = children_[node.childBase + slot++];
            if (c == kDagNullNode || (c & kDagAnimBit) != 0u) continue;
            if (lm & (1u << o)) markBrick(c); else markNode(c);
        }
    }

    // The maps first, by what they point at: what the key held is not rebuilt.
    std::erase_if(nodeMap_, [&](const auto& kv) { return !liveNode[kv.second]; });
    std::erase_if(brickMap_, [&](const auto& kv) { return !liveBrick[kv.second]; });

    std::vector<uint8_t> wasFree(nodes_.size(), 0u);
    for (uint32_t i : freeNodes_) wasFree[i] = 1u;
    if (freeChildren_.size() < 9) freeChildren_.resize(9);
    size_t freed = 0;
    for (uint32_t i = 0; i < nodes_.size(); ++i) {
        if (liveNode[i] || wasFree[i]) continue;
        DagNode& node = nodes_[i];
        const uint32_t count = dag_pack::popcount8(dag_pack::validMaskOf(node));
        if (count > 0u) freeChildren_[count].push_back(node.childBase);
        node = DagNode{};
        freeNodes_.push_back(i);
        ++freed;
    }
    std::vector<uint8_t> brickWasFree(bricks_.size(), 0u);
    for (uint32_t i : freeBricks_) brickWasFree[i] = 1u;
    if (freeAttrs_.size() < 65) freeAttrs_.resize(65);
    for (uint32_t i = 0; i < bricks_.size(); ++i) {
        if (liveBrick[i] || brickWasFree[i]) continue;
        DagBrick& b = bricks_[i];
        const uint32_t solid = dag_pack::popcount32(b.solidLo) + dag_pack::popcount32(b.solidHi);
        const uint32_t attrCount = solid == 0u ? 0u : ((b.flags & kBrickFlagUniformAttr) ? 1u : solid);
        if (attrCount > 0u) freeAttrs_[attrCount].push_back(b.attrOffset);
        b = DagBrick{};
        freeBricks_.push_back(i);
        ++freed;
    }
    freshNodes_.clear();
    freshBricks_.clear();
    recountRefs();
    if (freed > 0) ++generation_;
    return freed;
}

uint32_t VoxelDagPool::internBrick(
    uint32_t solidLo,
    uint32_t solidHi,
    uint32_t carveLo,
    uint32_t carveHi,
    const std::vector<DagVoxelAttr>& solidAttrs,
    const DagVoxelAttr* uniformAttr) {
    if ((solidLo | solidHi | carveLo | carveHi) == 0u) {
        return kDagNullNode;
    }

    const uint32_t solidCount = dag_pack::popcount32(solidLo) + dag_pack::popcount32(solidHi);
    uint32_t flags = 0u;
    if (carveLo == 0xFFFFFFFFu && carveHi == 0xFFFFFFFFu && solidCount == 0u) {
        flags |= kBrickFlagUniformCarve;
    }

    std::vector<DagVoxelAttr> storeAttrs;
    if (solidCount > 0u) {
        if (uniformAttr != nullptr) {
            flags |= kBrickFlagUniformAttr;
            storeAttrs.push_back(*uniformAttr);
        } else if (solidAttrs.size() == solidCount) {
            bool allSame = true;
            for (size_t i = 1; i < solidAttrs.size(); ++i) {
                if (solidAttrs[i].colorOpacity != solidAttrs[0].colorOpacity ||
                    solidAttrs[i].normalMat != solidAttrs[0].normalMat) {
                    allSame = false;
                    break;
                }
            }
            if (allSame) {
                flags |= kBrickFlagUniformAttr;
                storeAttrs.push_back(solidAttrs[0]);
            } else {
                storeAttrs = solidAttrs;
            }
        } else {
            flags |= kBrickFlagUniformAttr;
            storeAttrs.push_back(
                DagVoxelAttr{dag_pack::packRgba8(0.7f, 0.7f, 0.7f, 1.0f), dag_pack::packNormalMat(0.0f, 1.0f, 0.0f, 0.8f)});
        }
    }

    // Compute pre-filtered brick mip attribute & directional projection coverage.
    uint32_t mipColorOpacity = 0u;
    uint32_t mipNormalMat = 0u;
    if (solidCount > 0u) {
        double sumR = 0.0, sumG = 0.0, sumB = 0.0, sumW = 0.0;
        double sumNx = 0.0, sumNy = 0.0, sumNz = 0.0, sumRough = 0.0, maxEmis = 0.0;

        const size_t n = (flags & kBrickFlagUniformAttr) ? 1u : storeAttrs.size();
        for (size_t i = 0; i < n; ++i) {
            float r, g, b, a;
            dag_pack::unpackRgba8(storeAttrs[i].colorOpacity, r, g, b, a);
            const double w = std::max(0.01, static_cast<double>(a));
            sumR += r * w;
            sumG += g * w;
            sumB += b * w;
            sumW += w;

            float nx, ny, nz, rough, emis;
            dag_pack::unpackNormalMat(storeAttrs[i].normalMat, nx, ny, nz, rough, emis);
            sumNx += nx * w;
            sumNy += ny * w;
            sumNz += nz * w;
            sumRough += rough * w;
            maxEmis = std::max(maxEmis, static_cast<double>(emis));
        }

        const uint64_t solidBits = (static_cast<uint64_t>(solidHi) << 32u) | solidLo;
        uint32_t projX = 0, projY = 0, projZ = 0;
        for (int z = 0; z < 4; ++z) {
            for (int y = 0; y < 4; ++y) {
                for (int x = 0; x < 4; ++x) {
                    const int b = x + (y << 2) + (z << 4);
                    if ((solidBits >> b) & 1ull) {
                        projX |= (1u << (y + (z << 2)));
                        projY |= (1u << (x + (z << 2)));
                        projZ |= (1u << (x + (y << 2)));
                    }
                }
            }
        }
        const float px = static_cast<float>(dag_pack::popcount32(projX)) / 16.0f;
        const float py = static_cast<float>(dag_pack::popcount32(projY)) / 16.0f;
        const float pz = static_cast<float>(dag_pack::popcount32(projZ)) / 16.0f;
        const float coverage = std::clamp(0.5f * std::max({px, py, pz}) + 0.5f * ((px + py + pz) / 3.0f), 0.0f, 1.0f);
        const float meanA = static_cast<float>(sumW / static_cast<double>(n));
        const float outA = std::clamp(meanA * coverage, 1.0f / 255.0f, 1.0f);

        mipColorOpacity = dag_pack::packRgba8(
            static_cast<float>(sumR / sumW),
            static_cast<float>(sumG / sumW),
            static_cast<float>(sumB / sumW),
            outA);
        mipNormalMat = dag_pack::packNormalMat(
            static_cast<float>(sumNx),
            static_cast<float>(sumNy),
            static_cast<float>(sumNz),
            static_cast<float>(sumRough / sumW),
            static_cast<float>(maxEmis));
    }

    uint32_t attrH = 0u;
    for (const DagVoxelAttr& a : storeAttrs) {
        attrH = fnvMix(attrH, a.colorOpacity);
        attrH = fnvMix(attrH, a.normalMat);
    }

    BrickKey key;
    key.solidLo = solidLo;
    key.solidHi = solidHi;
    key.carveLo = carveLo;
    key.carveHi = carveHi;
    key.packedColorOpacity = mipColorOpacity;
    key.packedNormalMat = mipNormalMat;
    key.flags = flags;
    key.attrHash = attrH;
    key.inlineAttrs = storeAttrs;

    auto it = brickMap_.find(key);
    if (it != brickMap_.end()) {
        return it->second;
    }

    DagBrick brick;
    brick.solidLo = solidLo;
    brick.solidHi = solidHi;
    brick.carveLo = carveLo;
    brick.carveHi = carveHi;
    brick.attrOffset = 0u;
    brick.packedColorOpacity = mipColorOpacity;
    brick.packedNormalMat = mipNormalMat;
    brick.flags = flags;

    if (!storeAttrs.empty()) {
        brick.attrOffset = allocAttrs(static_cast<uint32_t>(storeAttrs.size()));
        std::copy(storeAttrs.begin(), storeAttrs.end(), attrs_.begin() + brick.attrOffset);
    }

    const uint32_t idx = allocBrick();
    bricks_[idx] = brick;
    brickMap_.emplace(std::move(key), idx);
    ++generation_;
    touch(Array::Bricks, idx, 1u);
    touch(Array::Attrs, brick.attrOffset, static_cast<uint32_t>(storeAttrs.size()));
    freshBricks_.push_back(idx);
    return idx;
}

void VoxelDagPool::describeRef(
    uint32_t rawRef,
    bool isLeaf,
    uint8_t& outHasSolid,
    uint8_t& outHasCarve,
    uint8_t& outUniformCarve,
    uint32_t& outColorOpacity,
    uint32_t& outNormalMat) const {
    outHasSolid = 0;
    outHasCarve = 0;
    outUniformCarve = 0;
    outColorOpacity = 0;
    outNormalMat = 0;
    if (rawRef == kDagNullNode) return;

    if (!isLeaf && (rawRef & kDagAnimBit) != 0u) {
        const uint32_t clipIdx = (rawRef >> 16u) & 0x7FFFu;
        if (clipIdx < animClips_.size() && animClips_[clipIdx].frameCount > 0u) {
            rawRef = animFrames_[animClips_[clipIdx].frameOffset];
        } else {
            return;
        }
    }

    if (isLeaf) {
        if (rawRef >= bricks_.size()) return;
        const DagBrick& b = bricks_[rawRef];
        outHasSolid = ((b.solidLo | b.solidHi) != 0u) ? 1u : 0u;
        outHasCarve = ((b.carveLo | b.carveHi) != 0u) ? 1u : 0u;
        outUniformCarve = ((b.flags & kBrickFlagUniformCarve) != 0u) ? 1u : 0u;
        outColorOpacity = b.packedColorOpacity;
        outNormalMat = b.packedNormalMat;
        return;
    }

    if (rawRef >= nodes_.size()) return;
    const DagNode& n = nodes_[rawRef];
    outHasSolid = (dag_pack::solidMaskOf(n) != 0u) ? 1u : 0u;
    outHasCarve = (dag_pack::carveMaskOf(n) != 0u) ? 1u : 0u;
    outUniformCarve = ((dag_pack::flagsOf(n) & kDagFlagUniformCarve) != 0u) ? 1u : 0u;
    outColorOpacity = n.packedColorOpacity;
    outNormalMat = n.packedNormalMat;
}

uint32_t VoxelDagPool::internNode(
    const uint32_t childRefs[8],
    const uint8_t childIsLeaf[8],
    const uint8_t childHasSolid[8],
    const uint8_t childHasCarve[8],
    const uint8_t childUniformCarve[8]) {
    uint32_t solidMask = 0u;
    uint32_t carveMask = 0u;
    uint32_t leafMask = 0u;
    uint32_t count = 0u;
    uint32_t packedChildren[8] = {};

    double sumR = 0.0, sumG = 0.0, sumB = 0.0, sumW = 0.0;
    double sumNx = 0.0, sumNy = 0.0, sumNz = 0.0, sumRough = 0.0, maxEmis = 0.0;
    uint32_t uniformCarveCount = 0u;

    for (uint32_t o = 0u; o < 8u; ++o) {
        if (childRefs[o] == kDagNullNode) continue;
        const bool isLeaf = childIsLeaf[o] != 0;
        const bool hasSolid = childHasSolid[o] != 0;
        const bool hasCarve = childHasCarve[o] != 0;
        if (!hasSolid && !hasCarve) continue;

        const uint32_t bit = 1u << o;
        if (hasSolid) solidMask |= bit;
        if (hasCarve) carveMask |= bit;
        if (isLeaf) leafMask |= bit;
        if (childUniformCarve[o] != 0 && !hasSolid) ++uniformCarveCount;
        packedChildren[count++] = childRefs[o];

        if (hasSolid) {
            uint8_t s, c, u;
            uint32_t colOp = 0, nrmMat = 0;
            describeRef(childRefs[o], isLeaf, s, c, u, colOp, nrmMat);
            float r, g, b, a;
            dag_pack::unpackRgba8(colOp, r, g, b, a);
            const double w = std::max(0.01, static_cast<double>(a));
            sumR += r * w;
            sumG += g * w;
            sumB += b * w;
            sumW += w;

            float nx, ny, nz, rough, emis;
            dag_pack::unpackNormalMat(nrmMat, nx, ny, nz, rough, emis);
            sumNx += nx * w;
            sumNy += ny * w;
            sumNz += nz * w;
            sumRough += rough * w;
            maxEmis = std::max(maxEmis, static_cast<double>(emis));
        }
    }

    if (count == 0u) return kDagNullNode;

    uint32_t flags = 0u;
    if (solidMask == 0u && uniformCarveCount == 8u) {
        flags |= kDagFlagUniformCarve;
    }
    if (maxEmis > 0.01) {
        flags |= kDagFlagEmissive;
    }

    const uint32_t masks = dag_pack::makeNodeMasks(solidMask, carveMask, leafMask, flags);

    NodeKey key;
    key.masks = masks;
    key.count = count;
    for (uint32_t i = 0; i < count; ++i) key.childRefs[i] = packedChildren[i];

    auto it = nodeMap_.find(key);
    if (it != nodeMap_.end()) {
        return it->second;
    }

    DagNode node;
    node.masks = masks;
    node.childBase = allocChildren(count);
    if (sumW > 0.0) {
        // Divide accumulated child opacity by 4 rather than 8 so planar surfaces remain opaque across mip levels.
        const float outA = std::clamp(static_cast<float>(sumW / 4.0), 1.0f / 255.0f, 1.0f);
        node.packedColorOpacity = dag_pack::packRgba8(
            static_cast<float>(sumR / sumW),
            static_cast<float>(sumG / sumW),
            static_cast<float>(sumB / sumW),
            outA);
        node.packedNormalMat = dag_pack::packNormalMat(
            static_cast<float>(sumNx),
            static_cast<float>(sumNy),
            static_cast<float>(sumNz),
            static_cast<float>(sumRough / sumW),
            static_cast<float>(maxEmis));
    }

    std::copy(packedChildren, packedChildren + count, children_.begin() + node.childBase);
    const uint32_t idx = allocNode();
    nodes_[idx] = node;
    nodeMap_.emplace(key, idx);
    // The node holds each child (an animated child's frames are held by its clip).
    uint32_t slot = 0u;
    for (uint32_t o = 0u; o < 8u; ++o) {
        if ((dag_pack::validMaskOf(node) & (1u << o)) == 0u) continue;
        const uint32_t c = packedChildren[slot++];
        if (c == kDagNullNode || (c & kDagAnimBit) != 0u) continue;
        if (leafMask & (1u << o)) ++brickRefs_[c]; else ++nodeRefs_[c];
    }
    ++generation_;
    touch(Array::Nodes, idx, 1u);
    touch(Array::Children, node.childBase, count);
    freshNodes_.push_back(idx);
    return idx;
}

uint32_t VoxelDagPool::internUniformCarveNode(int depthRemaining) {
    if (depthRemaining <= 0) return kDagNullNode;
    const uint64_t key = static_cast<uint64_t>(depthRemaining);
    if (auto it = uniformCarveNodes_.find(key); it != uniformCarveNodes_.end()) {
        return it->second;
    }

    uint32_t child = kDagNullNode;
    bool isLeaf = false;
    if (depthRemaining == 1) {
        child = internBrick(0u, 0u, 0xFFFFFFFFu, 0xFFFFFFFFu, {});
        isLeaf = true;
    } else {
        child = internUniformCarveNode(depthRemaining - 1);
    }

    uint32_t refs[8];
    uint8_t leaves[8], solids[8], carves[8], uniCarve[8];
    for (int o = 0; o < 8; ++o) {
        refs[o] = child;
        leaves[o] = isLeaf ? 1u : 0u;
        solids[o] = 0u;
        carves[o] = 1u;
        uniCarve[o] = 1u;
    }
    const uint32_t nodeIdx = internNode(refs, leaves, solids, carves, uniCarve);
    uniformCarveNodes_.emplace(key, nodeIdx);
    retain(nodeIdx, false);    // kept for good, as the map is
    return nodeIdx;
}

uint32_t VoxelDagPool::wrapBrickAsNode(uint32_t brickIdx, bool fillAllOctants) {
    if (brickIdx == kDagNullNode || brickIdx >= bricks_.size()) return kDagNullNode;
    uint8_t s, c, u;
    uint32_t co, nm;
    describeRef(brickIdx, true, s, c, u, co, nm);
    uint32_t refs[8];
    uint8_t leaves[8] = {}, solids[8] = {}, carves[8] = {}, uni[8] = {};
    for (int o = 0; o < 8; ++o) {
        if (o == 0 || fillAllOctants) {
            refs[o] = brickIdx;
            leaves[o] = 1u;
            solids[o] = s;
            carves[o] = c;
            uni[o] = u;
        } else {
            refs[o] = kDagNullNode;
        }
    }
    return internNode(refs, leaves, solids, carves, uni);
}

uint32_t VoxelDagPool::registerAnimClip(const std::vector<uint32_t>& frameRoots) {
    DagAnimClip clip;
    clip.frameOffset = static_cast<uint32_t>(animFrames_.size());
    clip.frameCount = static_cast<uint32_t>(frameRoots.size());
    animFrames_.insert(animFrames_.end(), frameRoots.begin(), frameRoots.end());
    // A frame is read as a node or a brick, as it is placed: hold it as both
    // that exists (keeping too much is safe, freeing a frame is not).
    for (uint32_t f : frameRoots) {
        if (f < nodes_.size() && dag_pack::validMaskOf(nodes_[f]) != 0u) retain(f, false);
        if (f < bricks_.size()) {
            const DagBrick& b = bricks_[f];
            if ((b.solidLo | b.solidHi | b.carveLo | b.carveHi) != 0u) retain(f, true);
        }
    }
    const uint32_t idx = static_cast<uint32_t>(animClips_.size());
    animClips_.push_back(clip);
    ++generation_;
    return idx;
}

uint32_t VoxelDagPool::makeAnimChildRef(uint32_t clipIndex, uint32_t phase) const {
    return kDagAnimBit | ((clipIndex & 0x7FFFu) << 16u) | (phase & 0xFFFFu);
}

uint32_t VoxelDagPool::resolveChildRef(uint32_t rawRef, uint32_t animTick) const noexcept {
    if (rawRef == kDagNullNode || (rawRef & kDagAnimBit) == 0u) return rawRef;
    const uint32_t clipIdx = (rawRef >> 16u) & 0x7FFFu;
    const uint32_t phase = rawRef & 0xFFFFu;
    if (clipIdx >= animClips_.size()) return kDagNullNode;
    const DagAnimClip& clip = animClips_[clipIdx];
    if (clip.frameCount == 0u) return kDagNullNode;
    const uint32_t f = (animTick + phase) % clip.frameCount;
    return animFrames_[clip.frameOffset + f];
}

uint32_t VoxelDagPool::setRecursive(
    uint32_t curNode,
    int curDepth,
    int targetDepth,
    int ix,
    int iy,
    int iz,
    uint32_t childRef,
    bool childIsLeaf) {
    const int sh = targetDepth - 1 - curDepth;
    const uint32_t oct = static_cast<uint32_t>(((ix >> sh) & 1) | (((iy >> sh) & 1) << 1) | (((iz >> sh) & 1) << 2));

    uint32_t refs[8];
    uint8_t leaves[8] = {}, solids[8] = {}, carves[8] = {}, uni[8] = {};
    for (int o = 0; o < 8; ++o) refs[o] = kDagNullNode;

    if (curNode != kDagNullNode && curNode < nodes_.size()) {
        const DagNode& n = nodes_[curNode];
        const uint32_t vm = dag_pack::validMaskOf(n);
        const uint32_t sm = dag_pack::solidMaskOf(n);
        const uint32_t cm = dag_pack::carveMaskOf(n);
        const uint32_t lm = dag_pack::leafMaskOf(n);
        uint32_t cursor = 0u;
        for (uint32_t o = 0u; o < 8u; ++o) {
            if ((vm & (1u << o)) == 0u) continue;
            refs[o] = children_[n.childBase + cursor++];
            leaves[o] = (lm & (1u << o)) ? 1u : 0u;
            uint8_t s, c, u;
            uint32_t co, nm;
            describeRef(refs[o], leaves[o] != 0, s, c, u, co, nm);
            solids[o] = (sm & (1u << o)) ? 1u : s;
            carves[o] = (cm & (1u << o)) ? 1u : c;
            uni[o] = u;
        }
    }

    if (curDepth + 1 == targetDepth) {
        refs[oct] = childRef;
        leaves[oct] = childIsLeaf ? 1u : 0u;
        uint8_t s, c, u;
        uint32_t co, nm;
        describeRef(childRef, childIsLeaf, s, c, u, co, nm);
        solids[oct] = s;
        carves[oct] = c;
        uni[oct] = u;
    } else {
        uint32_t existing = refs[oct];
        if (leaves[oct] != 0) existing = promoteBrickToNode(existing);
        const uint32_t updated = setRecursive(existing, curDepth + 1, targetDepth, ix, iy, iz, childRef, childIsLeaf);
        refs[oct] = updated;
        leaves[oct] = 0u;
        uint8_t s, c, u;
        uint32_t co, nm;
        describeRef(updated, false, s, c, u, co, nm);
        solids[oct] = s;
        carves[oct] = c;
        uni[oct] = u;
    }

    return internNode(refs, leaves, solids, carves, uni);
}

uint32_t VoxelDagPool::setSubtreeAt(
    uint32_t rootNode,
    int /*maxDepth*/,
    int targetDepth,
    int ix,
    int iy,
    int iz,
    uint32_t childRef,
    bool childIsLeaf) {
    if (targetDepth <= 0) return childRef;
    return setRecursive(rootNode, 0, targetDepth, ix, iy, iz, childRef, childIsLeaf);
}

uint32_t VoxelDagPool::batchSetRecursive(
    uint32_t curNode,
    int curDepth,
    int targetDepth,
    const DagPlacement* items,
    const uint64_t* keys,
    size_t count) {
    if (count == 0u) return curNode;

    uint32_t refs[8];
    uint8_t leaves[8] = {}, solids[8] = {}, carves[8] = {}, uni[8] = {};
    for (int o = 0; o < 8; ++o) refs[o] = kDagNullNode;

    if (curNode != kDagNullNode && curNode < nodes_.size()) {
        const DagNode& n = nodes_[curNode];
        const uint32_t vm = dag_pack::validMaskOf(n);
        const uint32_t sm = dag_pack::solidMaskOf(n);
        const uint32_t cm = dag_pack::carveMaskOf(n);
        const uint32_t lm = dag_pack::leafMaskOf(n);
        uint32_t cursor = 0u;
        for (uint32_t o = 0u; o < 8u; ++o) {
            if ((vm & (1u << o)) == 0u) continue;
            refs[o] = children_[n.childBase + cursor++];
            leaves[o] = (lm & (1u << o)) ? 1u : 0u;
            uint8_t s, c, u;
            uint32_t co, nm;
            describeRef(refs[o], leaves[o] != 0, s, c, u, co, nm);
            solids[o] = (sm & (1u << o)) ? 1u : s;
            carves[o] = (cm & (1u << o)) ? 1u : c;
            uni[o] = u;
        }
    }

    const int sh = targetDepth - 1 - curDepth;
    if (curDepth + 1 == targetDepth) {
        for (size_t i = 0; i < count; ++i) {
            const DagPlacement& p = items[i];
            const uint32_t oct = static_cast<uint32_t>(
                ((p.ix >> sh) & 1) | (((p.iy >> sh) & 1) << 1) | (((p.iz >> sh) & 1) << 2));
            refs[oct] = p.childRef;
            leaves[oct] = p.isLeaf ? 1u : 0u;
            uint8_t s, c, u;
            uint32_t co, nm;
            describeRef(p.childRef, p.isLeaf, s, c, u, co, nm);
            solids[oct] = s;
            carves[oct] = c;
            uni[oct] = u;
        }
    } else {
        // Sorted by key, so an octant's placements are one run.
        const uint32_t keyShift = static_cast<uint32_t>(3 * sh);
        for (size_t b = 0; b < count;) {
            const uint32_t oct = static_cast<uint32_t>((keys[b] >> keyShift) & 7u);
            size_t e = b + 1;
            while (e < count && static_cast<uint32_t>((keys[e] >> keyShift) & 7u) == oct) ++e;
            uint32_t existing = refs[oct];
            if (leaves[oct] != 0) existing = promoteBrickToNode(existing);
            const uint32_t updated = batchSetRecursive(
                existing, curDepth + 1, targetDepth, items + b, keys + b, e - b);
            b = e;
            refs[oct] = updated;
            leaves[oct] = 0u;
            uint8_t s, c, u;
            uint32_t co, nm;
            describeRef(updated, false, s, c, u, co, nm);
            solids[oct] = s;
            carves[oct] = c;
            uni[oct] = u;
        }
    }

    return internNode(refs, leaves, solids, carves, uni);
}

uint32_t VoxelDagPool::setSubtreesBatch(
    uint32_t rootNode,
    int targetDepth,
    const std::vector<DagPlacement>& placements) {
    if (placements.empty() || targetDepth <= 0) return rootNode;
    // Each cell's octants from the root down, 3 bits a level: sorted, a
    // subtree's placements are one run at every level. Stable, so a cell
    // placed twice keeps the later placement.
    const size_t n = placements.size();
    std::vector<uint64_t> cellKeys(n);
    for (size_t i = 0; i < n; ++i) {
        const DagPlacement& p = placements[i];
        uint64_t k = 0u;
        for (int s = targetDepth - 1; s >= 0; --s) {
            const uint64_t oct = static_cast<uint64_t>(
                ((p.ix >> s) & 1) | (((p.iy >> s) & 1) << 1) | (((p.iz >> s) & 1) << 2));
            k = (k << 3u) | oct;
        }
        cellKeys[i] = k;
    }
    std::vector<uint32_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = static_cast<uint32_t>(i);
    std::stable_sort(order.begin(), order.end(),
                     [&](uint32_t a, uint32_t b) { return cellKeys[a] < cellKeys[b]; });
    std::vector<DagPlacement> items(n);
    std::vector<uint64_t> keys(n);
    for (size_t i = 0; i < n; ++i) {
        items[i] = placements[order[i]];
        keys[i] = cellKeys[order[i]];
    }
    return batchSetRecursive(rootNode, 0, targetDepth, items.data(), keys.data(), n);
}

uint32_t VoxelDagPool::promoteBrickToNode(uint32_t brickIdx) {
    if (brickIdx == kDagNullNode || brickIdx >= bricks_.size()) return kDagNullNode;
    const DagBrick b = bricks_[brickIdx];
    const uint64_t solidBits = (static_cast<uint64_t>(b.solidHi) << 32u) | b.solidLo;
    const uint64_t carveBits = (static_cast<uint64_t>(b.carveHi) << 32u) | b.carveLo;
    const bool uniAttr = (b.flags & kBrickFlagUniformAttr) != 0u;

    uint32_t refs[8];
    uint8_t leaves[8] = {}, solids[8] = {}, carves[8] = {}, uni[8] = {};
    for (int o = 0; o < 8; ++o) {
        refs[o] = kDagNullNode;
        const int ox = (o & 1) << 1;
        const int oy = ((o >> 1) & 1) << 1;
        const int oz = ((o >> 2) & 1) << 1;

        uint64_t subSolid = 0u, subCarve = 0u;
        std::vector<DagVoxelAttr> subAttrs;
        for (int lz = 0; lz < 4; ++lz) {
            for (int ly = 0; ly < 4; ++ly) {
                for (int lx = 0; lx < 4; ++lx) {
                    const int srcBit = (ox + (lx >> 1)) + ((oy + (ly >> 1)) << 2) + ((oz + (lz >> 1)) << 4);
                    const int dstBit = lx + (ly << 2) + (lz << 4);
                    if ((carveBits >> srcBit) & 1ull) subCarve |= (1ull << dstBit);
                    if ((solidBits >> srcBit) & 1ull) {
                        subSolid |= (1ull << dstBit);
                        const uint32_t rank = uniAttr ? 0u : dag_pack::popcount64Below(b.solidLo, b.solidHi, srcBit);
                        subAttrs.push_back(attrs_[b.attrOffset + rank]);
                    }
                }
            }
        }
        const uint32_t childBrick = internBrick(
            static_cast<uint32_t>(subSolid),
            static_cast<uint32_t>(subSolid >> 32u),
            static_cast<uint32_t>(subCarve),
            static_cast<uint32_t>(subCarve >> 32u),
            subAttrs);
        if (childBrick != kDagNullNode) {
            refs[o] = childBrick;
            leaves[o] = 1u;
            uint8_t s, c, u;
            uint32_t co, nm;
            describeRef(childBrick, true, s, c, u, co, nm);
            solids[o] = s;
            carves[o] = c;
            uni[o] = u;
        }
    }
    return internNode(refs, leaves, solids, carves, uni);
}

uint32_t VoxelDagPool::mergeSubtrees(uint32_t baseRef, bool baseIsLeaf, uint32_t overlayRef, bool overlayIsLeaf) {
    if (overlayRef == kDagNullNode) return baseRef;
    if (baseRef == kDagNullNode) return overlayRef;
    if (baseRef == overlayRef && baseIsLeaf == overlayIsLeaf) return baseRef;

    if (baseIsLeaf && overlayIsLeaf) {
        const DagBrick a = bricks_[baseRef];
        const DagBrick b = bricks_[overlayRef];
        const uint64_t aSolid = (static_cast<uint64_t>(a.solidHi) << 32u) | a.solidLo;
        const uint64_t bSolid = (static_cast<uint64_t>(b.solidHi) << 32u) | b.solidLo;
        const uint64_t mSolid = aSolid | bSolid;
        const uint64_t mCarve = ((static_cast<uint64_t>(a.carveHi) << 32u) | a.carveLo) |
                                ((static_cast<uint64_t>(b.carveHi) << 32u) | b.carveLo);

        const bool aUni = (a.flags & kBrickFlagUniformAttr) != 0u;
        const bool bUni = (b.flags & kBrickFlagUniformAttr) != 0u;
        std::vector<DagVoxelAttr> mergedAttrs;
        for (int bit = 0; bit < 64; ++bit) {
            if (((mSolid >> bit) & 1ull) == 0ull) continue;
            if ((bSolid >> bit) & 1ull) {
                const uint32_t r = bUni ? 0u : dag_pack::popcount64Below(b.solidLo, b.solidHi, bit);
                mergedAttrs.push_back(attrs_[b.attrOffset + r]);
            } else {
                const uint32_t r = aUni ? 0u : dag_pack::popcount64Below(a.solidLo, a.solidHi, bit);
                mergedAttrs.push_back(attrs_[a.attrOffset + r]);
            }
        }
        return internBrick(
            static_cast<uint32_t>(mSolid),
            static_cast<uint32_t>(mSolid >> 32u),
            static_cast<uint32_t>(mCarve),
            static_cast<uint32_t>(mCarve >> 32u),
            mergedAttrs);
    }

    if (baseIsLeaf) {
        baseRef = promoteBrickToNode(baseRef);
    }
    if (overlayIsLeaf) {
        overlayRef = promoteBrickToNode(overlayRef);
    }

    if (baseRef >= nodes_.size() || overlayRef >= nodes_.size()) return kDagNullNode;
    const DagNode na = nodes_[baseRef];
    const DagNode nb = nodes_[overlayRef];

    uint32_t refs[8];
    uint8_t leaves[8] = {}, solids[8] = {}, carves[8] = {}, uni[8] = {};
    uint32_t cursorA = 0u, cursorB = 0u;
    const uint32_t vmA = dag_pack::validMaskOf(na);
    const uint32_t vmB = dag_pack::validMaskOf(nb);
    const uint32_t lmA = dag_pack::leafMaskOf(na);
    const uint32_t lmB = dag_pack::leafMaskOf(nb);

    for (uint32_t o = 0u; o < 8u; ++o) {
        const bool hasA = (vmA & (1u << o)) != 0u;
        const bool hasB = (vmB & (1u << o)) != 0u;
        const uint32_t refA = hasA ? children_[na.childBase + cursorA++] : kDagNullNode;
        const uint32_t refB = hasB ? children_[nb.childBase + cursorB++] : kDagNullNode;
        const bool leafA = (lmA & (1u << o)) != 0u;
        const bool leafB = (lmB & (1u << o)) != 0u;

        if (!hasA && !hasB) {
            refs[o] = kDagNullNode;
            continue;
        }
        if (!hasA) {
            refs[o] = refB;
            leaves[o] = leafB ? 1u : 0u;
        } else if (!hasB) {
            refs[o] = refA;
            leaves[o] = leafA ? 1u : 0u;
        } else {
            refs[o] = mergeSubtrees(refA, leafA, refB, leafB);
            leaves[o] = (leafA && leafB) ? 1u : 0u;
        }

        uint8_t s, c, u;
        uint32_t co, nm;
        describeRef(refs[o], leaves[o] != 0u, s, c, u, co, nm);
        solids[o] = s;
        carves[o] = c;
        uni[o] = u;
    }

    return internNode(refs, leaves, solids, carves, uni);
}

uint32_t VoxelDagPool::stampPairAt(
    uint32_t rootNode,
    int maxDepth,
    int targetDepth,
    int ix,
    int iy,
    int iz,
    uint32_t carveRoot,
    uint32_t solidRoot) {
    uint32_t combined = kDagNullNode;
    if (carveRoot != kDagNullNode && solidRoot != kDagNullNode) {
        combined = mergeSubtrees(carveRoot, false, solidRoot, false);
    } else if (carveRoot != kDagNullNode) {
        combined = carveRoot;
    } else {
        combined = solidRoot;
    }
    if (combined == kDagNullNode) return rootNode;
    return setSubtreeAt(rootNode, maxDepth, targetDepth, ix, iy, iz, combined, false);
}

VoxelState VoxelDagPool::queryPoint(
    uint32_t rootNode,
    float ux,
    float uy,
    float uz,
    DagVoxelAttr* outAttr,
    uint32_t animTick) const {
    if (rootNode == kDagNullNode || ux < 0.0f || ux >= 1.0f || uy < 0.0f || uy >= 1.0f || uz < 0.0f || uz >= 1.0f) {
        return VoxelState::Empty;
    }
    uint32_t curRef = resolveChildRef(rootNode, animTick);
    bool curIsLeaf = false;
    float px = ux, py = uy, pz = uz;

    for (int iter = 0; iter < 16; ++iter) {
        if (curRef == kDagNullNode) return VoxelState::Empty;
        if (curIsLeaf) {
            if (curRef >= bricks_.size()) return VoxelState::Empty;
            const DagBrick& b = bricks_[curRef];
            const int lx = std::clamp(static_cast<int>(std::floor(px * 4.0f)), 0, 3);
            const int ly = std::clamp(static_cast<int>(std::floor(py * 4.0f)), 0, 3);
            const int lz = std::clamp(static_cast<int>(std::floor(pz * 4.0f)), 0, 3);
            const int bit = lx + (ly << 2) + (lz << 4);
            const uint64_t solidBits = (static_cast<uint64_t>(b.solidHi) << 32u) | b.solidLo;
            if ((solidBits >> bit) & 1ull) {
                if (outAttr != nullptr) {
                    const uint32_t rank = (b.flags & kBrickFlagUniformAttr)
                        ? 0u
                        : dag_pack::popcount64Below(b.solidLo, b.solidHi, bit);
                    *outAttr = attrs_[b.attrOffset + rank];
                }
                return VoxelState::Solid;
            }
            const uint64_t carveBits = (static_cast<uint64_t>(b.carveHi) << 32u) | b.carveLo;
            return ((carveBits >> bit) & 1ull) ? VoxelState::Carve : VoxelState::Empty;
        }

        if (curRef >= nodes_.size()) return VoxelState::Empty;
        const DagNode& n = nodes_[curRef];
        const int ox = (px >= 0.5f) ? 1 : 0;
        const int oy = (py >= 0.5f) ? 1 : 0;
        const int oz = (pz >= 0.5f) ? 1 : 0;
        const uint32_t oct = static_cast<uint32_t>(ox | (oy << 1) | (oz << 2));
        const uint32_t bit = 1u << oct;
        const uint32_t vm = dag_pack::validMaskOf(n);
        if ((vm & bit) == 0u) return VoxelState::Empty;

        const uint32_t slot = dag_pack::popcount8(vm & (bit - 1u));
        curIsLeaf = (dag_pack::leafMaskOf(n) & bit) != 0u;
        curRef = children_[n.childBase + slot];
        if (!curIsLeaf) curRef = resolveChildRef(curRef, animTick);
        px = px * 2.0f - static_cast<float>(ox);
        py = py * 2.0f - static_cast<float>(oy);
        pz = pz * 2.0f - static_cast<float>(oz);
    }
    return VoxelState::Empty;
}

} // namespace bromesh
