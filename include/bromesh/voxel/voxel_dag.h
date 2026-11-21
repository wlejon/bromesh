#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/voxel/voxel_chunk.h"
#include <bromath/aabb.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

namespace bromesh {

enum class VoxelState : uint8_t {
    Empty = 0,
    Solid = 1,
    Carve = 2,
};

inline constexpr uint32_t kDagNullNode = 0xFFFFFFFFu;
inline constexpr uint32_t kDagAnimBit = 0x80000000u;

inline constexpr uint32_t kDagFlagUniformCarve = 1u << 0;
inline constexpr uint32_t kDagFlagEmissive = 1u << 1;

inline constexpr uint32_t kBrickFlagUniformAttr = 1u << 0;
inline constexpr uint32_t kBrickFlagUniformCarve = 1u << 1;

/// Internal 2x2x2 octree node stored flat in GPU SSBO (`16` bytes = `vec4<u32>`).
/// `masks`:
///   bits 0..7   = solidMask (octant contains Solid voxels)
///   bits 8..15  = carveMask (octant contains Carve voxels)
///   bits 16..23 = leafMask  (valid octant indexes `bricks_`, else `nodes_` or anim clip)
///   bits 24..31 = flags     (`kDagFlag*`)
struct DagNode {
    uint32_t masks = 0;
    uint32_t childBase = 0;
    uint32_t packedColorOpacity = 0; // RGBA8: RGB albedo, A = coverage/opacity (0..255)
    uint32_t packedNormalMat = 0;    // octX:8 | octY:8 | roughness:8 | emissiveOrAo:8
};
static_assert(sizeof(DagNode) == 16, "DagNode must be 16 bytes for GPU vec4<u32> layout");

/// Leaf 4x4x4 brick stored in GPU SSBO (`32` bytes = 2 x `vec4<u32>`).
/// Voxel (lx, ly, lz) in [0..3]^3 maps to bit `lx + (ly << 2) + (lz << 4)`.
struct DagBrick {
    uint32_t solidLo = 0;
    uint32_t solidHi = 0;
    uint32_t carveLo = 0;
    uint32_t carveHi = 0;
    uint32_t attrOffset = 0;
    uint32_t packedColorOpacity = 0;
    uint32_t packedNormalMat = 0;
    uint32_t flags = 0;
};
static_assert(sizeof(DagBrick) == 32, "DagBrick must be 32 bytes for GPU layout");

/// Per-solid-voxel attribute (`8` bytes = `vec2<u32>`).
struct DagVoxelAttr {
    uint32_t colorOpacity = 0;
    uint32_t normalMat = 0;
};
static_assert(sizeof(DagVoxelAttr) == 8, "DagVoxelAttr must be 8 bytes");

struct DagPlacement {
    int ix = 0;
    int iy = 0;
    int iz = 0;
    uint32_t childRef = kDagNullNode;
    bool isLeaf = false;
};

/// Cyclic animated subtree descriptor (`16` bytes = `vec4<u32>`).
struct DagAnimClip {
    uint32_t frameOffset = 0;
    uint32_t frameCount = 0;
    uint32_t flags = 0;
    uint32_t reserved = 0;
};
static_assert(sizeof(DagAnimClip) == 16, "DagAnimClip must be 16 bytes");

struct DagVoxelizeOptions {
    int resolution = 32; // Power of 2 in [8, 256]
    float defaultColor[4] = {0.72f, 0.72f, 0.72f, 1.0f};
    float roughness = 0.8f;
    float emissive = 0.0f;
    bool fillInterior = false;
    VoxelState state = VoxelState::Solid;
    bool fitToCube = true;
    bool unitBounds = false; // When true, map [-1, +1]^3 directly onto [0, res]^3 (matches voxelizeSdfToDag).
};

namespace dag_pack {

uint32_t packRgba8(float r, float g, float b, float a);
void unpackRgba8(uint32_t packed, float& r, float& g, float& b, float& a);

uint16_t packOctNormal(float nx, float ny, float nz);
void unpackOctNormal(uint16_t packed, float& nx, float& ny, float& nz);

uint32_t packNormalMat(float nx, float ny, float nz, float roughness, float emissive = 0.0f);
void unpackNormalMat(uint32_t packed, float& nx, float& ny, float& nz, float& roughness, float& emissive);

inline uint32_t popcount8(uint32_t v) noexcept {
    uint32_t c = 0;
    v &= 0xFFu;
    while (v) {
        v &= (v - 1u);
        ++c;
    }
    return c;
}

inline uint32_t popcount32(uint32_t v) noexcept {
    v = v - ((v >> 1u) & 0x55555555u);
    v = (v & 0x33333333u) + ((v >> 2u) & 0x33333333u);
    return (((v + (v >> 4u)) & 0x0F0F0F0Fu) * 0x01010101u) >> 24u;
}

inline uint32_t popcount64Below(uint32_t lo, uint32_t hi, int bit) noexcept {
    if (bit <= 0) return 0u;
    if (bit < 32) {
        return popcount32(lo & ((1u << static_cast<uint32_t>(bit)) - 1u));
    }
    const uint32_t hiMask = (bit >= 64) ? 0xFFFFFFFFu : ((1u << static_cast<uint32_t>(bit - 32)) - 1u);
    return popcount32(lo) + popcount32(hi & hiMask);
}

inline uint32_t solidMaskOf(const DagNode& n) noexcept { return n.masks & 0xFFu; }
inline uint32_t carveMaskOf(const DagNode& n) noexcept { return (n.masks >> 8u) & 0xFFu; }
inline uint32_t leafMaskOf(const DagNode& n) noexcept { return (n.masks >> 16u) & 0xFFu; }
inline uint32_t flagsOf(const DagNode& n) noexcept { return (n.masks >> 24u) & 0xFFu; }
inline uint32_t validMaskOf(const DagNode& n) noexcept { return (n.masks & 0xFFu) | ((n.masks >> 8u) & 0xFFu); }

inline uint32_t makeNodeMasks(uint32_t solid, uint32_t carve, uint32_t leaf, uint32_t flags) noexcept {
    return (solid & 0xFFu) | ((carve & 0xFFu) << 8u) | ((leaf & 0xFFu) << 16u) | ((flags & 0xFFu) << 24u);
}

} // namespace dag_pack

/// Content-addressed Sparse Voxel DAG storage pool supporting:
///   - hash-consed brick & internal node deduplication across identical subtrees
///   - stacked additive (`Solid`) and subtractive (`Carve`) bitmasks in one hierarchy
///   - O(depth) path-copied pointer swaps (`setSubtreeAt`, `stampPairAt`)
///   - O(1) cyclic animation via frame-root clips (`kDagAnimBit` child pointers)
class VoxelDagPool {
public:
    VoxelDagPool() = default;

    void clear();

    /// Intern a 4x4x4 leaf brick. `solidAttrs` has `popcount64(solidLo, solidHi)`
    /// entries in bit order (or empty if all solid voxels share `uniformAttr`).
    uint32_t internBrick(
        uint32_t solidLo,
        uint32_t solidHi,
        uint32_t carveLo,
        uint32_t carveHi,
        const std::vector<DagVoxelAttr>& solidAttrs,
        const DagVoxelAttr* uniformAttr = nullptr);

    /// Intern an internal 2x2x2 octree node given its 8 octant entries (`childRefs[o] == kDagNullNode` if empty).
    uint32_t internNode(
        const uint32_t childRefs[8],
        const uint8_t childIsLeaf[8],
        const uint8_t childHasSolid[8],
        const uint8_t childHasCarve[8],
        const uint8_t childUniformCarve[8]);

    /// Intern a pure uniform-carve internal subtree of depth `depthRemaining` (>= 1).
    uint32_t internUniformCarveNode(int depthRemaining);

    /// Wrap a single leaf brick into a depth=1 DagNode located in octant 0 (or replicated into all 8 octants).
    uint32_t wrapBrickAsNode(uint32_t brickIdx, bool fillAllOctants = false);

    /// Register a cyclic animation clip from frame `DagNode` root indices.
    /// Returns clip index for `makeAnimChildRef`.
    uint32_t registerAnimClip(const std::vector<uint32_t>& frameRoots);

    /// Construct a tagged child reference selecting `clipIndex` with starting frame offset `phase`.
    uint32_t makeAnimChildRef(uint32_t clipIndex, uint32_t phase = 0) const;

    /// Path-copy `rootNode` replacing octree cell `(ix, iy, iz)` at `targetDepth` (1 <= targetDepth <= maxDepth)
    /// with `childRef` (or clearing the cell when `childRef == kDagNullNode`).
    /// Returns the new root node index (or `kDagNullNode` if the whole tree became empty).
    uint32_t setSubtreeAt(
        uint32_t rootNode,
        int maxDepth,
        int targetDepth,
        int ix,
        int iy,
        int iz,
        uint32_t childRef,
        bool childIsLeaf = false);

    /// Apply a batch of cell updates at `targetDepth` in a single top-down deduplicated pass
    /// without creating intermediate transient spine nodes.
    uint32_t setSubtreesBatch(
        uint32_t rootNode,
        int targetDepth,
        const std::vector<DagPlacement>& placements);

    /// Stamp both a subtractive `Carve` subtree and an additive `Solid` subtree into cell `(ix, iy, iz)`
    /// at `targetDepth`, merging them cleanly so overhangs and carved voids coexist in the same octants.
    uint32_t stampPairAt(
        uint32_t rootNode,
        int maxDepth,
        int targetDepth,
        int ix,
        int iy,
        int iz,
        uint32_t carveRoot,
        uint32_t solidRoot);

    /// Recursively merge `overlay` onto `base` (either nodes or leaf bricks), keeping shared untouched subtrees
    /// in O(1). Where both have solid voxels, `overlay` wins; `Carve` masks union.
    uint32_t mergeSubtrees(uint32_t baseRef, bool baseIsLeaf, uint32_t overlayRef, bool overlayIsLeaf);

    /// CPU reference query at normalized coordinates `[0, 1)^3` within `rootNode` (`depth = maxDepth`,
    /// where `rootNode` may resolve animated references at `animTick`).
    VoxelState queryPoint(
        uint32_t rootNode,
        float ux,
        float uy,
        float uz,
        DagVoxelAttr* outAttr = nullptr,
        uint32_t animTick = 0) const;

    uint32_t resolveChildRef(uint32_t rawRef, uint32_t animTick) const noexcept;

    /// Hold `ref` (a node, or a brick when `isLeaf`): it and all it reaches
    /// stay until it is released as many times (a region's root, an archetype).
    void retain(uint32_t ref, bool isLeaf);

    /// Drop a hold. What then nothing holds and no node reaches is freed at
    /// once, so replacing a root costs what it alone reached.
    void release(uint32_t ref, bool isLeaf);

    /// Keep `ref` for good: `retain` never released.
    void pin(uint32_t ref, bool isLeaf) { retain(ref, isLeaf); }

    /// Free what interning made since the last sweep that nothing holds and
    /// no node reaches (the spines an edit replaced, a result never kept).
    /// Call after an edit has its result held. Returns nodes and bricks freed.
    size_t sweep();

    /// Free what is reachable from none of `rootNodes`, the held refs, the
    /// animation clips' frames, whatever the holds say (a full mark and
    /// sweep; `release` and `sweep` already free as they go). Slots are
    /// reused by later interning, so every index still reached stays valid
    /// (nothing moves). Returns the nodes and bricks freed.
    size_t collect(const std::vector<uint32_t>& rootNodes);

    /// The arrays a copy of the pool (a GPU buffer) mirrors, by page.
    enum class Array : uint32_t { Nodes, Children, Bricks, Attrs };
    static constexpr uint32_t kPageElems = 4096;
    /// The generation each page (`kPageElems` elements) of an array was last
    /// written at: a copy made at generation g needs only the pages newer
    /// than g (and the array's growth).
    const std::vector<uint64_t>& pageGenerations(Array a) const noexcept {
        return pageGen_[static_cast<uint32_t>(a)];
    }

    /// Nodes and bricks in use (the arrays' sizes less their free slots).
    size_t liveNodes() const noexcept { return nodes_.size() - freeNodes_.size(); }
    size_t liveBricks() const noexcept { return bricks_.size() - freeBricks_.size(); }

    const std::vector<DagNode>& nodes() const noexcept { return nodes_; }
    const std::vector<uint32_t>& children() const noexcept { return children_; }
    const std::vector<DagBrick>& bricks() const noexcept { return bricks_; }
    const std::vector<DagVoxelAttr>& attrs() const noexcept { return attrs_; }
    const std::vector<DagAnimClip>& animClips() const noexcept { return animClips_; }
    const std::vector<uint32_t>& animFrames() const noexcept { return animFrames_; }

    uint64_t generation() const noexcept { return generation_; }

private:
    struct BrickKey {
        uint32_t solidLo = 0;
        uint32_t solidHi = 0;
        uint32_t carveLo = 0;
        uint32_t carveHi = 0;
        uint32_t packedColorOpacity = 0;
        uint32_t packedNormalMat = 0;
        uint32_t flags = 0;
        uint32_t attrHash = 0;
        std::vector<DagVoxelAttr> inlineAttrs;
        bool operator==(const BrickKey& o) const noexcept;
    };
    struct BrickKeyHash {
        size_t operator()(const BrickKey& k) const noexcept;
    };

    struct NodeKey {
        uint32_t masks = 0;
        uint32_t count = 0;
        uint32_t childRefs[8] = {};
        bool operator==(const NodeKey& o) const noexcept;
    };
    struct NodeKeyHash {
        size_t operator()(const NodeKey& k) const noexcept;
    };

    uint32_t promoteBrickToNode(uint32_t brickIdx);
    void describeRef(
        uint32_t rawRef,
        bool isLeaf,
        uint8_t& outHasSolid,
        uint8_t& outHasCarve,
        uint8_t& outUniformCarve,
        uint32_t& outColorOpacity,
        uint32_t& outNormalMat) const;
    uint32_t setRecursive(
        uint32_t curNode,
        int curDepth,
        int targetDepth,
        int ix,
        int iy,
        int iz,
        uint32_t childRef,
        bool childIsLeaf);
    // `items` sorted by `keys`, the cell's octants from the root down (3 bits
    // a level), so each octant's placements are one run at every level.
    uint32_t batchSetRecursive(
        uint32_t curNode,
        int curDepth,
        int targetDepth,
        const DagPlacement* items,
        const uint64_t* keys,
        size_t count);

    std::vector<DagNode> nodes_;
    std::vector<uint32_t> children_;
    std::vector<DagBrick> bricks_;
    std::vector<DagVoxelAttr> attrs_;
    std::vector<DagAnimClip> animClips_;
    std::vector<uint32_t> animFrames_;

    std::unordered_map<BrickKey, uint32_t, BrickKeyHash> brickMap_;
    std::unordered_map<NodeKey, uint32_t, NodeKeyHash> nodeMap_;
    std::unordered_map<uint64_t, uint32_t> uniformCarveNodes_;
    uint64_t generation_ = 1;

    // What was freed, reused by interning: node and brick slots, and runs of
    // `children_` (by length 1..8) and `attrs_` (by length 1..64).
    uint32_t allocNode();
    uint32_t allocBrick();
    uint32_t allocChildren(uint32_t count);
    uint32_t allocAttrs(uint32_t count);
    std::vector<uint32_t> freeNodes_;
    std::vector<uint32_t> freeBricks_;
    std::vector<std::vector<uint32_t>> freeChildren_;
    std::vector<std::vector<uint32_t>> freeAttrs_;

    // References to each node and brick: from the nodes that hold it, and
    // the outside holds (`retain`), which are also kept by ref in
    // `heldNodes_` / `heldBricks_` for `collect`. A slot is free when its
    // node has no children (its brick no bits).
    std::vector<uint32_t> nodeRefs_;
    std::vector<uint32_t> brickRefs_;
    std::unordered_map<uint32_t, uint32_t> heldNodes_;
    std::unordered_map<uint32_t, uint32_t> heldBricks_;
    // Made since the last sweep.
    std::vector<uint32_t> freshNodes_;
    std::vector<uint32_t> freshBricks_;
    size_t freeUnreferenced(uint32_t ref, bool isLeaf);
    void freeNodeSlot(uint32_t i);
    void freeBrickSlot(uint32_t i);
    void recountRefs();

    std::vector<uint64_t> pageGen_[4];
    void touch(Array a, uint32_t first, uint32_t count);
};

/// Build a `VoxelDagPool` subtree from a dense `VoxelChunk`.
uint32_t buildDagFromChunk(
    VoxelDagPool& pool,
    const VoxelChunk& chunk,
    const float* palette = nullptr,
    int paletteCount = 0,
    int* outDepth = nullptr);

/// Voxelize an indexed triangle mesh into a deduplicated `VoxelDagPool` subtree.
uint32_t voxelizeMeshToDag(
    VoxelDagPool& pool,
    const MeshData& mesh,
    const DagVoxelizeOptions& opts = {},
    int* outDepth = nullptr);

/// Voxelize an implicit SDF field (`sdf(x, y, z) <= 0` inside) over `bounds` into a `VoxelDagPool` subtree.
uint32_t voxelizeSdfToDag(
    VoxelDagPool& pool,
    const std::function<float(float, float, float)>& sdfFn,
    const bromath::AABB3& bounds,
    const DagVoxelizeOptions& opts = {},
    int* outDepth = nullptr);

/// Parse an in-memory MagicaVoxel `.vox` binary buffer into a `VoxelDagPool` subtree (Y-up aligned).
uint32_t loadVoxToDag(
    VoxelDagPool& pool,
    const uint8_t* data,
    size_t size,
    int* outDepth = nullptr,
    float roughness = 0.75f);

} // namespace bromesh
