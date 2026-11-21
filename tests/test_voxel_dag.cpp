#include "test_framework.h"
#include "bromesh/voxel/voxel_dag.h"
#include "bromesh/voxel/voxel_chunk.h"
#include "bromesh/primitives/primitives.h"

#include <cmath>
#include <cstdint>
#include <vector>

namespace {

void appendU32(std::vector<uint8_t>& out, uint32_t v) {
    out.push_back(static_cast<uint8_t>(v & 0xFFu));
    out.push_back(static_cast<uint8_t>((v >> 8u) & 0xFFu));
    out.push_back(static_cast<uint8_t>((v >> 16u) & 0xFFu));
    out.push_back(static_cast<uint8_t>((v >> 24u) & 0xFFu));
}

std::vector<uint8_t> buildSyntheticVox() {
    std::vector<uint8_t> vox = {'V', 'O', 'X', ' '};
    appendU32(vox, 150);

    std::vector<uint8_t> children;
    // SIZE chunk (4x4x4)
    children.insert(children.end(), {'S', 'I', 'Z', 'E'});
    appendU32(children, 12);
    appendU32(children, 0);
    appendU32(children, 4);
    appendU32(children, 4);
    appendU32(children, 4);

    // XYZI chunk (2 voxels)
    children.insert(children.end(), {'X', 'Y', 'Z', 'I'});
    appendU32(children, 4 + 2 * 4);
    appendU32(children, 0);
    appendU32(children, 2);
    children.push_back(1);
    children.push_back(1);
    children.push_back(0);
    children.push_back(10);
    children.push_back(2);
    children.push_back(1);
    children.push_back(2);
    children.push_back(25);

    // RGBA chunk (256 * 4 bytes)
    children.insert(children.end(), {'R', 'G', 'B', 'A'});
    appendU32(children, 256 * 4);
    appendU32(children, 0);
    for (int i = 0; i < 256; ++i) {
        if (i == 9) { // palette index 10
            children.push_back(200);
            children.push_back(40);
            children.push_back(40);
            children.push_back(255);
        } else if (i == 24) { // palette index 25
            children.push_back(40);
            children.push_back(180);
            children.push_back(220);
            children.push_back(255);
        } else {
            children.push_back(128);
            children.push_back(128);
            children.push_back(128);
            children.push_back(255);
        }
    }

    vox.insert(vox.end(), {'M', 'A', 'I', 'N'});
    appendU32(vox, 0);
    appendU32(vox, static_cast<uint32_t>(children.size()));
    vox.insert(vox.end(), children.begin(), children.end());
    return vox;
}

} // namespace

TEST(voxel_dag_dedup_and_pointer_instancing) {
    bromesh::VoxelDagPool pool;

    // Build a 16x16x16 sphere archetype via SDF.
    bromesh::DagVoxelizeOptions opts;
    opts.resolution = 16;
    opts.defaultColor[0] = 0.2f;
    opts.defaultColor[1] = 0.7f;
    opts.defaultColor[2] = 0.25f;
    opts.defaultColor[3] = 1.0f;

    int archDepth = 0;
    const uint32_t treeArch = bromesh::voxelizeSdfToDag(
        pool,
        [](float x, float y, float z) {
            return std::sqrt(x * x + y * y + z * z) - 0.85f;
        },
        bromath::AABB3{{-1.0f, -1.0f, -1.0f}, {1.0f, 1.0f, 1.0f}},
        opts,
        &archDepth);
    ASSERT(treeArch != bromesh::kDagNullNode, "archetype DAG built");
    ASSERT(archDepth == 2, "16^3 resolution maps to depth=2 above 4^3 bricks");

    const size_t baseBricks = pool.bricks().size();
    const size_t baseAttrs = pool.attrs().size();

    // Place 1,000 copies into a depth=6 (64x64x64 cell slots) world octree via batch pointer reuse.
    std::vector<bromesh::DagPlacement> batch;
    batch.reserve(1000);
    for (int z = 0; z < 10; ++z) {
        for (int y = 0; y < 10; ++y) {
            for (int x = 0; x < 10; ++x) {
                batch.push_back({x * 2, y, z * 2, treeArch, false});
            }
        }
    }
    uint32_t worldRoot = pool.setSubtreesBatch(bromesh::kDagNullNode, 4, batch);
    ASSERT(batch.size() == 1000, "placed 1000 pointer instances");
    ASSERT(worldRoot != bromesh::kDagNullNode, "world root valid");
    // No new leaf bricks or voxel attrs should have been created!
    ASSERT(pool.bricks().size() == baseBricks, "zero duplicate bricks across 1000 instances");
    ASSERT(pool.attrs().size() == baseAttrs, "zero duplicate voxel attributes across 1000 instances");
    ASSERT(pool.nodes().size() < 300, "hash-consing keeps internal nodes compact across 1000 instances");

    // Incremental single-cell pointer swap must still keep bricks/attrs deduplicated
    const uint32_t singleRoot = pool.setSubtreeAt(worldRoot, 6, 4, 15, 0, 15, treeArch, false);
    ASSERT(singleRoot != bromesh::kDagNullNode && pool.bricks().size() == baseBricks,
           "incremental setSubtreeAt preserves brick deduplication");

    // Verify center of instance (0, 0, 0) in cell (0,0,0) at targetDepth=4 is Solid:
    // cell size = 1/16, center at (1/32, 1/32, 1/32).
    bromesh::DagVoxelAttr attr{};
    const auto st = pool.queryPoint(worldRoot, 1.0f / 32.0f, 1.0f / 32.0f, 1.0f / 32.0f, &attr);
    ASSERT(st == bromesh::VoxelState::Solid, "instance center voxel queries Solid");

    float r, g, b, a;
    bromesh::dag_pack::unpackRgba8(attr.colorOpacity, r, g, b, a);
    ASSERT(g > r && g > b && a > 0.9f, "unpacked green color matches archetype");
}

TEST(voxel_dag_stacked_carve_and_solid_overhang) {
    bromesh::VoxelDagPool pool;

    // Carve cylinder tunnel along Z at y in [-0.5, 0.2], x in [-0.5, 0.5]
    bromesh::DagVoxelizeOptions carveOpts;
    carveOpts.resolution = 16;
    carveOpts.state = bromesh::VoxelState::Carve;
    const uint32_t carveSub = bromesh::voxelizeSdfToDag(
        pool,
        [](float x, float y, float z) {
            (void)z;
            const float dy = y + 0.15f;
            return std::sqrt(x * x + dy * dy) - 0.45f;
        },
        bromath::AABB3{{-1.0f, -1.0f, -1.0f}, {1.0f, 1.0f, 1.0f}},
        carveOpts);
    ASSERT(carveSub != bromesh::kDagNullNode, "carve subtree created");

    // Solid overhanging rock rim above y >= 0.45
    bromesh::DagVoxelizeOptions solidOpts;
    solidOpts.resolution = 16;
    solidOpts.state = bromesh::VoxelState::Solid;
    solidOpts.defaultColor[0] = 0.6f;
    solidOpts.defaultColor[1] = 0.45f;
    solidOpts.defaultColor[2] = 0.3f;
    const uint32_t solidSub = bromesh::voxelizeSdfToDag(
        pool,
        [](float x, float y, float z) {
            (void)z;
            const float dx = std::abs(x) - 0.65f;
            const float dy = std::abs(y - 0.65f) - 0.22f;
            return std::max(dx, dy);
        },
        bromath::AABB3{{-1.0f, -1.0f, -1.0f}, {1.0f, 1.0f, 1.0f}},
        solidOpts);
    ASSERT(solidSub != bromesh::kDagNullNode, "solid overhang subtree created");

    // Stamp both into cell (1, 1, 1) at targetDepth=2 of a depth=4 region -> cell [0.25..0.5)^3
    uint32_t root = pool.stampPairAt(bromesh::kDagNullNode, 4, 2, 1, 1, 1, carveSub, solidSub);
    ASSERT(root != bromesh::kDagNullNode, "stamped pair root valid");

    const auto& rootNode = pool.nodes()[root];
    ASSERT(bromesh::dag_pack::carveMaskOf(rootNode) != 0u, "carveMask propagated to root");
    ASSERT(bromesh::dag_pack::solidMaskOf(rootNode) != 0u, "solidMask propagated to root");

    // Inside cell [0.25..0.5):
    // Cave interior at local (0.5, 0.425, 0.5) -> world (0.375, 0.25 + 0.25*0.425, 0.375)
    ASSERT(pool.queryPoint(root, 0.375f, 0.356f, 0.375f) == bromesh::VoxelState::Carve,
           "cave tunnel interior queries Carve");
    // Overhang at local (0.5, 0.825, 0.5) -> world (0.375, 0.25 + 0.25*0.825, 0.375)
    ASSERT(pool.queryPoint(root, 0.375f, 0.456f, 0.375f) == bromesh::VoxelState::Solid,
           "overhanging rock rim queries Solid");
    // Unstamped octant
    ASSERT(pool.queryPoint(root, 0.8f, 0.8f, 0.8f) == bromesh::VoxelState::Empty,
           "empty octant queries Empty");
}

TEST(voxel_dag_mesh_voxelize_chunk_and_vox_loader) {
    bromesh::VoxelDagPool pool;

    // 1. Mesh voxelizer on sphere
    bromesh::MeshData sph = bromesh::sphere(1.0f, 16, 12);
    bromesh::DagVoxelizeOptions opts;
    opts.resolution = 16;
    opts.defaultColor[0] = 0.9f;
    opts.defaultColor[1] = 0.5f;
    opts.defaultColor[2] = 0.2f;
    const uint32_t sphRoot = bromesh::voxelizeMeshToDag(pool, sph, opts);
    ASSERT(sphRoot != bromesh::kDagNullNode, "sphere mesh voxelized to DAG");
    ASSERT(pool.bricks().size() > 0 && pool.nodes().size() > 0, "sphere DAG has bricks and nodes");

    // Top pole of sphere (~0.5, ~0.95, ~0.5) should have an upward-pointing normal
    bromesh::DagVoxelAttr poleAttr{};
    bool foundTop = false;
    for (int step = 0; step < 16 && !foundTop; ++step) {
        const float uy = 0.96f - static_cast<float>(step) * 0.03f;
        if (pool.queryPoint(sphRoot, 0.52f, uy, 0.52f, &poleAttr) == bromesh::VoxelState::Solid) {
            float nx, ny, nz, rough, emis;
            bromesh::dag_pack::unpackNormalMat(poleAttr.normalMat, nx, ny, nz, rough, emis);
            ASSERT(ny > 0.5f, "top pole voxel normal points upward (+Y)");
            foundTop = true;
        }
    }
    ASSERT(foundTop, "found top pole solid voxel in voxelized sphere");

    // 2. VoxelChunk builder
    bromesh::VoxelChunk chunk(8, 8, 8, 1.0f);
    chunk.setVoxel(4, 4, 4, 1);
    int chunkDepth = 0;
    const uint32_t chunkRoot = bromesh::buildDagFromChunk(pool, chunk, nullptr, 0, &chunkDepth);
    ASSERT(chunkRoot != bromesh::kDagNullNode && chunkDepth == 1, "chunk converted to depth=1 DAG");

    // 3. Synthetic .vox RIFF loader + animation clip
    const std::vector<uint8_t> voxBytes = buildSyntheticVox();
    int voxDepth = 0;
    const uint32_t voxRoot = bromesh::loadVoxToDag(pool, voxBytes.data(), voxBytes.size(), &voxDepth);
    ASSERT(voxRoot != bromesh::kDagNullNode && voxDepth == 1, "synthetic .vox loaded into DAG");

    const uint32_t clip = pool.registerAnimClip({sphRoot, voxRoot});
    const uint32_t animRef = pool.makeAnimChildRef(clip, 1);
    ASSERT((animRef & bromesh::kDagAnimBit) != 0u, "anim child ref sets kDagAnimBit");
    ASSERT(pool.resolveChildRef(animRef, 0) == voxRoot, "phase=1 tick=0 resolves frame 1");
    ASSERT(pool.resolveChildRef(animRef, 1) == sphRoot, "phase=1 tick=1 wraps to frame 0");
}
