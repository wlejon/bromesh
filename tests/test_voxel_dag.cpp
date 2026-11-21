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

// The brick table grows (incrementally) through several sizes: every brick
// interned before and during a growth is still found.
TEST(voxel_dag_brick_table_grows_without_losing_bricks) {
    bromesh::VoxelDagPool pool;
    const bromesh::DagVoxelAttr attr{};
    const uint32_t n = 200000u;
    std::vector<uint32_t> ids(n);
    for (uint32_t i = 0; i < n; ++i) ids[i] = pool.internBrick(i + 1u, i ^ 0x5a5a5a5au, 0u, 0u, {}, &attr);
    bool same = true, distinct = true;
    for (uint32_t i = 0; i < n; ++i) {
        same = same && pool.internBrick(i + 1u, i ^ 0x5a5a5a5au, 0u, 0u, {}, &attr) == ids[i];
        distinct = distinct && (i == 0 || ids[i] != ids[i - 1]);
    }
    ASSERT(same, "each of 200000 bricks interned again finds its first id");
    ASSERT(distinct, "distinct bricks get distinct ids");
    ASSERT(pool.bricks().size() == n, "re-interning adds no bricks");
}

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

TEST(voxel_dag_collect_reuses_unreached_slots) {
    bromesh::VoxelDagPool pool;
    bromesh::DagVoxelizeOptions opts;
    opts.resolution = 16;
    int archDepth = 0;
    const uint32_t arch = bromesh::voxelizeSdfToDag(
        pool,
        [](float x, float y, float z) { return std::sqrt(x * x + y * y + z * z) - 0.85f; },
        bromath::AABB3{{-1.0f, -1.0f, -1.0f}, {1.0f, 1.0f, 1.0f}},
        opts,
        &archDepth);
    pool.pin(arch, false);

    // A world rebuilt again and again at shifting placements, as a ring
    // recentring on a moving camera is: the pool must stop growing.
    uint32_t root = bromesh::kDagNullNode;
    size_t sizeAfterFirst = 0;
    for (int build = 0; build < 12; ++build) {
        std::vector<bromesh::DagPlacement> batch;
        for (int z = 0; z < 12; ++z) {
            for (int x = 0; x < 12; ++x) {
                if (((x * 7 + z * 3 + build) % 5) == 0) continue;
                batch.push_back({x + build, (x + z + build) & 3, z, arch, false});
            }
        }
        root = pool.setSubtreesBatch(bromesh::kDagNullNode, 5, batch);
        pool.collect({root});
        if (build == 5) sizeAfterFirst = pool.nodes().size() + pool.children().size();
    }
    ASSERT(pool.nodes().size() + pool.children().size() <= sizeAfterFirst + sizeAfterFirst / 10,
           "rebuilt worlds reuse collected slots: the pool plateaus instead of growing");
    ASSERT(pool.liveNodes() < pool.nodes().size(), "the last collect left free node slots");

    // What is placed now still reads back, and the pinned archetype survived.
    const float cell = 1.0f / 32.0f;
    const int lastBuild = 11;
    const float cx = (static_cast<float>(1 + lastBuild) + 0.5f) * cell;
    const float cy = (static_cast<float>((1 + 0 + lastBuild) & 3) + 0.5f) * cell;
    const float cz = 0.5f * cell;
    ASSERT(pool.queryPoint(root, cx, cy, cz) == bromesh::VoxelState::Solid,
           "a placement of the last build queries Solid after collects");
    const uint32_t again = pool.setSubtreesBatch(bromesh::kDagNullNode, 5, {{0, 0, 0, arch, false}});
    ASSERT(pool.queryPoint(again, 0.5f * cell, 0.5f * cell, 0.5f * cell) == bromesh::VoxelState::Solid,
           "the pinned archetype still places after collects");
}

TEST(voxel_dag_release_and_sweep_free_as_they_go) {
    bromesh::VoxelDagPool pool;
    bromesh::DagVoxelizeOptions opts;
    opts.resolution = 16;
    const uint32_t arch = bromesh::voxelizeSdfToDag(
        pool,
        [](float x, float y, float z) { return std::sqrt(x * x + y * y + z * z) - 0.85f; },
        bromath::AABB3{{-1.0f, -1.0f, -1.0f}, {1.0f, 1.0f, 1.0f}},
        opts);
    pool.pin(arch, false);
    pool.sweep();
    const size_t archLive = pool.liveNodes();

    // Three regions rebuilt in turn, each root held, its old root released:
    // no collect, yet the pool keeps to its working size.
    uint32_t roots[3] = {bromesh::kDagNullNode, bromesh::kDagNullNode, bromesh::kDagNullNode};
    size_t sizeAt12 = 0;
    for (int build = 0; build < 24; ++build) {
        std::vector<bromesh::DagPlacement> batch;
        for (int z = 0; z < 16; ++z) {
            for (int x = 0; x < 16; ++x) {
                if (((x * 7 + z * 3 + build) % 5) == 0) continue;
                batch.push_back({x + build % 7, (x + z + build) & 3, z, arch, false});
            }
        }
        // Built in two batches over the old root, as an incremental edit is.
        const int r = build % 3;
        const size_t half = batch.size() / 2;
        uint32_t root = pool.setSubtreesBatch(bromesh::kDagNullNode, 6,
            std::vector<bromesh::DagPlacement>(batch.begin(), batch.begin() + half));
        root = pool.setSubtreesBatch(root, 6, std::vector<bromesh::DagPlacement>(batch.begin() + half, batch.end()));
        pool.retain(root, false);
        pool.release(roots[r], false);
        roots[r] = root;
        pool.sweep();
        // Children runs are reused by length, so they settle a little later.
        if (build == 12) sizeAt12 = pool.nodes().size() + pool.children().size();
    }
    ASSERT(pool.nodes().size() + pool.children().size() <= sizeAt12 + sizeAt12 / 10,
           "released roots and swept spines are reused: the pool plateaus without collect");

    // The counts agree with a full mark: collect finds nothing more to free,
    // and every region still reads back.
    std::vector<float> before;
    const float cell = 1.0f / 64.0f;
    for (uint32_t root : roots)
        for (int z = 0; z < 16; ++z)
            for (int x = 0; x < 24; ++x)
                for (int y = 0; y < 4; ++y)
                    before.push_back(static_cast<float>(pool.queryPoint(
                        root, (x + 0.5f) * cell, (y + 0.5f) * cell, (z + 0.5f) * cell)));
    ASSERT(pool.collect({}) == 0u, "collect frees nothing release and sweep left behind");
    size_t k = 0;
    bool same = true;
    for (uint32_t root : roots)
        for (int z = 0; z < 16; ++z)
            for (int x = 0; x < 24; ++x)
                for (int y = 0; y < 4; ++y)
                    same &= before[k++] == static_cast<float>(pool.queryPoint(
                        root, (x + 0.5f) * cell, (y + 0.5f) * cell, (z + 0.5f) * cell));
    ASSERT(same, "held regions read back the same after a full collect");

    // Releasing every region leaves the archetype alone.
    for (uint32_t& root : roots) {
        pool.release(root, false);
        root = bromesh::kDagNullNode;
    }
    ASSERT(pool.liveNodes() == archLive, "released regions free all they alone reached");
}

TEST(voxel_dag_batch_places_as_single_sets_do) {
    bromesh::VoxelDagPool pool;
    const bromesh::DagVoxelAttr red{bromesh::dag_pack::packRgba8(1.0f, 0.1f, 0.1f, 1.0f),
                                    bromesh::dag_pack::packNormalMat(0.0f, 1.0f, 0.0f, 0.8f)};
    const bromesh::DagVoxelAttr blue{bromesh::dag_pack::packRgba8(0.1f, 0.1f, 1.0f, 1.0f),
                                     bromesh::dag_pack::packNormalMat(0.0f, 1.0f, 0.0f, 0.8f)};
    const uint32_t a = pool.internBrick(0xFFFFFFFFu, 0xFFFFFFFFu, 0u, 0u, {}, &red);
    const uint32_t b = pool.internBrick(0x0000FFFFu, 0u, 0u, 0u, {}, &blue);
    std::vector<bromesh::DagPlacement> batch;
    uint32_t single = bromesh::kDagNullNode;
    for (int i = 0; i < 400; ++i) {
        const int x = (i * 37) % 32, y = (i * 11) % 32, z = (i * 23) % 32;
        const bool useA = (i % 3) != 0;
        batch.push_back({x, y, z, useA ? a : b, true});
        single = pool.setSubtreeAt(single, 5, 5, x, y, z, useA ? a : b, true);
    }
    // A cell placed twice keeps its later placement, in either.
    batch.push_back({1, 2, 3, b, true});
    single = pool.setSubtreeAt(single, 5, 5, 1, 2, 3, b, true);
    batch.push_back({1, 2, 3, a, true});
    single = pool.setSubtreeAt(single, 5, 5, 1, 2, 3, a, true);
    const uint32_t batched = pool.setSubtreesBatch(bromesh::kDagNullNode, 5, batch);
    ASSERT(batched == single, "a sorted batch interns the very tree single sets make");
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

TEST(voxel_dag_mesh_fill_interior_and_split_commit) {
    bromesh::VoxelDagPool pool;
    bromesh::MeshData sph = bromesh::sphere(0.8f, 24, 16);
    bromesh::DagVoxelizeOptions opts;
    opts.resolution = 64;
    opts.unitBounds = true;
    const uint32_t shell = bromesh::voxelizeMeshToDag(pool, sph, opts);
    ASSERT(shell != bromesh::kDagNullNode, "shell voxelized");
    ASSERT(pool.queryPoint(shell, 0.5f, 0.5f, 0.5f) == bromesh::VoxelState::Empty, "a shell is hollow");
    const size_t shellBricks = pool.liveBricks();

    opts.fillInterior = true;
    // Made off the pool (any thread), committed after: the same as in one call.
    bromesh::DagBrickSet set = bromesh::voxelizeMeshToBricks(sph, opts);
    ASSERT(!set.empty() && set.depth() == 4, "brick set made without a pool");
    const uint32_t solid = bromesh::commitVoxelBricks(pool, set);
    ASSERT(solid == bromesh::voxelizeMeshToDag(pool, sph, opts), "split commit interns the same tree");
    ASSERT(pool.queryPoint(solid, 0.5f, 0.5f, 0.5f) == bromesh::VoxelState::Solid, "filled at the centre");
    ASSERT(pool.queryPoint(solid, 0.5f, 0.9f, 0.5f) == bromesh::VoxelState::Solid, "filled below the top");
    ASSERT(pool.queryPoint(solid, 0.02f, 0.02f, 0.02f) == bromesh::VoxelState::Empty, "outside stays empty");
    // The inside's full bricks are one colour, so they are one brick.
    ASSERT(pool.liveBricks() < shellBricks * 2 + 8, "interior bricks dedupe");
}

TEST(voxel_dag_brick_commit_in_steps) {
    bromesh::MeshData sph = bromesh::sphere(1.0f, 16, 12);
    bromesh::DagVoxelizeOptions opts;
    opts.resolution = 64;
    opts.fillInterior = true;

    // A commit a few bricks or nodes a step, a sweep between every two
    // steps: what it has interned survives, and the tree is the one a whole
    // commit makes (the pool interns the same tree to the same root).
    bromesh::VoxelDagPool pool;
    bromesh::DagBrickCommit commit(bromesh::voxelizeMeshToBricks(sph, opts));
    ASSERT(!commit.done() && commit.depth() == 4, "a commit starts undone");
    int steps = 0;
    while (!commit.step(pool, 7)) {
        pool.sweep();
        ++steps;
    }
    ASSERT(steps > 20, "a small budget takes many steps");
    pool.sweep();
    const uint32_t root = commit.root();
    ASSERT(root != bromesh::kDagNullNode, "the steps end at a root");
    const size_t bricks = pool.liveBricks(), nodes = pool.liveNodes();
    const uint32_t whole = bromesh::commitVoxelBricks(pool, bromesh::voxelizeMeshToBricks(sph, opts));
    ASSERT(whole == root, "stepped and whole commits intern one tree");
    ASSERT(pool.liveBricks() == bricks && pool.liveNodes() == nodes, "nothing interned twice");
    ASSERT(pool.queryPoint(root, 0.5f, 0.5f, 0.5f) == bromesh::VoxelState::Solid, "filled at the centre");

    // Held elsewhere, the commit's holds go and the tree stays; released
    // too, nothing holds it and it is freed.
    pool.retain(root, false);
    commit.release(pool);
    pool.sweep();
    ASSERT(pool.liveBricks() == bricks, "the root's hold keeps the tree");
    pool.release(root, false);
    ASSERT(pool.liveBricks() < bricks / 4, "released, the tree is freed");

    // Given up halfway, a commit's holds drop and what it made goes.
    bromesh::VoxelDagPool pool2;
    bromesh::DagBrickCommit half(bromesh::voxelizeMeshToBricks(sph, opts));
    for (int i = 0; i < 30; ++i) half.step(pool2, 7);
    ASSERT(!half.done() && pool2.liveBricks() > 0, "half a commit holds its bricks");
    half.release(pool2);
    pool2.sweep();
    ASSERT(pool2.liveBricks() <= 1 && half.done(), "a commit given up frees what it made");
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
