#include "test_framework.h"

TEST(greedy_mesh_single_voxel) {
    // 3x3x3 grid, only the center voxel (1,1,1) is solid
    uint8_t voxels[27] = {};
    voxels[1 * 3 * 3 + 1 * 3 + 1] = 1; // z=1, y=1, x=1
    auto m = bromesh::greedyMesh(voxels, 3, 3, 3);
    ASSERT(!m.empty(), "single voxel mesh should be non-empty");
    // 6 faces, each face = 2 triangles = 12 triangles total
    ASSERT(m.triangleCount() == 12, "single voxel should have 12 triangles (6 faces)");
    // 6 faces * 4 vertices = 24 vertices
    ASSERT(m.vertexCount() == 24, "single voxel should have 24 vertices");
    ASSERT(m.hasNormals(), "single voxel mesh should have normals");
    ASSERT(m.hasUVs(), "single voxel mesh should have UVs");
    ASSERT(m.hasColors(), "single voxel mesh should have colors");
}

TEST(greedy_mesh_full_block) {
    // 4x4x4 completely solid grid
    uint8_t voxels[64];
    std::memset(voxels, 1, sizeof(voxels));
    auto m = bromesh::greedyMesh(voxels, 4, 4, 4);
    ASSERT(!m.empty(), "full block mesh should be non-empty");
    // Greedy meshing should merge each face of the cube into a single quad.
    // 6 faces * 2 triangles = 12 triangles
    ASSERT(m.triangleCount() == 12, "full 4x4x4 block should have 12 triangles (6 merged faces)");
    ASSERT(m.vertexCount() == 24, "full 4x4x4 block should have 24 vertices");
}

TEST(greedy_mesh_with_palette) {
    // 2x2x2 grid, all material 2
    uint8_t voxels[8];
    std::memset(voxels, 2, sizeof(voxels));
    // Palette: entry 0 unused, entry 1 = red, entry 2 = green
    float palette[] = {
        0,0,0,0,       // 0: unused
        1,0,0,1,       // 1: red
        0,1,0,1        // 2: green
    };
    auto m = bromesh::greedyMesh(voxels, 2, 2, 2, 1.0f, palette, 3);
    ASSERT(!m.empty(), "palette mesh should be non-empty");
    ASSERT(m.hasColors(), "palette mesh should have colors");
    // All colors should be green (material 2)
    bool allGreen = true;
    for (size_t v = 0; v < m.vertexCount(); ++v) {
        float r = m.colors[v * 4 + 0];
        float g = m.colors[v * 4 + 1];
        float b = m.colors[v * 4 + 2];
        if (std::fabs(r) > 0.001f || std::fabs(g - 1.0f) > 0.001f || std::fabs(b) > 0.001f) {
            allGreen = false;
            break;
        }
    }
    ASSERT(allGreen, "all vertex colors should be green for material 2");
}

