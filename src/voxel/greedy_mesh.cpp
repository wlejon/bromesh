#include "bromesh/voxel/greedy_mesh.h"

#include <cstring>
#include <vector>

namespace bromesh {

MeshData greedyMesh(const uint8_t* voxels, int gridX, int gridY, int gridZ,
                    float cellSize, const float* palette, int paletteCount) {
    MeshData mesh;
    if (!voxels || gridX <= 0 || gridY <= 0 || gridZ <= 0) return mesh;

    // Helper to sample voxel safely (returns 0 for out-of-bounds)
    auto voxelAt = [&](int x, int y, int z) -> uint8_t {
        if (x < 0 || x >= gridX || y < 0 || y >= gridY || z < 0 || z >= gridZ)
            return 0;
        return voxels[z * gridY * gridX + y * gridX + x];
    };

    // Emit a quad as 2 triangles (4 verts, 6 indices)
    // p0--p1
    // |  / |
    // p2--p3
    // But we need correct winding for the normal direction.
    auto emitQuad = [&](const float p[4][3], const float n[3],
                        float uWidth, float vHeight, uint8_t matId) {
        uint32_t base = static_cast<uint32_t>(mesh.positions.size() / 3);

        for (int i = 0; i < 4; ++i) {
            mesh.positions.push_back(p[i][0]);
            mesh.positions.push_back(p[i][1]);
            mesh.positions.push_back(p[i][2]);
            mesh.normals.push_back(n[0]);
            mesh.normals.push_back(n[1]);
            mesh.normals.push_back(n[2]);
        }

        // UVs: map quad dimensions
        mesh.uvs.push_back(0.0f);       mesh.uvs.push_back(0.0f);
        mesh.uvs.push_back(uWidth);     mesh.uvs.push_back(0.0f);
        mesh.uvs.push_back(0.0f);       mesh.uvs.push_back(vHeight);
        mesh.uvs.push_back(uWidth);     mesh.uvs.push_back(vHeight);

        // Colors from palette
        float r = 1.0f, g = 1.0f, b = 1.0f, a = 1.0f;
        if (palette && paletteCount > 0 && matId < paletteCount) {
            r = palette[matId * 4 + 0];
            g = palette[matId * 4 + 1];
            b = palette[matId * 4 + 2];
            a = palette[matId * 4 + 3];
        }
        for (int i = 0; i < 4; ++i) {
            mesh.colors.push_back(r);
            mesh.colors.push_back(g);
            mesh.colors.push_back(b);
            mesh.colors.push_back(a);
        }

        // Two triangles: 0-1-2, 1-3-2
        mesh.indices.push_back(base + 0);
        mesh.indices.push_back(base + 1);
        mesh.indices.push_back(base + 2);
        mesh.indices.push_back(base + 1);
        mesh.indices.push_back(base + 3);
        mesh.indices.push_back(base + 2);
    };

    // For each axis (0=X, 1=Y, 2=Z) and direction (0=negative, 1=positive)
    for (int axis = 0; axis < 3; ++axis) {
        // u and v are the two axes perpendicular to the sweep axis
        int u = (axis + 1) % 3;
        int v = (axis + 2) % 3;

        int dims[3] = { gridX, gridY, gridZ };
        int dimAxis = dims[axis];
        int dimU = dims[u];
        int dimV = dims[v];

        for (int dir = 0; dir < 2; ++dir) {
            // Normal: points in +axis or -axis direction
            float normal[3] = { 0, 0, 0 };
            normal[axis] = (dir == 1) ? 1.0f : -1.0f;

            // Sweep slices perpendicular to axis
            std::vector<uint8_t> mask(static_cast<size_t>(dimU) * dimV, 0);

            for (int slice = 0; slice < dimAxis; ++slice) {
                // Build mask for this slice
                std::memset(mask.data(), 0, mask.size());

                for (int iv = 0; iv < dimV; ++iv) {
                    for (int iu = 0; iu < dimU; ++iu) {
                        int pos[3] = { 0, 0, 0 };
                        pos[axis] = slice;
                        pos[u] = iu;
                        pos[v] = iv;

                        uint8_t cur = voxelAt(pos[0], pos[1], pos[2]);

                        // Neighbor in the direction we're checking
                        int npos[3] = { pos[0], pos[1], pos[2] };
                        npos[axis] += (dir == 1) ? 1 : -1;
                        uint8_t neighbor = voxelAt(npos[0], npos[1], npos[2]);

                        // Face is exposed if current is solid and neighbor is empty
                        if (cur != 0 && neighbor == 0) {
                            mask[static_cast<size_t>(iv) * dimU + iu] = cur;
                        }
                    }
                }

                // Greedy merge the mask into quads
                for (int iv = 0; iv < dimV; ++iv) {
                    for (int iu = 0; iu < dimU; ) {
                        uint8_t mat = mask[static_cast<size_t>(iv) * dimU + iu];
                        if (mat == 0) {
                            ++iu;
                            continue;
                        }

                        // Find width (run in U direction)
                        int w = 1;
                        while (iu + w < dimU &&
                               mask[static_cast<size_t>(iv) * dimU + iu + w] == mat) {
                            ++w;
                        }

                        // Find height (extend in V direction)
                        int h = 1;
                        bool done = false;
                        while (iv + h < dimV && !done) {
                            for (int k = 0; k < w; ++k) {
                                if (mask[static_cast<size_t>(iv + h) * dimU + iu + k] != mat) {
                                    done = true;
                                    break;
                                }
                            }
                            if (!done) ++h;
                        }

                        // Clear mask for merged region
                        for (int jj = 0; jj < h; ++jj) {
                            for (int ii = 0; ii < w; ++ii) {
                                mask[static_cast<size_t>(iv + jj) * dimU + iu + ii] = 0;
                            }
                        }

                        // Emit quad
                        // The quad spans [iu, iu+w] x [iv, iv+h] on the slice plane
                        // Position along the axis:
                        // For positive face (dir==1): the face is at slice+1
                        // For negative face (dir==0): the face is at slice
                        float slicePos = (dir == 1) ? (slice + 1.0f) : (float)slice;

                        // Build 4 corner positions
                        float p[4][3];
                        for (int c = 0; c < 4; ++c) {
                            float cu = (c & 1) ? (float)(iu + w) : (float)iu;
                            float cv = (c & 2) ? (float)(iv + h) : (float)iv;
                            p[c][axis] = slicePos * cellSize;
                            p[c][u] = cu * cellSize;
                            p[c][v] = cv * cellSize;
                        }

                        // Ensure correct winding order for the face normal.
                        // The default winding (0,1,2 / 1,3,2) produces a normal
                        // via cross(p1-p0, p2-p0). We need this to match our
                        // intended normal direction. If it doesn't, swap p1 and p2
                        // (equivalently, swap the two triangle windings).
                        // Cross product of edge01 x edge02:
                        float e1[3] = { p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2] };
                        float e2[3] = { p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2] };
                        float cx = e1[1]*e2[2] - e1[2]*e2[1];
                        float cy = e1[2]*e2[0] - e1[0]*e2[2];
                        float cz = e1[0]*e2[1] - e1[1]*e2[0];
                        float dot = cx*normal[0] + cy*normal[1] + cz*normal[2];
                        if (dot < 0) {
                            // Swap p1 and p2 to flip winding
                            float tmp[3];
                            tmp[0] = p[1][0]; tmp[1] = p[1][1]; tmp[2] = p[1][2];
                            p[1][0] = p[2][0]; p[1][1] = p[2][1]; p[1][2] = p[2][2];
                            p[2][0] = tmp[0]; p[2][1] = tmp[1]; p[2][2] = tmp[2];
                            // Also swap p3 corners accordingly -- actually we need
                            // a full re-layout. Swapping p1<->p2 means:
                            // p0=(iu,iv), p1=(iu,iv+h), p2=(iu+w,iv), p3=(iu+w,iv+h)
                            // But we also need p3 consistent. Let's just swap and fix p3.
                            // Actually, the quad is p0-p1-p2-p3 where triangles are
                            // 0-1-2 and 1-3-2. After swapping p1<->p2:
                            // tri1: 0-p2_old-p1_old, tri2: p2_old-p3-p1_old
                            // This should flip the normal. p3 stays the same.
                        }

                        emitQuad(p, normal, (float)w, (float)h, mat);

                        iu += w;
                    }
                }
            }
        }
    }

    return mesh;
}

} // namespace bromesh
