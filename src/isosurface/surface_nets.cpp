#include "bromesh/isosurface/surface_nets.h"
#include "bromesh/manipulation/normals.h"

#include <cmath>
#include <vector>

namespace bromesh {

MeshData surfaceNets(const float* field, int gridX, int gridY, int gridZ,
                     float isoLevel, float cellSize) {
    MeshData mesh;
    if (!field || gridX < 2 || gridY < 2 || gridZ < 2)
        return mesh;

    // A cell (x,y,z) has corners at grid points (x..x+1, y..y+1, z..z+1).
    // Number of cells in each dimension:
    int cellsX = gridX - 1;
    int cellsY = gridY - 1;
    int cellsZ = gridZ - 1;

    // Map from cell index to vertex index (-1 means no vertex).
    std::vector<int> cellToVertex(cellsX * cellsY * cellsZ, -1);

    auto fieldIdx = [&](int x, int y, int z) -> int {
        return z * gridY * gridX + y * gridX + x;
    };

    auto cellIdx = [&](int x, int y, int z) -> int {
        return z * cellsY * cellsX + y * cellsX + x;
    };

    // 8 corners of a cell at (cx, cy, cz) as offsets from (cx, cy, cz)
    static const int cornerOffsets[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
    };

    // 12 edges of a cube, each connecting two corners
    static const int edges[12][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7}, // X-aligned
        {0, 2}, {1, 3}, {4, 6}, {5, 7}, // Y-aligned
        {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Z-aligned
    };

    // Phase 1: Place vertices
    for (int cz = 0; cz < cellsZ; ++cz) {
        for (int cy = 0; cy < cellsY; ++cy) {
            for (int cx = 0; cx < cellsX; ++cx) {
                // Get corner values
                float cornerVals[8];
                for (int i = 0; i < 8; ++i) {
                    int gx = cx + cornerOffsets[i][0];
                    int gy = cy + cornerOffsets[i][1];
                    int gz = cz + cornerOffsets[i][2];
                    cornerVals[i] = field[fieldIdx(gx, gy, gz)];
                }

                // Check for sign change
                bool hasPositive = false, hasNegative = false;
                for (int i = 0; i < 8; ++i) {
                    if (cornerVals[i] >= isoLevel) hasPositive = true;
                    else hasNegative = true;
                }
                if (!hasPositive || !hasNegative)
                    continue; // No sign change in this cell

                // Find edge crossings and average their positions
                float avgX = 0, avgY = 0, avgZ = 0;
                int crossingCount = 0;

                for (int e = 0; e < 12; ++e) {
                    int c0 = edges[e][0];
                    int c1 = edges[e][1];
                    float v0 = cornerVals[c0];
                    float v1 = cornerVals[c1];

                    // Check if this edge has a sign crossing
                    bool s0 = (v0 >= isoLevel);
                    bool s1 = (v1 >= isoLevel);
                    if (s0 == s1) continue;

                    // Interpolate crossing point
                    float t = (isoLevel - v0) / (v1 - v0);

                    float px = cx + cornerOffsets[c0][0] + t * (cornerOffsets[c1][0] - cornerOffsets[c0][0]);
                    float py = cy + cornerOffsets[c0][1] + t * (cornerOffsets[c1][1] - cornerOffsets[c0][1]);
                    float pz = cz + cornerOffsets[c0][2] + t * (cornerOffsets[c1][2] - cornerOffsets[c0][2]);

                    avgX += px;
                    avgY += py;
                    avgZ += pz;
                    crossingCount++;
                }

                if (crossingCount > 0) {
                    float inv = 1.0f / crossingCount;
                    avgX *= inv;
                    avgY *= inv;
                    avgZ *= inv;

                    int vidx = static_cast<int>(mesh.positions.size() / 3);
                    cellToVertex[cellIdx(cx, cy, cz)] = vidx;
                    mesh.positions.push_back(avgX * cellSize);
                    mesh.positions.push_back(avgY * cellSize);
                    mesh.positions.push_back(avgZ * cellSize);
                }
            }
        }
    }

    if (mesh.positions.empty())
        return mesh;

    // Phase 2: Connect faces
    // For each grid edge with a sign change, the 4 cells sharing that edge form a quad.
    // We iterate over grid vertices and check the 3 edges in +X, +Y, +Z directions.

    for (int gz = 0; gz < gridZ; ++gz) {
        for (int gy = 0; gy < gridY; ++gy) {
            for (int gx = 0; gx < gridX; ++gx) {
                float v0 = field[fieldIdx(gx, gy, gz)];

                // X-edge: (gx,gy,gz) -> (gx+1,gy,gz)
                // Shared by cells: (gx,gy,gz), (gx,gy-1,gz), (gx,gy,gz-1), (gx,gy-1,gz-1)
                if (gx + 1 < gridX) {
                    float v1 = field[fieldIdx(gx + 1, gy, gz)];
                    bool s0 = (v0 >= isoLevel);
                    bool s1 = (v1 >= isoLevel);
                    if (s0 != s1) {
                        // 4 cells sharing this edge
                        int cy0 = gy, cy1 = gy - 1;
                        int cz0 = gz, cz1 = gz - 1;
                        // Cell coords must be valid (0..cells-1)
                        if (cy0 >= 0 && cy0 < cellsY && cy1 >= 0 && cy1 < cellsY &&
                            cz0 >= 0 && cz0 < cellsZ && cz1 >= 0 && cz1 < cellsZ &&
                            gx < cellsX) {
                            int a = cellToVertex[cellIdx(gx, cy0, cz0)];
                            int b = cellToVertex[cellIdx(gx, cy1, cz0)];
                            int c = cellToVertex[cellIdx(gx, cy1, cz1)];
                            int d = cellToVertex[cellIdx(gx, cy0, cz1)];
                            if (a >= 0 && b >= 0 && c >= 0 && d >= 0) {
                                // Winding order based on sign
                                if (s0) {
                                    // v0 >= isoLevel: normal should point in +X
                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(b);
                                    mesh.indices.push_back(c);

                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(c);
                                    mesh.indices.push_back(d);
                                } else {
                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(c);
                                    mesh.indices.push_back(b);

                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(d);
                                    mesh.indices.push_back(c);
                                }
                            }
                        }
                    }
                }

                // Y-edge: (gx,gy,gz) -> (gx,gy+1,gz)
                // Shared by cells: (gx,gy,gz), (gx-1,gy,gz), (gx,gy,gz-1), (gx-1,gy,gz-1)
                if (gy + 1 < gridY) {
                    float v1 = field[fieldIdx(gx, gy + 1, gz)];
                    bool s0 = (v0 >= isoLevel);
                    bool s1 = (v1 >= isoLevel);
                    if (s0 != s1) {
                        int cx0 = gx, cx1 = gx - 1;
                        int cz0 = gz, cz1 = gz - 1;
                        if (cx0 >= 0 && cx0 < cellsX && cx1 >= 0 && cx1 < cellsX &&
                            cz0 >= 0 && cz0 < cellsZ && cz1 >= 0 && cz1 < cellsZ &&
                            gy < cellsY) {
                            int a = cellToVertex[cellIdx(cx0, gy, cz0)];
                            int b = cellToVertex[cellIdx(cx1, gy, cz0)];
                            int c = cellToVertex[cellIdx(cx1, gy, cz1)];
                            int d = cellToVertex[cellIdx(cx0, gy, cz1)];
                            if (a >= 0 && b >= 0 && c >= 0 && d >= 0) {
                                if (s0) {
                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(c);
                                    mesh.indices.push_back(b);

                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(d);
                                    mesh.indices.push_back(c);
                                } else {
                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(b);
                                    mesh.indices.push_back(c);

                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(c);
                                    mesh.indices.push_back(d);
                                }
                            }
                        }
                    }
                }

                // Z-edge: (gx,gy,gz) -> (gx,gy,gz+1)
                // Shared by cells: (gx,gy,gz), (gx-1,gy,gz), (gx,gy-1,gz), (gx-1,gy-1,gz)
                if (gz + 1 < gridZ) {
                    float v1 = field[fieldIdx(gx, gy, gz + 1)];
                    bool s0 = (v0 >= isoLevel);
                    bool s1 = (v1 >= isoLevel);
                    if (s0 != s1) {
                        int cx0 = gx, cx1 = gx - 1;
                        int cy0 = gy, cy1 = gy - 1;
                        if (cx0 >= 0 && cx0 < cellsX && cx1 >= 0 && cx1 < cellsX &&
                            cy0 >= 0 && cy0 < cellsY && cy1 >= 0 && cy1 < cellsY &&
                            gz < cellsZ) {
                            int a = cellToVertex[cellIdx(cx0, cy0, gz)];
                            int b = cellToVertex[cellIdx(cx1, cy0, gz)];
                            int c = cellToVertex[cellIdx(cx1, cy1, gz)];
                            int d = cellToVertex[cellIdx(cx0, cy1, gz)];
                            if (a >= 0 && b >= 0 && c >= 0 && d >= 0) {
                                if (s0) {
                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(b);
                                    mesh.indices.push_back(c);

                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(c);
                                    mesh.indices.push_back(d);
                                } else {
                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(c);
                                    mesh.indices.push_back(b);

                                    mesh.indices.push_back(a);
                                    mesh.indices.push_back(d);
                                    mesh.indices.push_back(c);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Compute normals
    if (!mesh.indices.empty()) {
        computeNormals(mesh);
    }

    return mesh;
}

} // namespace bromesh
