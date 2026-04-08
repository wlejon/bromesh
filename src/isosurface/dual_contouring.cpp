#include "bromesh/isosurface/dual_contouring.h"
#include "bromesh/manipulation/normals.h"

#include <cmath>
#include <vector>

namespace bromesh {

// Solve a 3x3 symmetric positive-definite system (A^T A + lambda*I)x = A^T b
// using Cramer's rule. Returns false if the determinant is near zero.
static bool solve3x3(const float m[3][3], const float rhs[3], float out[3]) {
    // Compute determinant
    float det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
              - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
              + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if (std::fabs(det) < 1e-10f)
        return false;

    float invDet = 1.0f / det;

    // Cramer's rule
    out[0] = (rhs[0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (rhs[1] * m[2][2] - m[1][2] * rhs[2])
            + m[0][2] * (rhs[1] * m[2][1] - m[1][1] * rhs[2])) * invDet;

    out[1] = (m[0][0] * (rhs[1] * m[2][2] - m[1][2] * rhs[2])
            - rhs[0] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * rhs[2] - rhs[1] * m[2][0])) * invDet;

    out[2] = (m[0][0] * (m[1][1] * rhs[2] - rhs[1] * m[2][1])
            - m[0][1] * (m[1][0] * rhs[2] - rhs[1] * m[2][0])
            + rhs[0] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])) * invDet;

    return true;
}

MeshData dualContour(const float* field, int gridX, int gridY, int gridZ,
                     float isoLevel, float cellSize) {
    MeshData mesh;
    if (!field || gridX < 2 || gridY < 2 || gridZ < 2)
        return mesh;

    int cellsX = gridX - 1;
    int cellsY = gridY - 1;
    int cellsZ = gridZ - 1;

    std::vector<int> cellToVertex(cellsX * cellsY * cellsZ, -1);

    auto fieldIdx = [&](int x, int y, int z) -> int {
        return z * gridY * gridX + y * gridX + x;
    };

    auto cellIdx = [&](int x, int y, int z) -> int {
        return z * cellsY * cellsX + y * cellsX + x;
    };

    auto sampleField = [&](int x, int y, int z) -> float {
        return field[fieldIdx(x, y, z)];
    };

    // Gradient via central differences, clamped to grid boundaries
    auto gradient = [&](int x, int y, int z, float& gx, float& gy, float& gz) {
        int xm = (x > 0) ? x - 1 : x;
        int xp = (x < gridX - 1) ? x + 1 : x;
        int ym = (y > 0) ? y - 1 : y;
        int yp = (y < gridY - 1) ? y + 1 : y;
        int zm = (z > 0) ? z - 1 : z;
        int zp = (z < gridZ - 1) ? z + 1 : z;

        float dx = (float)(xp - xm);
        float dy = (float)(yp - ym);
        float dz = (float)(zp - zm);
        if (dx == 0) dx = 1.0f;
        if (dy == 0) dy = 1.0f;
        if (dz == 0) dz = 1.0f;

        gx = (sampleField(xp, y, z) - sampleField(xm, y, z)) / dx;
        gy = (sampleField(x, yp, z) - sampleField(x, ym, z)) / dy;
        gz = (sampleField(x, y, zp) - sampleField(x, y, zm)) / dz;
    };

    // Gradient at an arbitrary fractional position via interpolation of
    // the two endpoint gradients
    auto gradientLerp = [&](int x0, int y0, int z0, int x1, int y1, int z1,
                            float t, float& gx, float& gy, float& gz) {
        float g0x, g0y, g0z, g1x, g1y, g1z;
        gradient(x0, y0, z0, g0x, g0y, g0z);
        gradient(x1, y1, z1, g1x, g1y, g1z);
        gx = g0x + t * (g1x - g0x);
        gy = g0y + t * (g1y - g0y);
        gz = g0z + t * (g1z - g0z);
    };

    // 8 corners of a cell
    static const int cornerOffsets[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
    };

    // 12 edges of a cube
    static const int edges[12][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7}, // X-aligned
        {0, 2}, {1, 3}, {4, 6}, {5, 7}, // Y-aligned
        {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Z-aligned
    };

    const float lambda = 0.01f; // Tikhonov regularization

    // Phase 1: Place vertices using QEF minimization
    for (int cz = 0; cz < cellsZ; ++cz) {
        for (int cy = 0; cy < cellsY; ++cy) {
            for (int cx = 0; cx < cellsX; ++cx) {
                float cornerVals[8];
                for (int i = 0; i < 8; ++i) {
                    int gx = cx + cornerOffsets[i][0];
                    int gy = cy + cornerOffsets[i][1];
                    int gz = cz + cornerOffsets[i][2];
                    cornerVals[i] = field[fieldIdx(gx, gy, gz)];
                }

                bool hasPositive = false, hasNegative = false;
                for (int i = 0; i < 8; ++i) {
                    if (cornerVals[i] >= isoLevel) hasPositive = true;
                    else hasNegative = true;
                }
                if (!hasPositive || !hasNegative)
                    continue;

                // Collect edge crossing points and normals for QEF
                // Also compute the average for fallback
                float avgX = 0, avgY = 0, avgZ = 0;
                int crossingCount = 0;

                // QEF accumulators: A^T A (symmetric 3x3) and A^T b (3-vector)
                float ata[3][3] = {};
                float atb[3] = {};

                for (int e = 0; e < 12; ++e) {
                    int c0 = edges[e][0];
                    int c1 = edges[e][1];
                    float v0 = cornerVals[c0];
                    float v1 = cornerVals[c1];

                    bool s0 = (v0 >= isoLevel);
                    bool s1 = (v1 >= isoLevel);
                    if (s0 == s1) continue;

                    float t = (isoLevel - v0) / (v1 - v0);

                    float px = (float)cx + cornerOffsets[c0][0] + t * (cornerOffsets[c1][0] - cornerOffsets[c0][0]);
                    float py = (float)cy + cornerOffsets[c0][1] + t * (cornerOffsets[c1][1] - cornerOffsets[c0][1]);
                    float pz = (float)cz + cornerOffsets[c0][2] + t * (cornerOffsets[c1][2] - cornerOffsets[c0][2]);

                    // Gradient (normal) at crossing point
                    int gx0 = cx + cornerOffsets[c0][0];
                    int gy0 = cy + cornerOffsets[c0][1];
                    int gz0 = cz + cornerOffsets[c0][2];
                    int gx1 = cx + cornerOffsets[c1][0];
                    int gy1 = cy + cornerOffsets[c1][1];
                    int gz1 = cz + cornerOffsets[c1][2];

                    float nx, ny, nz;
                    gradientLerp(gx0, gy0, gz0, gx1, gy1, gz1, t, nx, ny, nz);

                    // Normalize the gradient
                    float len = std::sqrt(nx * nx + ny * ny + nz * nz);
                    if (len > 1e-8f) {
                        nx /= len;
                        ny /= len;
                        nz /= len;
                    }

                    // Accumulate A^T A (outer product of normal)
                    ata[0][0] += nx * nx;
                    ata[0][1] += nx * ny;
                    ata[0][2] += nx * nz;
                    ata[1][0] += ny * nx;
                    ata[1][1] += ny * ny;
                    ata[1][2] += ny * nz;
                    ata[2][0] += nz * nx;
                    ata[2][1] += nz * ny;
                    ata[2][2] += nz * nz;

                    // Accumulate A^T b where b_i = dot(n_i, p_i)
                    float d = nx * px + ny * py + nz * pz;
                    atb[0] += nx * d;
                    atb[1] += ny * d;
                    atb[2] += nz * d;

                    avgX += px;
                    avgY += py;
                    avgZ += pz;
                    crossingCount++;
                }

                if (crossingCount == 0) continue;

                float inv = 1.0f / crossingCount;
                avgX *= inv;
                avgY *= inv;
                avgZ *= inv;

                // Add Tikhonov regularization: (A^T A + lambda * I)
                ata[0][0] += lambda;
                ata[1][1] += lambda;
                ata[2][2] += lambda;

                // Regularization also shifts the solution toward the origin,
                // so add lambda * massPoint to atb to bias toward the average
                atb[0] += lambda * avgX;
                atb[1] += lambda * avgY;
                atb[2] += lambda * avgZ;

                float qefPos[3];
                bool solved = solve3x3(ata, atb, qefPos);

                float finalX, finalY, finalZ;
                if (solved) {
                    finalX = qefPos[0];
                    finalY = qefPos[1];
                    finalZ = qefPos[2];

                    // Clamp to cell bounding box
                    float minX = (float)cx;
                    float maxX = (float)(cx + 1);
                    float minY = (float)cy;
                    float maxY = (float)(cy + 1);
                    float minZ = (float)cz;
                    float maxZ = (float)(cz + 1);

                    if (finalX < minX) finalX = minX;
                    if (finalX > maxX) finalX = maxX;
                    if (finalY < minY) finalY = minY;
                    if (finalY > maxY) finalY = maxY;
                    if (finalZ < minZ) finalZ = minZ;
                    if (finalZ > maxZ) finalZ = maxZ;
                } else {
                    // Fallback to surface nets average
                    finalX = avgX;
                    finalY = avgY;
                    finalZ = avgZ;
                }

                int vidx = static_cast<int>(mesh.positions.size() / 3);
                cellToVertex[cellIdx(cx, cy, cz)] = vidx;
                mesh.positions.push_back(finalX * cellSize);
                mesh.positions.push_back(finalY * cellSize);
                mesh.positions.push_back(finalZ * cellSize);
            }
        }
    }

    if (mesh.positions.empty())
        return mesh;

    // Phase 2: Connect faces (identical to surface nets)
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
                        int cy0 = gy, cy1 = gy - 1;
                        int cz0 = gz, cz1 = gz - 1;
                        if (cy0 >= 0 && cy0 < cellsY && cy1 >= 0 && cy1 < cellsY &&
                            cz0 >= 0 && cz0 < cellsZ && cz1 >= 0 && cz1 < cellsZ &&
                            gx < cellsX) {
                            int a = cellToVertex[cellIdx(gx, cy0, cz0)];
                            int b = cellToVertex[cellIdx(gx, cy1, cz0)];
                            int c = cellToVertex[cellIdx(gx, cy1, cz1)];
                            int d = cellToVertex[cellIdx(gx, cy0, cz1)];
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
