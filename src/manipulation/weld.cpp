#include "bromesh/manipulation/weld.h"

#include <cmath>
#include <unordered_map>
#include <vector>

namespace bromesh {

// Spatial hash key from quantized grid cell coordinates
struct CellKey {
    int x, y, z;
    bool operator==(const CellKey& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

struct CellKeyHash {
    size_t operator()(const CellKey& k) const {
        // Large primes for spatial hashing
        size_t h = static_cast<size_t>(k.x) * 73856093u;
        h ^= static_cast<size_t>(k.y) * 19349663u;
        h ^= static_cast<size_t>(k.z) * 83492791u;
        return h;
    }
};

MeshData weldVertices(const MeshData& mesh, float epsilon) {
    if (mesh.empty()) return {};

    const size_t vCount = mesh.vertexCount();
    const bool hasN = mesh.hasNormals();
    const bool hasUV = mesh.hasUVs();
    const bool hasC = mesh.hasColors();

    float cellSize = epsilon > 0.0f ? epsilon : 1e-5f;
    float invCell = 1.0f / cellSize;

    // Map from grid cell to list of vertex indices that landed in that cell
    std::unordered_map<CellKey, std::vector<uint32_t>, CellKeyHash> grid;

    // oldToNew[i] = new index for old vertex i
    std::vector<uint32_t> oldToNew(vCount);

    MeshData result;
    // We'll build the new vertex buffer incrementally
    uint32_t newCount = 0;

    for (size_t i = 0; i < vCount; ++i) {
        float px = mesh.positions[i * 3 + 0];
        float py = mesh.positions[i * 3 + 1];
        float pz = mesh.positions[i * 3 + 2];

        int cx = static_cast<int>(std::floor(px * invCell));
        int cy = static_cast<int>(std::floor(py * invCell));
        int cz = static_cast<int>(std::floor(pz * invCell));

        // Check this cell and 26 neighbors for a match
        uint32_t merged = UINT32_MAX;
        for (int dx = -1; dx <= 1 && merged == UINT32_MAX; ++dx) {
            for (int dy = -1; dy <= 1 && merged == UINT32_MAX; ++dy) {
                for (int dz = -1; dz <= 1 && merged == UINT32_MAX; ++dz) {
                    CellKey nk{cx + dx, cy + dy, cz + dz};
                    auto it = grid.find(nk);
                    if (it == grid.end()) continue;
                    for (uint32_t ni : it->second) {
                        float ex = mesh.positions[ni * 3 + 0] - px;
                        float ey = mesh.positions[ni * 3 + 1] - py;
                        float ez = mesh.positions[ni * 3 + 2] - pz;
                        if (ex * ex + ey * ey + ez * ez <= epsilon * epsilon) {
                            merged = oldToNew[ni];
                            break;
                        }
                    }
                }
            }
        }

        if (merged != UINT32_MAX) {
            oldToNew[i] = merged;
        } else {
            // Add new vertex
            oldToNew[i] = newCount;
            result.positions.push_back(px);
            result.positions.push_back(py);
            result.positions.push_back(pz);
            if (hasN) {
                result.normals.push_back(mesh.normals[i * 3 + 0]);
                result.normals.push_back(mesh.normals[i * 3 + 1]);
                result.normals.push_back(mesh.normals[i * 3 + 2]);
            }
            if (hasUV) {
                result.uvs.push_back(mesh.uvs[i * 2 + 0]);
                result.uvs.push_back(mesh.uvs[i * 2 + 1]);
            }
            if (hasC) {
                result.colors.push_back(mesh.colors[i * 4 + 0]);
                result.colors.push_back(mesh.colors[i * 4 + 1]);
                result.colors.push_back(mesh.colors[i * 4 + 2]);
                result.colors.push_back(mesh.colors[i * 4 + 3]);
            }
            newCount++;
        }

        CellKey ck{cx, cy, cz};
        grid[ck].push_back(static_cast<uint32_t>(i));
    }

    // Remap indices, skip degenerate triangles
    result.indices.reserve(mesh.indices.size());
    for (size_t t = 0; t < mesh.indices.size(); t += 3) {
        uint32_t a = oldToNew[mesh.indices[t + 0]];
        uint32_t b = oldToNew[mesh.indices[t + 1]];
        uint32_t c = oldToNew[mesh.indices[t + 2]];
        if (a != b && b != c && a != c) {
            result.indices.push_back(a);
            result.indices.push_back(b);
            result.indices.push_back(c);
        }
    }

    return result;
}

} // namespace bromesh
