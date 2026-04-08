#include "bromesh/manipulation/split_components.h"

#include <queue>
#include <unordered_map>
#include <vector>

namespace bromesh {

std::vector<MeshData> splitConnectedComponents(const MeshData& mesh) {
    if (mesh.empty()) return {};

    const size_t triCount = mesh.triangleCount();
    const size_t vCount = mesh.vertexCount();
    const bool hasN = mesh.hasNormals();
    const bool hasUV = mesh.hasUVs();
    const bool hasC = mesh.hasColors();

    // Build adjacency: for each vertex, which triangles reference it
    std::vector<std::vector<uint32_t>> vertToTris(vCount);
    for (size_t t = 0; t < triCount; ++t) {
        vertToTris[mesh.indices[t * 3 + 0]].push_back(static_cast<uint32_t>(t));
        vertToTris[mesh.indices[t * 3 + 1]].push_back(static_cast<uint32_t>(t));
        vertToTris[mesh.indices[t * 3 + 2]].push_back(static_cast<uint32_t>(t));
    }

    // BFS to find connected components of triangles
    std::vector<bool> visited(triCount, false);
    std::vector<std::vector<uint32_t>> components; // each is a list of triangle indices

    for (size_t t = 0; t < triCount; ++t) {
        if (visited[t]) continue;

        std::vector<uint32_t> comp;
        std::queue<uint32_t> q;
        q.push(static_cast<uint32_t>(t));
        visited[t] = true;

        while (!q.empty()) {
            uint32_t tri = q.front();
            q.pop();
            comp.push_back(tri);

            // For each vertex of this triangle, visit all adjacent triangles
            for (int vi = 0; vi < 3; ++vi) {
                uint32_t v = mesh.indices[tri * 3 + vi];
                for (uint32_t adjTri : vertToTris[v]) {
                    if (!visited[adjTri]) {
                        visited[adjTri] = true;
                        q.push(adjTri);
                    }
                }
            }
        }

        components.push_back(std::move(comp));
    }

    // Build a MeshData for each component with remapped vertex indices
    std::vector<MeshData> result;
    result.reserve(components.size());

    for (auto& comp : components) {
        MeshData out;
        // Map old vertex index -> new vertex index in this component
        std::unordered_map<uint32_t, uint32_t> vertMap;
        uint32_t newIdx = 0;

        out.indices.reserve(comp.size() * 3);

        for (uint32_t tri : comp) {
            for (int vi = 0; vi < 3; ++vi) {
                uint32_t oldV = mesh.indices[tri * 3 + vi];
                auto it = vertMap.find(oldV);
                uint32_t nv;
                if (it == vertMap.end()) {
                    nv = newIdx++;
                    vertMap[oldV] = nv;

                    out.positions.push_back(mesh.positions[oldV * 3 + 0]);
                    out.positions.push_back(mesh.positions[oldV * 3 + 1]);
                    out.positions.push_back(mesh.positions[oldV * 3 + 2]);
                    if (hasN) {
                        out.normals.push_back(mesh.normals[oldV * 3 + 0]);
                        out.normals.push_back(mesh.normals[oldV * 3 + 1]);
                        out.normals.push_back(mesh.normals[oldV * 3 + 2]);
                    }
                    if (hasUV) {
                        out.uvs.push_back(mesh.uvs[oldV * 2 + 0]);
                        out.uvs.push_back(mesh.uvs[oldV * 2 + 1]);
                    }
                    if (hasC) {
                        out.colors.push_back(mesh.colors[oldV * 4 + 0]);
                        out.colors.push_back(mesh.colors[oldV * 4 + 1]);
                        out.colors.push_back(mesh.colors[oldV * 4 + 2]);
                        out.colors.push_back(mesh.colors[oldV * 4 + 3]);
                    }
                } else {
                    nv = it->second;
                }
                out.indices.push_back(nv);
            }
        }

        result.push_back(std::move(out));
    }

    return result;
}

} // namespace bromesh
