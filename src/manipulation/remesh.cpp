#include "bromesh/manipulation/remesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bromesh {

// Isotropic remeshing via iterative edge operations:
// 1. Split edges longer than 4/3 * target
// 2. Collapse edges shorter than 4/5 * target
// 3. Flip edges to improve vertex valence (target=6 for interior)
// 4. Tangential relaxation

struct HalfEdgeMesh {
    struct Vertex {
        float pos[3];
        bool boundary = false;
    };
    std::vector<Vertex> verts;
    std::vector<std::array<uint32_t, 3>> tris; // vertex indices per triangle

    float edgeLength(uint32_t a, uint32_t b) const {
        float dx = verts[a].pos[0] - verts[b].pos[0];
        float dy = verts[a].pos[1] - verts[b].pos[1];
        float dz = verts[a].pos[2] - verts[b].pos[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

MeshData remeshIsotropic(const MeshData& mesh, float targetEdgeLength,
                         int iterations) {
    if (mesh.empty() || mesh.indices.empty() || iterations <= 0) return mesh;

    // Build working copy
    HalfEdgeMesh hm;
    size_t vertCount = mesh.vertexCount();
    hm.verts.resize(vertCount);
    for (size_t v = 0; v < vertCount; ++v) {
        hm.verts[v].pos[0] = mesh.positions[v * 3 + 0];
        hm.verts[v].pos[1] = mesh.positions[v * 3 + 1];
        hm.verts[v].pos[2] = mesh.positions[v * 3 + 2];
    }

    size_t triCount = mesh.triangleCount();
    hm.tris.resize(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        hm.tris[t] = {mesh.indices[t*3+0], mesh.indices[t*3+1], mesh.indices[t*3+2]};
    }

    // Compute target edge length if not specified
    if (targetEdgeLength <= 0.0f) {
        double totalLen = 0.0;
        int edgeCount = 0;
        for (const auto& tri : hm.tris) {
            for (int e = 0; e < 3; ++e) {
                totalLen += hm.edgeLength(tri[e], tri[(e+1)%3]);
                edgeCount++;
            }
        }
        targetEdgeLength = static_cast<float>(totalLen / edgeCount);
    }

    float splitThreshold = targetEdgeLength * (4.0f / 3.0f);
    float collapseThreshold = targetEdgeLength * (4.0f / 5.0f);

    struct EdgeKey {
        uint32_t a, b;
        EdgeKey(uint32_t x, uint32_t y)
            : a(std::min(x, y)), b(std::max(x, y)) {}
        bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
    };
    struct EdgeKeyHash {
        size_t operator()(const EdgeKey& k) const {
            return std::hash<uint64_t>()(
                (static_cast<uint64_t>(k.a) << 32) | k.b);
        }
    };

    for (int iter = 0; iter < iterations; ++iter) {
        // --- Step 1: Split long edges ---
        {
            size_t numTris = hm.tris.size();
            for (size_t t = 0; t < numTris; ++t) {
                // Find longest edge
                float maxLen = 0;
                int maxEdge = -1;
                for (int e = 0; e < 3; ++e) {
                    float len = hm.edgeLength(hm.tris[t][e], hm.tris[t][(e+1)%3]);
                    if (len > maxLen) { maxLen = len; maxEdge = e; }
                }
                if (maxLen > splitThreshold && maxEdge >= 0) {
                    uint32_t a = hm.tris[t][maxEdge];
                    uint32_t b = hm.tris[t][(maxEdge+1)%3];
                    uint32_t c = hm.tris[t][(maxEdge+2)%3];

                    // Create midpoint vertex
                    HalfEdgeMesh::Vertex mid;
                    for (int ch = 0; ch < 3; ++ch)
                        mid.pos[ch] = (hm.verts[a].pos[ch] + hm.verts[b].pos[ch]) * 0.5f;
                    uint32_t midIdx = static_cast<uint32_t>(hm.verts.size());
                    hm.verts.push_back(mid);

                    // Replace triangle t with (a, mid, c), add (mid, b, c)
                    hm.tris[t] = {a, midIdx, c};
                    hm.tris.push_back({midIdx, b, c});
                }
            }
        }

        // --- Step 2: Collapse short edges ---
        {
            // Build edge-to-face map
            std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash> edgeFaces;
            for (size_t t = 0; t < hm.tris.size(); ++t) {
                for (int e = 0; e < 3; ++e) {
                    EdgeKey key(hm.tris[t][e], hm.tris[t][(e+1)%3]);
                    edgeFaces[key].push_back(t);
                }
            }

            std::unordered_set<size_t> removedTris;
            std::vector<bool> removedVerts(hm.verts.size(), false);

            for (auto& [key, faces] : edgeFaces) {
                if (removedVerts[key.a] || removedVerts[key.b]) continue;
                float len = hm.edgeLength(key.a, key.b);
                if (len < collapseThreshold) {
                    // Collapse b into a (move a to midpoint)
                    for (int ch = 0; ch < 3; ++ch)
                        hm.verts[key.a].pos[ch] =
                            (hm.verts[key.a].pos[ch] + hm.verts[key.b].pos[ch]) * 0.5f;
                    removedVerts[key.b] = true;

                    // Replace all references to b with a
                    for (auto& tri : hm.tris) {
                        for (int i = 0; i < 3; ++i) {
                            if (tri[i] == key.b) tri[i] = key.a;
                        }
                    }

                    // Mark degenerate triangles for removal
                    for (size_t fi : faces) {
                        removedTris.insert(fi);
                    }
                }
            }

            // Remove degenerate triangles (two or more same indices)
            std::vector<std::array<uint32_t, 3>> newTris;
            newTris.reserve(hm.tris.size());
            for (size_t t = 0; t < hm.tris.size(); ++t) {
                if (removedTris.count(t)) continue;
                const auto& tri = hm.tris[t];
                if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) continue;
                newTris.push_back(tri);
            }
            hm.tris = std::move(newTris);
        }

        // --- Step 3: Edge flip to improve valence ---
        {
            // Build vertex valence
            std::vector<int> valence(hm.verts.size(), 0);
            for (const auto& tri : hm.tris) {
                valence[tri[0]]++;
                valence[tri[1]]++;
                valence[tri[2]]++;
            }

            // Build edge-to-face for flip candidates
            std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash> edgeFaces;
            for (size_t t = 0; t < hm.tris.size(); ++t) {
                for (int e = 0; e < 3; ++e) {
                    EdgeKey key(hm.tris[t][e], hm.tris[t][(e+1)%3]);
                    edgeFaces[key].push_back(t);
                }
            }

            for (auto& [key, faces] : edgeFaces) {
                if (faces.size() != 2) continue;

                uint32_t a = key.a, b = key.b;
                // Find opposite vertices
                uint32_t c = UINT32_MAX, d = UINT32_MAX;
                for (int i = 0; i < 3; ++i) {
                    uint32_t v = hm.tris[faces[0]][i];
                    if (v != a && v != b) c = v;
                    v = hm.tris[faces[1]][i];
                    if (v != a && v != b) d = v;
                }
                if (c == UINT32_MAX || d == UINT32_MAX || c == d) continue;

                // Check if flip improves valence deviation from 6
                int devBefore = std::abs(valence[a] - 6) + std::abs(valence[b] - 6)
                              + std::abs(valence[c] - 6) + std::abs(valence[d] - 6);
                int devAfter = std::abs(valence[a] - 1 - 6) + std::abs(valence[b] - 1 - 6)
                             + std::abs(valence[c] + 1 - 6) + std::abs(valence[d] + 1 - 6);

                if (devAfter < devBefore) {
                    // Flip: replace (a,b,c) + (a,d,b) with (c,d,a) + (c,b,d)
                    hm.tris[faces[0]] = {c, d, a};
                    hm.tris[faces[1]] = {c, b, d};
                    valence[a]--;
                    valence[b]--;
                    valence[c]++;
                    valence[d]++;
                }
            }
        }

        // --- Step 4: Tangential relaxation ---
        {
            // Build adjacency
            std::vector<std::vector<uint32_t>> adj(hm.verts.size());
            for (const auto& tri : hm.tris) {
                for (int i = 0; i < 3; ++i) {
                    uint32_t a = tri[i], b = tri[(i+1)%3];
                    auto addUnique = [](std::vector<uint32_t>& list, uint32_t v) {
                        for (uint32_t x : list) if (x == v) return;
                        list.push_back(v);
                    };
                    addUnique(adj[a], b);
                    addUnique(adj[b], a);
                }
            }

            // Compute vertex normals
            std::vector<float> normals(hm.verts.size() * 3, 0.0f);
            for (const auto& tri : hm.tris) {
                float e1[3], e2[3];
                for (int c = 0; c < 3; ++c) {
                    e1[c] = hm.verts[tri[1]].pos[c] - hm.verts[tri[0]].pos[c];
                    e2[c] = hm.verts[tri[2]].pos[c] - hm.verts[tri[0]].pos[c];
                }
                float nx = e1[1]*e2[2] - e1[2]*e2[1];
                float ny = e1[2]*e2[0] - e1[0]*e2[2];
                float nz = e1[0]*e2[1] - e1[1]*e2[0];
                for (int i = 0; i < 3; ++i) {
                    normals[tri[i]*3+0] += nx;
                    normals[tri[i]*3+1] += ny;
                    normals[tri[i]*3+2] += nz;
                }
            }
            for (size_t v = 0; v < hm.verts.size(); ++v) {
                float* n = &normals[v*3];
                float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
                if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
            }

            // Relax: move toward neighbor centroid, then project back to tangent plane
            for (size_t v = 0; v < hm.verts.size(); ++v) {
                if (adj[v].empty()) continue;

                float cx = 0, cy = 0, cz = 0;
                for (uint32_t n : adj[v]) {
                    cx += hm.verts[n].pos[0];
                    cy += hm.verts[n].pos[1];
                    cz += hm.verts[n].pos[2];
                }
                float invN = 1.0f / adj[v].size();
                cx *= invN; cy *= invN; cz *= invN;

                // Displacement toward centroid
                float dx = cx - hm.verts[v].pos[0];
                float dy = cy - hm.verts[v].pos[1];
                float dz = cz - hm.verts[v].pos[2];

                // Project onto tangent plane
                float* n = &normals[v * 3];
                float dot = dx*n[0] + dy*n[1] + dz*n[2];
                dx -= dot * n[0];
                dy -= dot * n[1];
                dz -= dot * n[2];

                hm.verts[v].pos[0] += dx * 0.5f;
                hm.verts[v].pos[1] += dy * 0.5f;
                hm.verts[v].pos[2] += dz * 0.5f;
            }
        }
    }

    // Compact: remove unused vertices and build output MeshData
    std::vector<uint32_t> remap(hm.verts.size(), UINT32_MAX);
    MeshData result;

    for (const auto& tri : hm.tris) {
        for (int i = 0; i < 3; ++i) {
            if (remap[tri[i]] == UINT32_MAX) {
                remap[tri[i]] = static_cast<uint32_t>(result.positions.size() / 3);
                result.positions.push_back(hm.verts[tri[i]].pos[0]);
                result.positions.push_back(hm.verts[tri[i]].pos[1]);
                result.positions.push_back(hm.verts[tri[i]].pos[2]);
            }
        }
        result.indices.push_back(remap[tri[0]]);
        result.indices.push_back(remap[tri[1]]);
        result.indices.push_back(remap[tri[2]]);
    }

    // Compute normals
    size_t outVertCount = result.vertexCount();
    result.normals.assign(outVertCount * 3, 0.0f);
    for (size_t t = 0; t < result.triangleCount(); ++t) {
        uint32_t i0 = result.indices[t*3+0];
        uint32_t i1 = result.indices[t*3+1];
        uint32_t i2 = result.indices[t*3+2];
        float e1x = result.positions[i1*3+0]-result.positions[i0*3+0];
        float e1y = result.positions[i1*3+1]-result.positions[i0*3+1];
        float e1z = result.positions[i1*3+2]-result.positions[i0*3+2];
        float e2x = result.positions[i2*3+0]-result.positions[i0*3+0];
        float e2y = result.positions[i2*3+1]-result.positions[i0*3+1];
        float e2z = result.positions[i2*3+2]-result.positions[i0*3+2];
        float nx = e1y*e2z - e1z*e2y;
        float ny = e1z*e2x - e1x*e2z;
        float nz = e1x*e2y - e1y*e2x;
        for (uint32_t idx : {i0, i1, i2}) {
            result.normals[idx*3+0] += nx;
            result.normals[idx*3+1] += ny;
            result.normals[idx*3+2] += nz;
        }
    }
    for (size_t v = 0; v < outVertCount; ++v) {
        float* n = &result.normals[v*3];
        float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
    }

    return result;
}

} // namespace bromesh
