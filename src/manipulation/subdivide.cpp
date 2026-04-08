#include "bromesh/manipulation/subdivide.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace bromesh {

// Edge key: ordered pair of vertex indices
struct EdgeKey {
    uint32_t a, b;
    EdgeKey(uint32_t x, uint32_t y) : a(std::min(x, y)), b(std::max(x, y)) {}
    bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const {
        return std::hash<uint64_t>()(
            (static_cast<uint64_t>(k.a) << 32) | k.b);
    }
};

// ----- Midpoint (no smoothing) -----

static MeshData subdivideMidpointOnce(const MeshData& mesh) {
    if (mesh.empty() || mesh.indices.empty()) return mesh;

    size_t vertCount = mesh.vertexCount();
    size_t triCount = mesh.triangleCount();
    bool hasN = mesh.hasNormals();
    bool hasUV = mesh.hasUVs();
    bool hasC = mesh.hasColors();

    // Edge midpoint map: edge -> new vertex index
    std::unordered_map<EdgeKey, uint32_t, EdgeKeyHash> edgeMid;
    edgeMid.reserve(triCount * 3);

    MeshData out;
    // Start with copies of original vertices
    out.positions = mesh.positions;
    if (hasN) out.normals = mesh.normals;
    if (hasUV) out.uvs = mesh.uvs;
    if (hasC) out.colors = mesh.colors;
    out.indices.reserve(triCount * 4 * 3);

    auto getOrCreateMid = [&](uint32_t i0, uint32_t i1) -> uint32_t {
        EdgeKey key(i0, i1);
        auto it = edgeMid.find(key);
        if (it != edgeMid.end()) return it->second;

        uint32_t newIdx = static_cast<uint32_t>(out.positions.size() / 3);
        edgeMid[key] = newIdx;

        // Midpoint position
        for (int c = 0; c < 3; ++c)
            out.positions.push_back(
                (mesh.positions[i0 * 3 + c] + mesh.positions[i1 * 3 + c]) * 0.5f);
        if (hasN) {
            for (int c = 0; c < 3; ++c)
                out.normals.push_back(
                    (mesh.normals[i0 * 3 + c] + mesh.normals[i1 * 3 + c]) * 0.5f);
            // Renormalize
            float* n = &out.normals[newIdx * 3];
            float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            if (len > 1e-8f) { n[0] /= len; n[1] /= len; n[2] /= len; }
        }
        if (hasUV) {
            for (int c = 0; c < 2; ++c)
                out.uvs.push_back(
                    (mesh.uvs[i0 * 2 + c] + mesh.uvs[i1 * 2 + c]) * 0.5f);
        }
        if (hasC) {
            for (int c = 0; c < 4; ++c)
                out.colors.push_back(
                    (mesh.colors[i0 * 4 + c] + mesh.colors[i1 * 4 + c]) * 0.5f);
        }
        return newIdx;
    };

    // Split each triangle into 4
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v0 = mesh.indices[t * 3 + 0];
        uint32_t v1 = mesh.indices[t * 3 + 1];
        uint32_t v2 = mesh.indices[t * 3 + 2];

        uint32_t m01 = getOrCreateMid(v0, v1);
        uint32_t m12 = getOrCreateMid(v1, v2);
        uint32_t m20 = getOrCreateMid(v2, v0);

        // 4 sub-triangles
        out.indices.push_back(v0);  out.indices.push_back(m01); out.indices.push_back(m20);
        out.indices.push_back(m01); out.indices.push_back(v1);  out.indices.push_back(m12);
        out.indices.push_back(m20); out.indices.push_back(m12); out.indices.push_back(v2);
        out.indices.push_back(m01); out.indices.push_back(m12); out.indices.push_back(m20);
    }

    return out;
}

MeshData subdivideMidpoint(const MeshData& mesh, int iterations) {
    if (iterations <= 0) return mesh;
    MeshData result = mesh;
    for (int i = 0; i < iterations; ++i)
        result = subdivideMidpointOnce(result);
    return result;
}

// ----- Loop subdivision -----

static MeshData subdivideLoopOnce(const MeshData& mesh) {
    if (mesh.empty() || mesh.indices.empty()) return mesh;

    size_t vertCount = mesh.vertexCount();
    size_t triCount = mesh.triangleCount();
    bool hasUV = mesh.hasUVs();
    bool hasC = mesh.hasColors();

    // Build adjacency: for each edge, store the two opposite vertices
    struct EdgeAdj {
        uint32_t midIdx = UINT32_MAX;
        uint32_t opposite[2] = {UINT32_MAX, UINT32_MAX};
        int faceCount = 0;
    };
    std::unordered_map<EdgeKey, EdgeAdj, EdgeKeyHash> edges;
    edges.reserve(triCount * 3);

    // Also build vertex valence and neighbor list
    std::vector<std::vector<uint32_t>> vertNeighbors(vertCount);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {
            mesh.indices[t * 3 + 0],
            mesh.indices[t * 3 + 1],
            mesh.indices[t * 3 + 2]
        };
        for (int e = 0; e < 3; ++e) {
            uint32_t a = v[e], b = v[(e + 1) % 3], c = v[(e + 2) % 3];
            EdgeKey key(a, b);
            auto& adj = edges[key];
            if (adj.faceCount < 2) adj.opposite[adj.faceCount] = c;
            adj.faceCount++;

            // Track neighbors
            auto addNeighbor = [](std::vector<uint32_t>& list, uint32_t n) {
                for (uint32_t x : list) if (x == n) return;
                list.push_back(n);
            };
            addNeighbor(vertNeighbors[a], b);
            addNeighbor(vertNeighbors[b], a);
        }
    }

    // Allocate output: start with space for new vertex positions
    MeshData out;
    // We'll build positions in two phases: moved originals + edge points

    // Phase 1: Compute new positions for original vertices (Loop smoothing)
    std::vector<float> newPos(vertCount * 3);
    for (size_t v = 0; v < vertCount; ++v) {
        int n = static_cast<int>(vertNeighbors[v].size());
        if (n < 2) {
            // Boundary or isolated vertex: keep as-is
            for (int c = 0; c < 3; ++c)
                newPos[v * 3 + c] = mesh.positions[v * 3 + c];
            continue;
        }

        // Check if boundary vertex (has any boundary edge)
        bool boundary = false;
        for (uint32_t nb : vertNeighbors[v]) {
            EdgeKey key(static_cast<uint32_t>(v), nb);
            auto it = edges.find(key);
            if (it != edges.end() && it->second.faceCount == 1) {
                boundary = true;
                break;
            }
        }

        if (boundary) {
            // Boundary rule: 3/4 * v + 1/8 * (b0 + b1) where b0,b1 are boundary neighbors
            std::vector<uint32_t> bndNeighbors;
            for (uint32_t nb : vertNeighbors[v]) {
                EdgeKey key(static_cast<uint32_t>(v), nb);
                auto it = edges.find(key);
                if (it != edges.end() && it->second.faceCount == 1)
                    bndNeighbors.push_back(nb);
            }
            if (bndNeighbors.size() == 2) {
                for (int c = 0; c < 3; ++c)
                    newPos[v * 3 + c] = mesh.positions[v * 3 + c] * 0.75f
                        + mesh.positions[bndNeighbors[0] * 3 + c] * 0.125f
                        + mesh.positions[bndNeighbors[1] * 3 + c] * 0.125f;
            } else {
                for (int c = 0; c < 3; ++c)
                    newPos[v * 3 + c] = mesh.positions[v * 3 + c];
            }
        } else {
            // Interior Loop rule
            float beta;
            if (n == 3)
                beta = 3.0f / 16.0f;
            else
                beta = 3.0f / (8.0f * n);

            float selfWeight = 1.0f - n * beta;
            for (int c = 0; c < 3; ++c) {
                float sum = 0.0f;
                for (uint32_t nb : vertNeighbors[v])
                    sum += mesh.positions[nb * 3 + c];
                newPos[v * 3 + c] = mesh.positions[v * 3 + c] * selfWeight + sum * beta;
            }
        }
    }

    // Build output vertices starting with smoothed originals
    out.positions.assign(newPos.begin(), newPos.end());
    if (hasUV) out.uvs.assign(mesh.uvs.begin(), mesh.uvs.end());
    if (hasC) out.colors.assign(mesh.colors.begin(), mesh.colors.end());

    // Phase 2: Create edge vertices
    for (auto& [key, adj] : edges) {
        uint32_t a = key.a, b = key.b;
        adj.midIdx = static_cast<uint32_t>(out.positions.size() / 3);

        if (adj.faceCount == 2 && adj.opposite[0] != UINT32_MAX && adj.opposite[1] != UINT32_MAX) {
            // Interior edge: 3/8 * (a+b) + 1/8 * (c+d)
            uint32_t c = adj.opposite[0], d = adj.opposite[1];
            for (int ch = 0; ch < 3; ++ch)
                out.positions.push_back(
                    (mesh.positions[a * 3 + ch] + mesh.positions[b * 3 + ch]) * 0.375f +
                    (mesh.positions[c * 3 + ch] + mesh.positions[d * 3 + ch]) * 0.125f);
        } else {
            // Boundary edge: midpoint
            for (int ch = 0; ch < 3; ++ch)
                out.positions.push_back(
                    (mesh.positions[a * 3 + ch] + mesh.positions[b * 3 + ch]) * 0.5f);
        }

        // UVs: always linear interpolation (to preserve seams)
        if (hasUV) {
            for (int ch = 0; ch < 2; ++ch)
                out.uvs.push_back(
                    (mesh.uvs[a * 2 + ch] + mesh.uvs[b * 2 + ch]) * 0.5f);
        }
        if (hasC) {
            for (int ch = 0; ch < 4; ++ch)
                out.colors.push_back(
                    (mesh.colors[a * 4 + ch] + mesh.colors[b * 4 + ch]) * 0.5f);
        }
    }

    // Phase 3: Build new triangles
    out.indices.reserve(triCount * 4 * 3);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v0 = mesh.indices[t * 3 + 0];
        uint32_t v1 = mesh.indices[t * 3 + 1];
        uint32_t v2 = mesh.indices[t * 3 + 2];

        uint32_t m01 = edges[EdgeKey(v0, v1)].midIdx;
        uint32_t m12 = edges[EdgeKey(v1, v2)].midIdx;
        uint32_t m20 = edges[EdgeKey(v2, v0)].midIdx;

        out.indices.push_back(v0);  out.indices.push_back(m01); out.indices.push_back(m20);
        out.indices.push_back(m01); out.indices.push_back(v1);  out.indices.push_back(m12);
        out.indices.push_back(m20); out.indices.push_back(m12); out.indices.push_back(v2);
        out.indices.push_back(m01); out.indices.push_back(m12); out.indices.push_back(m20);
    }

    // Compute normals for the new mesh
    out.normals.assign(out.positions.size(), 0.0f);
    size_t newTriCount = out.indices.size() / 3;
    for (size_t t = 0; t < newTriCount; ++t) {
        uint32_t i0 = out.indices[t * 3 + 0];
        uint32_t i1 = out.indices[t * 3 + 1];
        uint32_t i2 = out.indices[t * 3 + 2];
        float e1x = out.positions[i1*3+0] - out.positions[i0*3+0];
        float e1y = out.positions[i1*3+1] - out.positions[i0*3+1];
        float e1z = out.positions[i1*3+2] - out.positions[i0*3+2];
        float e2x = out.positions[i2*3+0] - out.positions[i0*3+0];
        float e2y = out.positions[i2*3+1] - out.positions[i0*3+1];
        float e2z = out.positions[i2*3+2] - out.positions[i0*3+2];
        float nx = e1y*e2z - e1z*e2y;
        float ny = e1z*e2x - e1x*e2z;
        float nz = e1x*e2y - e1y*e2x;
        for (uint32_t idx : {i0, i1, i2}) {
            out.normals[idx*3+0] += nx;
            out.normals[idx*3+1] += ny;
            out.normals[idx*3+2] += nz;
        }
    }
    size_t outVertCount = out.positions.size() / 3;
    for (size_t v = 0; v < outVertCount; ++v) {
        float* n = &out.normals[v * 3];
        float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (len > 1e-8f) { n[0] /= len; n[1] /= len; n[2] /= len; }
    }

    return out;
}

MeshData subdivideLoop(const MeshData& mesh, int iterations) {
    if (iterations <= 0) return mesh;
    MeshData result = mesh;
    for (int i = 0; i < iterations; ++i)
        result = subdivideLoopOnce(result);
    return result;
}

// ----- Catmull-Clark subdivision -----

static MeshData subdivideCCOnce(const MeshData& mesh) {
    if (mesh.empty() || mesh.indices.empty()) return mesh;

    size_t vertCount = mesh.vertexCount();
    size_t triCount = mesh.triangleCount();
    bool hasUV = mesh.hasUVs();
    bool hasC = mesh.hasColors();

    // Build adjacency
    struct EdgeAdj {
        uint32_t midIdx = UINT32_MAX;
        std::vector<size_t> faces; // face indices using this edge
    };
    std::unordered_map<EdgeKey, EdgeAdj, EdgeKeyHash> edges;
    edges.reserve(triCount * 3);

    // Faces adjacent to each vertex
    std::vector<std::vector<size_t>> vertFaces(vertCount);
    // Edges adjacent to each vertex
    std::vector<std::vector<EdgeKey>> vertEdges(vertCount);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {
            mesh.indices[t * 3 + 0],
            mesh.indices[t * 3 + 1],
            mesh.indices[t * 3 + 2]
        };
        for (int i = 0; i < 3; ++i) {
            vertFaces[v[i]].push_back(t);
            EdgeKey key(v[i], v[(i + 1) % 3]);
            edges[key].faces.push_back(t);

            auto addEdge = [](std::vector<EdgeKey>& list, const EdgeKey& e) {
                for (const auto& x : list) if (x == e) return;
                list.push_back(e);
            };
            addEdge(vertEdges[v[i]], key);
            addEdge(vertEdges[v[(i + 1) % 3]], key);
        }
    }

    // Phase 1: Face points (centroid of each face)
    std::vector<float> facePoints(triCount * 3);
    std::vector<float> faceUVs, faceColors;
    if (hasUV) faceUVs.resize(triCount * 2);
    if (hasC) faceColors.resize(triCount * 4);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];
        for (int c = 0; c < 3; ++c)
            facePoints[t * 3 + c] = (mesh.positions[i0 * 3 + c] +
                                      mesh.positions[i1 * 3 + c] +
                                      mesh.positions[i2 * 3 + c]) / 3.0f;
        if (hasUV) {
            for (int c = 0; c < 2; ++c)
                faceUVs[t * 2 + c] = (mesh.uvs[i0 * 2 + c] +
                                       mesh.uvs[i1 * 2 + c] +
                                       mesh.uvs[i2 * 2 + c]) / 3.0f;
        }
        if (hasC) {
            for (int c = 0; c < 4; ++c)
                faceColors[t * 4 + c] = (mesh.colors[i0 * 4 + c] +
                                          mesh.colors[i1 * 4 + c] +
                                          mesh.colors[i2 * 4 + c]) / 3.0f;
        }
    }

    // Phase 2: Edge points
    // Interior edge: avg of edge midpoint + avg of adjacent face points
    // Boundary edge: midpoint
    struct EdgePoint {
        float pos[3];
        float uv[2];
        float col[4];
    };
    std::unordered_map<EdgeKey, EdgePoint, EdgeKeyHash> edgePoints;
    for (auto& [key, adj] : edges) {
        EdgePoint ep{};
        if (adj.faces.size() == 2) {
            for (int c = 0; c < 3; ++c)
                ep.pos[c] = (mesh.positions[key.a * 3 + c] +
                             mesh.positions[key.b * 3 + c] +
                             facePoints[adj.faces[0] * 3 + c] +
                             facePoints[adj.faces[1] * 3 + c]) * 0.25f;
        } else {
            for (int c = 0; c < 3; ++c)
                ep.pos[c] = (mesh.positions[key.a * 3 + c] +
                             mesh.positions[key.b * 3 + c]) * 0.5f;
        }
        if (hasUV) {
            for (int c = 0; c < 2; ++c)
                ep.uv[c] = (mesh.uvs[key.a * 2 + c] + mesh.uvs[key.b * 2 + c]) * 0.5f;
        }
        if (hasC) {
            for (int c = 0; c < 4; ++c)
                ep.col[c] = (mesh.colors[key.a * 4 + c] + mesh.colors[key.b * 4 + c]) * 0.5f;
        }
        edgePoints[key] = ep;
    }

    // Phase 3: Move original vertices
    // Interior: V' = (F + 2R + (n-3)V) / n
    //   F = avg of face points, R = avg of edge midpoints, n = valence
    // Boundary: V' = (1/8)*(B0 + B1) + (3/4)*V
    std::vector<float> movedPos(vertCount * 3);
    for (size_t v = 0; v < vertCount; ++v) {
        int n = static_cast<int>(vertEdges[v].size());
        if (n < 2) {
            for (int c = 0; c < 3; ++c)
                movedPos[v * 3 + c] = mesh.positions[v * 3 + c];
            continue;
        }

        // Check boundary
        bool boundary = false;
        for (const auto& ek : vertEdges[v]) {
            if (edges[ek].faces.size() == 1) { boundary = true; break; }
        }

        if (boundary) {
            std::vector<uint32_t> bndNeighbors;
            for (const auto& ek : vertEdges[v]) {
                if (edges[ek].faces.size() == 1) {
                    uint32_t other = (ek.a == v) ? ek.b : ek.a;
                    bndNeighbors.push_back(other);
                }
            }
            if (bndNeighbors.size() == 2) {
                for (int c = 0; c < 3; ++c)
                    movedPos[v * 3 + c] = mesh.positions[v * 3 + c] * 0.75f
                        + mesh.positions[bndNeighbors[0] * 3 + c] * 0.125f
                        + mesh.positions[bndNeighbors[1] * 3 + c] * 0.125f;
            } else {
                for (int c = 0; c < 3; ++c)
                    movedPos[v * 3 + c] = mesh.positions[v * 3 + c];
            }
        } else {
            // F: average of adjacent face points
            float F[3] = {0, 0, 0};
            for (size_t fi : vertFaces[v])
                for (int c = 0; c < 3; ++c)
                    F[c] += facePoints[fi * 3 + c];
            int nf = static_cast<int>(vertFaces[v].size());
            for (int c = 0; c < 3; ++c) F[c] /= nf;

            // R: average of edge midpoints
            float R[3] = {0, 0, 0};
            for (const auto& ek : vertEdges[v])
                for (int c = 0; c < 3; ++c)
                    R[c] += (mesh.positions[ek.a * 3 + c] +
                             mesh.positions[ek.b * 3 + c]) * 0.5f;
            for (int c = 0; c < 3; ++c) R[c] /= n;

            // V' = (F + 2R + (n-3)V) / n
            for (int c = 0; c < 3; ++c)
                movedPos[v * 3 + c] = (F[c] + 2.0f * R[c] +
                    (n - 3.0f) * mesh.positions[v * 3 + c]) / n;
        }
    }

    // Phase 4: Build output mesh
    // For each triangle, create 3 quads (each split into 2 triangles):
    //   face_point -- edge_point_01 -- v0 -- edge_point_20
    //   face_point -- edge_point_12 -- v1 -- edge_point_01
    //   face_point -- edge_point_20 -- v2 -- edge_point_12

    MeshData out;
    // Add moved original vertices
    out.positions.assign(movedPos.begin(), movedPos.end());
    if (hasUV) out.uvs = mesh.uvs; // keep original UVs for original verts
    if (hasC) out.colors = mesh.colors;

    // Add edge point vertices
    std::unordered_map<EdgeKey, uint32_t, EdgeKeyHash> edgeVertIdx;
    for (auto& [key, ep] : edgePoints) {
        uint32_t idx = static_cast<uint32_t>(out.positions.size() / 3);
        edgeVertIdx[key] = idx;
        out.positions.push_back(ep.pos[0]);
        out.positions.push_back(ep.pos[1]);
        out.positions.push_back(ep.pos[2]);
        if (hasUV) {
            out.uvs.push_back(ep.uv[0]);
            out.uvs.push_back(ep.uv[1]);
        }
        if (hasC) {
            out.colors.push_back(ep.col[0]);
            out.colors.push_back(ep.col[1]);
            out.colors.push_back(ep.col[2]);
            out.colors.push_back(ep.col[3]);
        }
    }

    // Add face point vertices
    std::vector<uint32_t> faceVertIdx(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t idx = static_cast<uint32_t>(out.positions.size() / 3);
        faceVertIdx[t] = idx;
        out.positions.push_back(facePoints[t * 3 + 0]);
        out.positions.push_back(facePoints[t * 3 + 1]);
        out.positions.push_back(facePoints[t * 3 + 2]);
        if (hasUV) {
            out.uvs.push_back(faceUVs[t * 2 + 0]);
            out.uvs.push_back(faceUVs[t * 2 + 1]);
        }
        if (hasC) {
            out.colors.push_back(faceColors[t * 4 + 0]);
            out.colors.push_back(faceColors[t * 4 + 1]);
            out.colors.push_back(faceColors[t * 4 + 2]);
            out.colors.push_back(faceColors[t * 4 + 3]);
        }
    }

    // Build triangulated quads for each original triangle
    out.indices.reserve(triCount * 6 * 3); // 3 quads * 2 tris each
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v0 = mesh.indices[t * 3 + 0];
        uint32_t v1 = mesh.indices[t * 3 + 1];
        uint32_t v2 = mesh.indices[t * 3 + 2];

        uint32_t fp = faceVertIdx[t];
        uint32_t e01 = edgeVertIdx[EdgeKey(v0, v1)];
        uint32_t e12 = edgeVertIdx[EdgeKey(v1, v2)];
        uint32_t e20 = edgeVertIdx[EdgeKey(v2, v0)];

        // Quad: v0, e01, fp, e20 -> 2 triangles
        out.indices.push_back(v0);  out.indices.push_back(e01); out.indices.push_back(fp);
        out.indices.push_back(v0);  out.indices.push_back(fp);  out.indices.push_back(e20);

        // Quad: v1, e12, fp, e01 -> 2 triangles
        out.indices.push_back(v1);  out.indices.push_back(e12); out.indices.push_back(fp);
        out.indices.push_back(v1);  out.indices.push_back(fp);  out.indices.push_back(e01);

        // Quad: v2, e20, fp, e12 -> 2 triangles
        out.indices.push_back(v2);  out.indices.push_back(e20); out.indices.push_back(fp);
        out.indices.push_back(v2);  out.indices.push_back(fp);  out.indices.push_back(e12);
    }

    // Compute normals
    size_t outVertCount = out.positions.size() / 3;
    out.normals.assign(outVertCount * 3, 0.0f);
    size_t newTriCount = out.indices.size() / 3;
    for (size_t t = 0; t < newTriCount; ++t) {
        uint32_t i0 = out.indices[t * 3 + 0];
        uint32_t i1 = out.indices[t * 3 + 1];
        uint32_t i2 = out.indices[t * 3 + 2];
        float e1x = out.positions[i1*3+0] - out.positions[i0*3+0];
        float e1y = out.positions[i1*3+1] - out.positions[i0*3+1];
        float e1z = out.positions[i1*3+2] - out.positions[i0*3+2];
        float e2x = out.positions[i2*3+0] - out.positions[i0*3+0];
        float e2y = out.positions[i2*3+1] - out.positions[i0*3+1];
        float e2z = out.positions[i2*3+2] - out.positions[i0*3+2];
        float nx = e1y*e2z - e1z*e2y;
        float ny = e1z*e2x - e1x*e2z;
        float nz = e1x*e2y - e1y*e2x;
        for (uint32_t idx : {i0, i1, i2}) {
            out.normals[idx*3+0] += nx;
            out.normals[idx*3+1] += ny;
            out.normals[idx*3+2] += nz;
        }
    }
    for (size_t v = 0; v < outVertCount; ++v) {
        float* n = &out.normals[v * 3];
        float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (len > 1e-8f) { n[0] /= len; n[1] /= len; n[2] /= len; }
    }

    return out;
}

MeshData subdivideCatmullClark(const MeshData& mesh, int iterations) {
    if (iterations <= 0) return mesh;
    MeshData result = mesh;
    for (int i = 0; i < iterations; ++i)
        result = subdivideCCOnce(result);
    return result;
}

} // namespace bromesh
