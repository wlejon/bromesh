#include "bromesh/manipulation/repair.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bromesh {

// Hash a 3-float position quantized to a grid for fast lookup
static uint64_t quantizePos(const float* p, float epsilon) {
    auto q = [epsilon](float v) -> int32_t {
        return static_cast<int32_t>(std::round(v / epsilon));
    };
    uint64_t x = static_cast<uint32_t>(q(p[0]));
    uint64_t y = static_cast<uint32_t>(q(p[1]));
    uint64_t z = static_cast<uint32_t>(q(p[2]));
    return (x * 73856093ULL) ^ (y * 19349663ULL) ^ (z * 83492791ULL);
}

MeshData removeDegenerateTriangles(const MeshData& mesh, float areaEpsilon) {
    if (mesh.empty() || mesh.indices.empty()) return mesh;

    // Identify non-degenerate triangles
    std::vector<bool> keep(mesh.triangleCount(), true);
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        // Degenerate if any two indices are the same
        if (i0 == i1 || i1 == i2 || i0 == i2) {
            keep[t] = false;
            continue;
        }

        // Check area via cross product magnitude
        const float* p0 = &mesh.positions[i0 * 3];
        const float* p1 = &mesh.positions[i1 * 3];
        const float* p2 = &mesh.positions[i2 * 3];

        float e1x = p1[0] - p0[0], e1y = p1[1] - p0[1], e1z = p1[2] - p0[2];
        float e2x = p2[0] - p0[0], e2y = p2[1] - p0[1], e2z = p2[2] - p0[2];

        float cx = e1y * e2z - e1z * e2y;
        float cy = e1z * e2x - e1x * e2z;
        float cz = e1x * e2y - e1y * e2x;

        float areaSq = cx * cx + cy * cy + cz * cz;
        if (areaSq < areaEpsilon * areaEpsilon) {
            keep[t] = false;
        }
    }

    // Build new index buffer with only kept triangles
    MeshData result = mesh; // copy all vertex data
    result.indices.clear();
    result.indices.reserve(mesh.indices.size());

    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        if (keep[t]) {
            result.indices.push_back(mesh.indices[t * 3 + 0]);
            result.indices.push_back(mesh.indices[t * 3 + 1]);
            result.indices.push_back(mesh.indices[t * 3 + 2]);
        }
    }

    return result;
}

MeshData removeDuplicateTriangles(const MeshData& mesh) {
    if (mesh.empty() || mesh.indices.empty()) return mesh;

    // For each triangle, create a canonical key from sorted quantized vertex positions
    struct TriKey {
        uint64_t h0, h1, h2;
        bool operator==(const TriKey& o) const { return h0 == o.h0 && h1 == o.h1 && h2 == o.h2; }
    };

    struct TriKeyHash {
        size_t operator()(const TriKey& k) const {
            return std::hash<uint64_t>()(k.h0) ^
                   (std::hash<uint64_t>()(k.h1) * 2654435761ULL) ^
                   (std::hash<uint64_t>()(k.h2) * 40343ULL);
        }
    };

    std::unordered_set<TriKey, TriKeyHash> seen;
    MeshData result = mesh;
    result.indices.clear();
    result.indices.reserve(mesh.indices.size());

    float eps = 1e-6f;
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        uint64_t h0 = quantizePos(&mesh.positions[i0 * 3], eps);
        uint64_t h1 = quantizePos(&mesh.positions[i1 * 3], eps);
        uint64_t h2 = quantizePos(&mesh.positions[i2 * 3], eps);

        // Sort hashes for canonical form
        if (h0 > h1) std::swap(h0, h1);
        if (h1 > h2) std::swap(h1, h2);
        if (h0 > h1) std::swap(h0, h1);

        TriKey key{h0, h1, h2};
        if (seen.insert(key).second) {
            result.indices.push_back(i0);
            result.indices.push_back(i1);
            result.indices.push_back(i2);
        }
    }

    return result;
}

MeshData fillHoles(const MeshData& mesh, int maxEdges) {
    if (mesh.empty() || mesh.indices.empty()) return mesh;

    // Build half-edge structure: for each directed edge (a->b), record which triangle uses it
    // A boundary edge has only one triangle (no opposite half-edge).
    struct Edge {
        uint32_t a, b;
        bool operator==(const Edge& o) const { return a == o.a && b == o.b; }
    };
    struct EdgeHash {
        size_t operator()(const Edge& e) const {
            return std::hash<uint64_t>()(static_cast<uint64_t>(e.a) << 32 | e.b);
        }
    };

    std::unordered_set<Edge, EdgeHash> halfEdges;
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];
        halfEdges.insert({i0, i1});
        halfEdges.insert({i1, i2});
        halfEdges.insert({i2, i0});
    }

    // Find boundary edges (a->b exists but b->a doesn't)
    // Store as adjacency: for each boundary vertex, where does the boundary go next?
    std::unordered_map<uint32_t, uint32_t> boundaryNext;
    for (const auto& e : halfEdges) {
        if (halfEdges.find({e.b, e.a}) == halfEdges.end()) {
            // e is a boundary edge; the boundary loop goes in the opposite direction (b->a)
            boundaryNext[e.b] = e.a;
        }
    }

    if (boundaryNext.empty()) return mesh; // No holes

    // Trace boundary loops
    std::unordered_set<uint32_t> visited;
    std::vector<std::vector<uint32_t>> loops;

    for (auto& [start, _] : boundaryNext) {
        if (visited.count(start)) continue;

        std::vector<uint32_t> loop;
        uint32_t cur = start;
        bool valid = true;
        while (true) {
            if (visited.count(cur)) {
                // If we've come back to start, we've completed the loop
                if (cur == start && !loop.empty()) break;
                // Otherwise this is a non-manifold boundary, skip
                valid = false;
                break;
            }
            visited.insert(cur);
            loop.push_back(cur);

            auto it = boundaryNext.find(cur);
            if (it == boundaryNext.end()) { valid = false; break; }
            cur = it->second;

            if (loop.size() > static_cast<size_t>(maxEdges)) { valid = false; break; }
        }

        if (valid && loop.size() >= 3 && static_cast<int>(loop.size()) <= maxEdges) {
            loops.push_back(std::move(loop));
        }
    }

    if (loops.empty()) return mesh;

    // Build result: copy original mesh and add fan triangles for each hole
    MeshData result = mesh;
    for (const auto& loop : loops) {
        // Simple fan triangulation from the first vertex
        uint32_t v0 = loop[0];
        for (size_t i = 1; i + 1 < loop.size(); ++i) {
            result.indices.push_back(v0);
            result.indices.push_back(loop[i]);
            result.indices.push_back(loop[i + 1]);
        }
    }

    return result;
}

} // namespace bromesh
