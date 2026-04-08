#include "bromesh/optimization/progressive.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bromesh {

// ---- Quadric Error Metric (QEM) implementation ----

// A symmetric 4x4 matrix stored as 10 unique values (upper triangle).
struct Quadric {
    double a[10] = {};  // a00, a01, a02, a03, a11, a12, a13, a22, a23, a33

    Quadric() = default;

    Quadric(double a00, double a01, double a02, double a03,
            double a11, double a12, double a13,
            double a22, double a23, double a33) {
        a[0]=a00; a[1]=a01; a[2]=a02; a[3]=a03;
        a[4]=a11; a[5]=a12; a[6]=a13;
        a[7]=a22; a[8]=a23; a[9]=a33;
    }

    // Build quadric from a plane equation ax+by+cz+d=0
    static Quadric fromPlane(double nx, double ny, double nz, double d) {
        return Quadric(
            nx*nx, nx*ny, nx*nz, nx*d,
                   ny*ny, ny*nz, ny*d,
                          nz*nz, nz*d,
                                 d*d
        );
    }

    Quadric operator+(const Quadric& o) const {
        Quadric r;
        for (int i = 0; i < 10; ++i) r.a[i] = a[i] + o.a[i];
        return r;
    }

    Quadric& operator+=(const Quadric& o) {
        for (int i = 0; i < 10; ++i) a[i] += o.a[i];
        return *this;
    }

    // Evaluate v^T Q v for point (x,y,z)
    double evaluate(double x, double y, double z) const {
        return a[0]*x*x + 2*a[1]*x*y + 2*a[2]*x*z + 2*a[3]*x
             + a[4]*y*y + 2*a[5]*y*z + 2*a[6]*y
             + a[7]*z*z + 2*a[8]*z
             + a[9];
    }

    // Try to find the optimal point by solving the 3x3 linear system.
    // Returns false if the system is singular (use midpoint instead).
    bool optimalPoint(double& ox, double& oy, double& oz) const {
        // Solve [a00 a01 a02; a01 a11 a12; a02 a12 a22] * [x;y;z] = -[a03;a13;a23]
        double m00 = a[0], m01 = a[1], m02 = a[2];
        double m10 = a[1], m11 = a[4], m12 = a[5];
        double m20 = a[2], m21 = a[5], m22 = a[7];
        double b0 = -a[3], b1 = -a[6], b2 = -a[8];

        double det = m00*(m11*m22-m12*m21) - m01*(m10*m22-m12*m20) + m02*(m10*m21-m11*m20);
        if (std::fabs(det) < 1e-12) return false;

        double invDet = 1.0 / det;
        ox = invDet * (b0*(m11*m22-m12*m21) - m01*(b1*m22-m12*b2) + m02*(b1*m21-m11*b2));
        oy = invDet * (m00*(b1*m22-m12*b2) - b0*(m10*m22-m12*m20) + m02*(m10*b2-b1*m20));
        oz = invDet * (m00*(m11*b2-b1*m21) - m01*(m10*b2-b1*m20) + b0*(m10*m21-m11*m20));
        return true;
    }
};

// ---- Edge pair hash ----

struct EdgeKey {
    uint32_t v0, v1;
    bool operator==(const EdgeKey& o) const { return v0 == o.v0 && v1 == o.v1; }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& e) const {
        return std::hash<uint64_t>()(((uint64_t)e.v0 << 32) | e.v1);
    }
};


// ---- Priority queue entry ----

struct CollapseCandidate {
    double error;
    uint32_t v0, v1;
    float targetPos[3];
    bool operator>(const CollapseCandidate& o) const { return error > o.error; }
};

// ---- Progressive mesh builder ----

ProgressiveMesh buildProgressiveMesh(const MeshData& mesh) {
    ProgressiveMesh pm;
    if (mesh.empty() || mesh.indices.empty()) return pm;

    pm.baseMesh = mesh; // Store full mesh

    const size_t vertCount = mesh.vertexCount();
    const size_t triCount = mesh.triangleCount();

    // Working copies of mesh data
    std::vector<float> positions(mesh.positions);
    std::vector<uint32_t> indices(mesh.indices);

    // Per-vertex quadrics
    std::vector<Quadric> quadrics(vertCount);

    // Build per-face quadrics and accumulate to vertices
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = indices[t*3+0];
        uint32_t i1 = indices[t*3+1];
        uint32_t i2 = indices[t*3+2];

        float* p0 = &positions[i0*3];
        float* p1 = &positions[i1*3];
        float* p2 = &positions[i2*3];

        float e1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        float e2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        double nx = e1[1]*e2[2] - e1[2]*e2[1];
        double ny = e1[2]*e2[0] - e1[0]*e2[2];
        double nz = e1[0]*e2[1] - e1[1]*e2[0];
        double len = std::sqrt(nx*nx + ny*ny + nz*nz);
        if (len > 1e-12) { nx/=len; ny/=len; nz/=len; }

        double d = -(nx*p0[0] + ny*p0[1] + nz*p0[2]);
        Quadric q = Quadric::fromPlane(nx, ny, nz, d);

        quadrics[i0] += q;
        quadrics[i1] += q;
        quadrics[i2] += q;
    }

    // Build adjacency: per-vertex set of neighboring vertices and triangles
    std::vector<std::unordered_set<uint32_t>> vertNeighbors(vertCount);
    std::vector<std::unordered_set<uint32_t>> vertTris(vertCount); // triangles touching vertex

    // Track which triangles and vertices are alive
    std::vector<bool> triAlive(triCount, true);
    std::vector<bool> vertAlive(vertCount, true);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {indices[t*3], indices[t*3+1], indices[t*3+2]};
        for (int i = 0; i < 3; ++i) {
            vertTris[v[i]].insert((uint32_t)t);
            vertNeighbors[v[i]].insert(v[(i+1)%3]);
            vertNeighbors[v[i]].insert(v[(i+2)%3]);
        }
    }

    // Evaluate collapse cost for an edge
    auto evaluateEdge = [&](uint32_t a, uint32_t b) -> CollapseCandidate {
        CollapseCandidate c;
        c.v0 = a; c.v1 = b;

        Quadric q = quadrics[a] + quadrics[b];

        double ox, oy, oz;
        if (q.optimalPoint(ox, oy, oz)) {
            c.error = q.evaluate(ox, oy, oz);
            c.targetPos[0] = (float)ox;
            c.targetPos[1] = (float)oy;
            c.targetPos[2] = (float)oz;
        } else {
            // Fall back to midpoint
            float mx = (positions[a*3]+positions[b*3])*0.5f;
            float my = (positions[a*3+1]+positions[b*3+1])*0.5f;
            float mz = (positions[a*3+2]+positions[b*3+2])*0.5f;
            c.error = q.evaluate(mx, my, mz);
            c.targetPos[0] = mx;
            c.targetPos[1] = my;
            c.targetPos[2] = mz;
        }

        if (c.error < 0) c.error = 0;
        return c;
    };

    // Build initial edge set and priority queue
    // Use a vector-based heap
    std::vector<CollapseCandidate> heap;
    std::unordered_set<uint64_t> edgeSet;

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {indices[t*3], indices[t*3+1], indices[t*3+2]};
        for (int i = 0; i < 3; ++i) {
            uint32_t a = v[i], b = v[(i+1)%3];
            uint64_t key = ((uint64_t)std::min(a,b) << 32) | std::max(a,b);
            if (edgeSet.insert(key).second) {
                heap.push_back(evaluateEdge(std::min(a,b), std::max(a,b)));
            }
        }
    }

    std::make_heap(heap.begin(), heap.end(), std::greater<CollapseCandidate>());

    size_t currentTriCount = triCount;

    // Iteratively collapse edges
    while (!heap.empty() && currentTriCount > 4) {
        std::pop_heap(heap.begin(), heap.end(), std::greater<CollapseCandidate>());
        CollapseCandidate best = heap.back();
        heap.pop_back();

        uint32_t va = best.v0, vb = best.v1;

        // Skip if either vertex is dead
        if (!vertAlive[va] || !vertAlive[vb]) continue;

        // Skip if they're no longer neighbors (stale entry)
        if (vertNeighbors[va].find(vb) == vertNeighbors[va].end()) continue;

        // Count triangles that will be removed (those containing both va and vb)
        uint32_t removedCount = 0;
        std::vector<uint32_t> removedTris;
        for (uint32_t t : vertTris[va]) {
            if (!triAlive[t]) continue;
            uint32_t i0 = indices[t*3], i1 = indices[t*3+1], i2 = indices[t*3+2];
            if (i0 == vb || i1 == vb || i2 == vb) {
                removedTris.push_back(t);
                removedCount++;
            }
        }

        if (removedCount == 0) continue;

        // Record the collapse (va collapses into vb)
        CollapseRecord rec;
        rec.vertexFrom = va;
        rec.vertexTo = vb;
        rec.targetPosition[0] = best.targetPos[0];
        rec.targetPosition[1] = best.targetPos[1];
        rec.targetPosition[2] = best.targetPos[2];
        rec.numTrianglesRemoved = removedCount;
        pm.collapses.push_back(rec);

        // Kill the degenerate triangles
        for (uint32_t t : removedTris) {
            triAlive[t] = false;
            currentTriCount--;
            // Remove triangle from vertex adjacency
            uint32_t v[3] = {indices[t*3], indices[t*3+1], indices[t*3+2]};
            for (int i = 0; i < 3; ++i) {
                vertTris[v[i]].erase(t);
            }
        }

        // Remap all references from va to vb
        for (uint32_t t : vertTris[va]) {
            if (!triAlive[t]) continue;
            for (int i = 0; i < 3; ++i) {
                if (indices[t*3+i] == va) {
                    indices[t*3+i] = vb;
                }
            }
            vertTris[vb].insert(t);

            // Check for degenerate triangles after remap
            uint32_t i0 = indices[t*3], i1 = indices[t*3+1], i2 = indices[t*3+2];
            if (i0 == i1 || i1 == i2 || i0 == i2) {
                triAlive[t] = false;
                currentTriCount--;
                vertTris[i0].erase(t);
                vertTris[i1].erase(t);
                vertTris[i2].erase(t);
            }
        }

        // Update vertex position to optimal
        positions[vb*3+0] = best.targetPos[0];
        positions[vb*3+1] = best.targetPos[1];
        positions[vb*3+2] = best.targetPos[2];

        // Merge quadrics
        quadrics[vb] += quadrics[va];

        // Transfer neighbors
        for (uint32_t n : vertNeighbors[va]) {
            if (n == vb) continue;
            if (!vertAlive[n]) continue;
            vertNeighbors[n].erase(va);
            vertNeighbors[n].insert(vb);
            vertNeighbors[vb].insert(n);
        }

        // Kill va
        vertAlive[va] = false;
        vertNeighbors[va].clear();
        vertTris[va].clear();

        // Re-evaluate affected edges (all edges incident to vb)
        for (uint32_t n : vertNeighbors[vb]) {
            if (!vertAlive[n]) continue;
            auto cand = evaluateEdge(std::min(vb, n), std::max(vb, n));
            heap.push_back(cand);
            std::push_heap(heap.begin(), heap.end(), std::greater<CollapseCandidate>());
        }
    }

    return pm;
}

// ---- LOD extraction ----

MeshData progressiveMeshAtTriangleCount(const ProgressiveMesh& pm,
                                         size_t targetTriangles) {
    if (pm.baseMesh.empty()) return {};

    size_t maxTri = pm.maxTriangles();
    size_t minTri = pm.minTriangles();
    if (targetTriangles >= maxTri) return pm.baseMesh;
    if (targetTriangles < minTri) targetTriangles = minTri;

    // Replay collapses on a working copy
    MeshData mesh = pm.baseMesh;
    const size_t vertCount = mesh.vertexCount();

    // Build vertex remap and position overrides
    std::vector<uint32_t> remap(vertCount);
    for (uint32_t i = 0; i < vertCount; ++i) remap[i] = i;

    auto resolveVertex = [&](uint32_t v) -> uint32_t {
        while (remap[v] != v) v = remap[v];
        return v;
    };

    size_t currentTri = maxTri;

    for (const auto& col : pm.collapses) {
        if (currentTri <= targetTriangles) break;

        uint32_t from = col.vertexFrom;
        uint32_t to = col.vertexTo;

        remap[from] = to;

        // Update target position
        mesh.positions[to*3+0] = col.targetPosition[0];
        mesh.positions[to*3+1] = col.targetPosition[1];
        mesh.positions[to*3+2] = col.targetPosition[2];

        currentTri -= col.numTrianglesRemoved;
    }

    // Apply remap to indices and remove degenerate triangles
    std::vector<uint32_t> newIndices;
    newIndices.reserve(mesh.indices.size());

    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        uint32_t i0 = resolveVertex(mesh.indices[t*3+0]);
        uint32_t i1 = resolveVertex(mesh.indices[t*3+1]);
        uint32_t i2 = resolveVertex(mesh.indices[t*3+2]);

        // Skip degenerate triangles
        if (i0 == i1 || i1 == i2 || i0 == i2) continue;

        newIndices.push_back(i0);
        newIndices.push_back(i1);
        newIndices.push_back(i2);
    }

    // Build compacted output
    std::vector<uint32_t> vertRemap(vertCount, UINT32_MAX);
    uint32_t newVertCount = 0;
    for (uint32_t idx : newIndices) {
        if (vertRemap[idx] == UINT32_MAX) {
            vertRemap[idx] = newVertCount++;
        }
    }

    MeshData result;
    result.positions.resize(newVertCount * 3);
    if (mesh.hasNormals()) result.normals.resize(newVertCount * 3);
    if (mesh.hasUVs()) result.uvs.resize(newVertCount * 2);
    if (mesh.hasColors()) result.colors.resize(newVertCount * 4);

    for (uint32_t old = 0; old < vertCount; ++old) {
        if (vertRemap[old] == UINT32_MAX) continue;
        uint32_t nv = vertRemap[old];
        result.positions[nv*3+0] = mesh.positions[old*3+0];
        result.positions[nv*3+1] = mesh.positions[old*3+1];
        result.positions[nv*3+2] = mesh.positions[old*3+2];
        if (mesh.hasNormals()) {
            result.normals[nv*3+0] = mesh.normals[old*3+0];
            result.normals[nv*3+1] = mesh.normals[old*3+1];
            result.normals[nv*3+2] = mesh.normals[old*3+2];
        }
        if (mesh.hasUVs()) {
            result.uvs[nv*2+0] = mesh.uvs[old*2+0];
            result.uvs[nv*2+1] = mesh.uvs[old*2+1];
        }
        if (mesh.hasColors()) {
            result.colors[nv*4+0] = mesh.colors[old*4+0];
            result.colors[nv*4+1] = mesh.colors[old*4+1];
            result.colors[nv*4+2] = mesh.colors[old*4+2];
            result.colors[nv*4+3] = mesh.colors[old*4+3];
        }
    }

    result.indices.resize(newIndices.size());
    for (size_t i = 0; i < newIndices.size(); ++i) {
        result.indices[i] = vertRemap[newIndices[i]];
    }

    return result;
}

MeshData progressiveMeshAtRatio(const ProgressiveMesh& pm, float ratio) {
    if (pm.baseMesh.empty()) return {};
    ratio = std::max(0.0f, std::min(1.0f, ratio));
    size_t maxTri = pm.maxTriangles();
    size_t minTri = pm.minTriangles();
    size_t target = minTri + (size_t)((maxTri - minTri) * ratio);
    return progressiveMeshAtTriangleCount(pm, target);
}

size_t ProgressiveMesh::minTriangles() const {
    size_t total = maxTriangles();
    for (const auto& c : collapses) {
        if (total <= c.numTrianglesRemoved) break;
        total -= c.numTrianglesRemoved;
    }
    return total;
}

// ---- Serialization ----

static void writeU32(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(v & 0xFF);
    buf.push_back((v >> 8) & 0xFF);
    buf.push_back((v >> 16) & 0xFF);
    buf.push_back((v >> 24) & 0xFF);
}

static void writeF32(std::vector<uint8_t>& buf, float v) {
    uint32_t u;
    std::memcpy(&u, &v, 4);
    writeU32(buf, u);
}

static uint32_t readU32(const uint8_t*& ptr) {
    uint32_t v = (uint32_t)ptr[0] | ((uint32_t)ptr[1] << 8) |
                 ((uint32_t)ptr[2] << 16) | ((uint32_t)ptr[3] << 24);
    ptr += 4;
    return v;
}

static float readF32(const uint8_t*& ptr) {
    uint32_t u = readU32(ptr);
    float v;
    std::memcpy(&v, &u, 4);
    return v;
}

std::vector<uint8_t> serializeProgressiveMesh(const ProgressiveMesh& pm) {
    std::vector<uint8_t> buf;
    if (pm.baseMesh.empty()) return buf;

    // Magic + version
    writeU32(buf, 0x50524D53); // 'PRMS'
    writeU32(buf, 1);          // version

    const auto& m = pm.baseMesh;
    uint32_t vertCount = (uint32_t)m.vertexCount();
    uint32_t idxCount = (uint32_t)m.indices.size();
    uint32_t collapseCount = (uint32_t)pm.collapses.size();

    // Flags: bit 0 = hasNormals, bit 1 = hasUVs, bit 2 = hasColors
    uint32_t flags = 0;
    if (m.hasNormals()) flags |= 1;
    if (m.hasUVs()) flags |= 2;
    if (m.hasColors()) flags |= 4;

    writeU32(buf, vertCount);
    writeU32(buf, idxCount);
    writeU32(buf, collapseCount);
    writeU32(buf, flags);

    // Positions
    for (size_t i = 0; i < vertCount * 3; ++i) writeF32(buf, m.positions[i]);
    // Normals
    if (flags & 1) for (size_t i = 0; i < vertCount * 3; ++i) writeF32(buf, m.normals[i]);
    // UVs
    if (flags & 2) for (size_t i = 0; i < vertCount * 2; ++i) writeF32(buf, m.uvs[i]);
    // Colors
    if (flags & 4) for (size_t i = 0; i < vertCount * 4; ++i) writeF32(buf, m.colors[i]);
    // Indices
    for (size_t i = 0; i < idxCount; ++i) writeU32(buf, m.indices[i]);

    // Collapse records
    for (const auto& c : pm.collapses) {
        writeU32(buf, c.vertexFrom);
        writeU32(buf, c.vertexTo);
        writeF32(buf, c.targetPosition[0]);
        writeF32(buf, c.targetPosition[1]);
        writeF32(buf, c.targetPosition[2]);
        writeU32(buf, c.numTrianglesRemoved);
    }

    return buf;
}

ProgressiveMesh deserializeProgressiveMesh(const uint8_t* data, size_t size) {
    ProgressiveMesh pm;
    if (!data || size < 24) return pm;

    const uint8_t* ptr = data;
    uint32_t magic = readU32(ptr);
    uint32_t version = readU32(ptr);
    if (magic != 0x50524D53 || version != 1) return pm;

    uint32_t vertCount = readU32(ptr);
    uint32_t idxCount = readU32(ptr);
    uint32_t collapseCount = readU32(ptr);
    uint32_t flags = readU32(ptr);

    // Compute expected size
    size_t expectedSize = 24; // header: 6 x uint32
    expectedSize += vertCount * 3 * 4; // positions
    if (flags & 1) expectedSize += vertCount * 3 * 4; // normals
    if (flags & 2) expectedSize += vertCount * 2 * 4; // uvs
    if (flags & 4) expectedSize += vertCount * 4 * 4; // colors
    expectedSize += idxCount * 4; // indices
    expectedSize += collapseCount * 24; // collapses (6 * 4 bytes each)
    if (size < expectedSize) return pm;

    auto& m = pm.baseMesh;
    m.positions.resize(vertCount * 3);
    for (size_t i = 0; i < vertCount * 3; ++i) m.positions[i] = readF32(ptr);

    if (flags & 1) {
        m.normals.resize(vertCount * 3);
        for (size_t i = 0; i < vertCount * 3; ++i) m.normals[i] = readF32(ptr);
    }
    if (flags & 2) {
        m.uvs.resize(vertCount * 2);
        for (size_t i = 0; i < vertCount * 2; ++i) m.uvs[i] = readF32(ptr);
    }
    if (flags & 4) {
        m.colors.resize(vertCount * 4);
        for (size_t i = 0; i < vertCount * 4; ++i) m.colors[i] = readF32(ptr);
    }

    m.indices.resize(idxCount);
    for (size_t i = 0; i < idxCount; ++i) m.indices[i] = readU32(ptr);

    pm.collapses.resize(collapseCount);
    for (size_t i = 0; i < collapseCount; ++i) {
        pm.collapses[i].vertexFrom = readU32(ptr);
        pm.collapses[i].vertexTo = readU32(ptr);
        pm.collapses[i].targetPosition[0] = readF32(ptr);
        pm.collapses[i].targetPosition[1] = readF32(ptr);
        pm.collapses[i].targetPosition[2] = readF32(ptr);
        pm.collapses[i].numTrianglesRemoved = readU32(ptr);
    }

    return pm;
}

} // namespace bromesh
