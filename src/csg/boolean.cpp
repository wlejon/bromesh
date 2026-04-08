#include "bromesh/csg/boolean.h"

#ifdef BROMESH_HAS_MANIFOLD
#include <manifold/manifold.h>
#include <unordered_map>
#include <cmath>
#include <cstring>
#endif

namespace bromesh {

#ifdef BROMESH_HAS_MANIFOLD

// Weld coincident positions and build a manifold-friendly MeshGL.
// Our primitives often have split vertices for normals/UVs, but manifold
// needs shared vertices for topology. We merge by position here.
static manifold::MeshGL toMeshGL(const MeshData& mesh) {
    // Quantize positions to find duplicates
    struct PosKey {
        int32_t x, y, z;
        bool operator==(const PosKey& o) const { return x == o.x && y == o.y && z == o.z; }
    };
    struct PosKeyHash {
        size_t operator()(const PosKey& k) const {
            return std::hash<int32_t>()(k.x) * 73856093ULL ^
                   std::hash<int32_t>()(k.y) * 19349663ULL ^
                   std::hash<int32_t>()(k.z) * 83492791ULL;
        }
    };

    const float quantScale = 1e5f;
    size_t origVertCount = mesh.vertexCount();

    std::unordered_map<PosKey, uint32_t, PosKeyHash> posMap;
    std::vector<uint32_t> remap(origVertCount);
    std::vector<float> uniquePositions;
    uint32_t uniqueCount = 0;

    for (size_t v = 0; v < origVertCount; ++v) {
        PosKey key{
            static_cast<int32_t>(std::round(mesh.positions[v * 3 + 0] * quantScale)),
            static_cast<int32_t>(std::round(mesh.positions[v * 3 + 1] * quantScale)),
            static_cast<int32_t>(std::round(mesh.positions[v * 3 + 2] * quantScale))
        };
        auto it = posMap.find(key);
        if (it != posMap.end()) {
            remap[v] = it->second;
        } else {
            remap[v] = uniqueCount;
            posMap[key] = uniqueCount;
            uniquePositions.push_back(mesh.positions[v * 3 + 0]);
            uniquePositions.push_back(mesh.positions[v * 3 + 1]);
            uniquePositions.push_back(mesh.positions[v * 3 + 2]);
            uniqueCount++;
        }
    }

    manifold::MeshGL gl;
    gl.numProp = 3;
    gl.vertProperties = std::move(uniquePositions);

    gl.triVerts.resize(mesh.indices.size());
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        gl.triVerts[i] = remap[mesh.indices[i]];
    }

    return gl;
}

// Convert manifold::MeshGL back to MeshData
static MeshData fromMeshGL(const manifold::MeshGL& gl) {
    MeshData mesh;
    size_t vertCount = gl.NumVert();
    size_t triCount = gl.NumTri();

    mesh.positions.resize(vertCount * 3);
    for (size_t v = 0; v < vertCount; ++v) {
        for (int c = 0; c < 3; ++c) {
            mesh.positions[v * 3 + c] = gl.vertProperties[v * gl.numProp + c];
        }
    }

    mesh.indices.resize(triCount * 3);
    for (size_t i = 0; i < triCount * 3; ++i) {
        mesh.indices[i] = gl.triVerts[i];
    }

    return mesh;
}

static manifold::OpType toOpType(BooleanOp op) {
    switch (op) {
        case BooleanOp::Union:        return manifold::OpType::Add;
        case BooleanOp::Difference:   return manifold::OpType::Subtract;
        case BooleanOp::Intersection: return manifold::OpType::Intersect;
    }
    return manifold::OpType::Add;
}

#endif // BROMESH_HAS_MANIFOLD

MeshData booleanOp(const MeshData& a, const MeshData& b, BooleanOp op) {
#ifdef BROMESH_HAS_MANIFOLD
    if (a.empty() || b.empty()) return {};

    try {
        manifold::Manifold mA(toMeshGL(a));
        if (mA.Status() != manifold::Manifold::Error::NoError) return {};

        manifold::Manifold mB(toMeshGL(b));
        if (mB.Status() != manifold::Manifold::Error::NoError) return {};

        manifold::Manifold result = mA.Boolean(mB, toOpType(op));
        if (result.Status() != manifold::Manifold::Error::NoError) return {};

        return fromMeshGL(result.GetMeshGL());
    } catch (...) {
        return {};
    }
#else
    (void)a; (void)b; (void)op;
    return {};
#endif
}

MeshData booleanUnion(const MeshData& a, const MeshData& b) {
    return booleanOp(a, b, BooleanOp::Union);
}

MeshData booleanDifference(const MeshData& a, const MeshData& b) {
    return booleanOp(a, b, BooleanOp::Difference);
}

MeshData booleanIntersection(const MeshData& a, const MeshData& b) {
    return booleanOp(a, b, BooleanOp::Intersection);
}

std::pair<MeshData, MeshData> splitByPlane(const MeshData& mesh,
                                           float nx, float ny, float nz,
                                           float offset) {
#ifdef BROMESH_HAS_MANIFOLD
    if (mesh.empty()) return {};

    try {
        manifold::Manifold m(toMeshGL(mesh));
        if (m.Status() != manifold::Manifold::Error::NoError) return {};

        auto [top, bottom] = m.SplitByPlane(
            manifold::vec3(nx, ny, nz), offset);

        MeshData topMesh, bottomMesh;

        if (top.Status() == manifold::Manifold::Error::NoError &&
            top.NumTri() > 0) {
            topMesh = fromMeshGL(top.GetMeshGL());
        }
        if (bottom.Status() == manifold::Manifold::Error::NoError &&
            bottom.NumTri() > 0) {
            bottomMesh = fromMeshGL(bottom.GetMeshGL());
        }

        return {topMesh, bottomMesh};
    } catch (...) {
        return {};
    }
#else
    (void)mesh; (void)nx; (void)ny; (void)nz; (void)offset;
    return {};
#endif
}

} // namespace bromesh
