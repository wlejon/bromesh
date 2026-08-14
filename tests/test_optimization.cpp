#include "test_framework.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

struct Point3D {
    float x, y, z;
    bool operator<(const Point3D& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
    bool operator==(const Point3D& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

struct CanonicalTri {
    Point3D p0, p1, p2;

    static CanonicalTri make(Point3D a, Point3D b, Point3D c) {
        Point3D pts[3] = {a, b, c};
        int best = 0;
        if (pts[1] < pts[best]) best = 1;
        if (pts[2] < pts[best]) best = 2;
        return {pts[best], pts[(best + 1) % 3], pts[(best + 2) % 3]};
    }

    bool operator<(const CanonicalTri& o) const {
        if (!(p0 == o.p0)) return p0 < o.p0;
        if (!(p1 == o.p1)) return p1 < o.p1;
        return p2 < o.p2;
    }

    bool operator==(const CanonicalTri& o) const {
        return p0 == o.p0 && p1 == o.p1 && p2 == o.p2;
    }
};

static std::vector<CanonicalTri> extractTriangles(const bromesh::MeshData& mesh) {
    std::vector<CanonicalTri> tris;
    size_t triCount = mesh.triangleCount();
    tris.reserve(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];
        Point3D p0{mesh.positions[i0 * 3], mesh.positions[i0 * 3 + 1], mesh.positions[i0 * 3 + 2]};
        Point3D p1{mesh.positions[i1 * 3], mesh.positions[i1 * 3 + 1], mesh.positions[i1 * 3 + 2]};
        Point3D p2{mesh.positions[i2 * 3], mesh.positions[i2 * 3 + 1], mesh.positions[i2 * 3 + 2]};
        tris.push_back(CanonicalTri::make(p0, p1, p2));
    }
    return tris;
}

static std::vector<CanonicalTri> extractTriangles(const std::vector<float>& positions,
                                                  const std::vector<uint32_t>& indices) {
    std::vector<CanonicalTri> tris;
    size_t triCount = indices.size() / 3;
    tris.reserve(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = indices[t * 3 + 0];
        uint32_t i1 = indices[t * 3 + 1];
        uint32_t i2 = indices[t * 3 + 2];
        Point3D p0{positions[i0 * 3], positions[i0 * 3 + 1], positions[i0 * 3 + 2]};
        Point3D p1{positions[i1 * 3], positions[i1 * 3 + 1], positions[i1 * 3 + 2]};
        Point3D p2{positions[i2 * 3], positions[i2 * 3 + 1], positions[i2 * 3 + 2]};
        tris.push_back(CanonicalTri::make(p0, p1, p2));
    }
    return tris;
}

static bool haveSameTriangleSet(const bromesh::MeshData& a, const bromesh::MeshData& b) {
    if (a.triangleCount() != b.triangleCount()) return false;
    auto trisA = extractTriangles(a);
    auto trisB = extractTriangles(b);
    std::sort(trisA.begin(), trisA.end());
    std::sort(trisB.begin(), trisB.end());
    return trisA == trisB;
}

static bool haveSameTriangleSet(const bromesh::MeshData& a,
                                const std::vector<uint32_t>& indicesB) {
    if (a.triangleCount() != indicesB.size() / 3) return false;
    auto trisA = extractTriangles(a);
    auto trisB = extractTriangles(a.positions, indicesB);
    std::sort(trisA.begin(), trisA.end());
    std::sort(trisB.begin(), trisB.end());
    return trisA == trisB;
}

} // namespace

TEST(optimize_vertex_cache) {
    auto b = bromesh::box(1, 1, 1);
    auto orig = b;
    size_t origVerts = b.vertexCount();
    size_t origTris = b.triangleCount();
    bromesh::optimizeVertexCache(b);
    ASSERT(b.vertexCount() == origVerts, "vertex cache opt should preserve vertex count");
    ASSERT(b.triangleCount() == origTris, "vertex cache opt should preserve triangle count");
    ASSERT(haveSameTriangleSet(orig, b), "vertex cache opt should preserve exact geometric triangles");
}

TEST(optimize_vertex_fetch) {
    auto b = bromesh::box(1, 1, 1);
    auto orig = b;
    size_t origVerts = b.vertexCount();
    size_t origTris = b.triangleCount();
    bromesh::optimizeVertexFetch(b);
    ASSERT(b.vertexCount() == origVerts, "vertex fetch opt should preserve vertex count");
    ASSERT(b.triangleCount() == origTris, "vertex fetch opt should preserve triangle count");
    ASSERT(b.hasNormals(), "vertex fetch opt should preserve normals");
    ASSERT(b.hasUVs(), "vertex fetch opt should preserve UVs");
    ASSERT(haveSameTriangleSet(orig, b), "vertex fetch opt should preserve exact geometric triangles");
}


#if BROMESH_HAS_MESHOPTIMIZER
TEST(meshlets_sphere) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    auto meshlets = bromesh::buildMeshlets(mesh);
    ASSERT(!meshlets.empty(), "meshlets_sphere: should produce meshlets");
    ASSERT(meshlets.size() >= 2, "meshlets_sphere: sphere should have multiple meshlets");

    // Verify each meshlet has valid data
    size_t totalTris = 0;
    for (const auto& ml : meshlets) {
        ASSERT(ml.vertexCount() > 0, "meshlets_sphere: meshlet should have vertices");
        ASSERT(ml.triangleCount() > 0, "meshlets_sphere: meshlet should have triangles");
        ASSERT(ml.vertexCount() <= 64, "meshlets_sphere: meshlet should respect maxVertices");
        ASSERT(ml.triangleCount() <= 124, "meshlets_sphere: meshlet should respect maxTriangles");
        ASSERT(ml.bounds.radius > 0, "meshlets_sphere: meshlet should have bounding sphere");
        totalTris += ml.triangleCount();
    }
    ASSERT(totalTris == mesh.triangleCount(), "meshlets_sphere: total triangles should match");
}

TEST(meshlets_custom_params) {
    auto mesh = bromesh::box(1, 1, 1);
    bromesh::MeshletParams params;
    params.maxVertices = 32;
    params.maxTriangles = 32;
    auto meshlets = bromesh::buildMeshlets(mesh, params);
    ASSERT(!meshlets.empty(), "meshlets_custom: should produce meshlets");
    for (const auto& ml : meshlets) {
        ASSERT(ml.vertexCount() <= 32, "meshlets_custom: should respect custom maxVertices");
        ASSERT(ml.triangleCount() <= 32, "meshlets_custom: should respect custom maxTriangles");
    }
}

TEST(analyze_vertex_cache) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    auto stats = bromesh::analyzeVertexCache(mesh);
    ASSERT(stats.verticesTransformed > 0, "analyze_vcache: should transform vertices");
    ASSERT(stats.acmr > 0, "analyze_vcache: ACMR should be positive");
    ASSERT(stats.atvr >= 1.0f, "analyze_vcache: ATVR should be >= 1.0");
}

TEST(analyze_vertex_fetch) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    auto stats = bromesh::analyzeVertexFetch(mesh);
    ASSERT(stats.bytesFetched > 0, "analyze_vfetch: should fetch bytes");
    ASSERT(stats.overfetch >= 1.0f, "analyze_vfetch: overfetch should be >= 1.0");
}

TEST(analyze_overdraw) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    auto stats = bromesh::analyzeOverdraw(mesh);
    ASSERT(stats.pixelsCovered > 0, "analyze_overdraw: should cover pixels");
    ASSERT(stats.overdraw >= 1.0f, "analyze_overdraw: overdraw should be >= 1.0");
}

TEST(analyze_before_after_optimize) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    auto before = bromesh::analyzeVertexCache(mesh);

    bromesh::optimizeVertexCache(mesh);
    auto after = bromesh::analyzeVertexCache(mesh);

    // After optimization, ACMR should be equal or better (lower)
    ASSERT(after.acmr <= before.acmr + 0.01f,
           "analyze_opt: ACMR should improve after vertex cache optimization");
}

TEST(spatial_sort_triangles) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    size_t origTriCount = mesh.triangleCount();
    size_t origVertCount = mesh.vertexCount();
    bromesh::spatialSortTriangles(mesh);
    ASSERT(mesh.triangleCount() == origTriCount, "spatial_tri: triangle count preserved");
    ASSERT(mesh.vertexCount() == origVertCount, "spatial_tri: vertex count preserved");
}

TEST(spatial_sort_vertices) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    size_t origTriCount = mesh.triangleCount();
    size_t origVertCount = mesh.vertexCount();
    bromesh::spatialSortVertices(mesh);
    ASSERT(mesh.triangleCount() == origTriCount, "spatial_vert: triangle count preserved");
    ASSERT(mesh.vertexCount() == origVertCount, "spatial_vert: vertex count preserved");
    ASSERT(mesh.hasNormals(), "spatial_vert: normals preserved");
}

TEST(shadow_index_buffer) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    auto shadow = bromesh::generateShadowIndexBuffer(mesh);
    ASSERT(shadow.size() == mesh.indices.size(), "shadow_ib: same index count");
    // Shadow indices should reference valid vertices
    for (uint32_t idx : shadow) {
        ASSERT(idx < mesh.vertexCount(), "shadow_ib: valid vertex reference");
    }
}

TEST(encode_decode_mesh_roundtrip) {
    auto mesh = bromesh::sphere(2.0f, 16, 12);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto encoded = bromesh::encodeMesh(mesh);
    ASSERT(!encoded.vertexData.empty(), "encode: vertex data should be non-empty");
    ASSERT(!encoded.indexData.empty(), "encode: index data should be non-empty");
    ASSERT(encoded.vertexCount == mesh.vertexCount(), "encode: vertex count preserved");
    ASSERT(encoded.indexCount == mesh.indices.size(), "encode: index count preserved");

    // Compressed should be smaller than raw
    size_t rawVertSize = mesh.vertexCount() * encoded.vertexSize;
    ASSERT(encoded.vertexData.size() < rawVertSize, "encode: vertex data should be compressed");

    auto decoded = bromesh::decodeMesh(encoded, true, true, false);
    ASSERT(decoded.vertexCount() == mesh.vertexCount(), "decode: vertex count matches");
    ASSERT(decoded.triangleCount() == mesh.triangleCount(), "decode: triangle count matches");
    ASSERT(decoded.hasNormals(), "decode: normals preserved");
    ASSERT(decoded.hasUVs(), "decode: UVs preserved");

    // Verify positions match
    bool posMatch = true;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        if (std::fabs(mesh.positions[i] - decoded.positions[i]) > 1e-5f) {
            posMatch = false; break;
        }
    }
    ASSERT(posMatch, "decode: positions match original");
}

TEST(encode_decode_index_buffer) {
    auto mesh = bromesh::box(1, 1, 1);
    auto encoded = bromesh::encodeIndexBuffer(mesh.indices, mesh.vertexCount());
    ASSERT(!encoded.empty(), "encode_ib: should produce data");
    ASSERT(encoded.size() < mesh.indices.size() * sizeof(uint32_t), "encode_ib: should compress");

    auto decoded = bromesh::decodeIndexBuffer(encoded, mesh.indices.size());
    ASSERT(decoded.size() == mesh.indices.size(), "decode_ib: size matches");
    ASSERT(decoded == mesh.indices, "decode_ib: indices match original");
}

TEST(stripify_unstripify_roundtrip) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    bromesh::optimizeVertexCache(mesh);

    auto strip = bromesh::stripify(mesh.indices, mesh.vertexCount());
    ASSERT(!strip.empty(), "stripify: should produce strip");

    auto restored = bromesh::unstripify(strip);
    ASSERT(!restored.empty(), "unstripify: should produce triangle list");
    // Restored triangle count should match original
    ASSERT(restored.size() / 3 == mesh.indices.size() / 3,
           "strip_roundtrip: triangle count preserved");
    ASSERT(haveSameTriangleSet(mesh, restored),
           "strip_roundtrip: unstripify restores exact same triangle set");
}

TEST(stripify_box) {
    auto mesh = bromesh::box(1, 1, 1);
    auto strip = bromesh::stripify(mesh.indices, mesh.vertexCount());
    ASSERT(!strip.empty(), "stripify_box: should produce strip");
    // Strip should be reasonably compact
    ASSERT(strip.size() <= mesh.indices.size() * 2,
           "stripify_box: strip shouldn't be much larger than triangle list");
    auto restored = bromesh::unstripify(strip);
    ASSERT(haveSameTriangleSet(mesh, restored),
           "stripify_box: unstripify restores exact same triangle set");
}

#endif // BROMESH_HAS_MESHOPTIMIZER

TEST(progressive_mesh_build) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    bromesh::computeNormals(mesh);

    auto pm = bromesh::buildProgressiveMesh(mesh);
    ASSERT(!pm.collapses.empty(), "pm_build: should have collapse records");
    ASSERT(pm.maxTriangles() == mesh.triangleCount(), "pm_build: max tris matches original");
    ASSERT(pm.minTriangles() < pm.maxTriangles(), "pm_build: can simplify below max");
}

TEST(progressive_mesh_full_resolution) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);

    auto pm = bromesh::buildProgressiveMesh(mesh);
    auto full = bromesh::progressiveMeshAtTriangleCount(pm, pm.maxTriangles());

    ASSERT(full.triangleCount() == mesh.triangleCount(),
           "pm_full: full resolution matches original triangle count");
}

TEST(progressive_mesh_half_resolution) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);

    auto pm = bromesh::buildProgressiveMesh(mesh);
    size_t halfTri = mesh.triangleCount() / 2;
    auto half = bromesh::progressiveMeshAtTriangleCount(pm, halfTri);

    // Should be approximately half (may not be exact due to collapse granularity)
    ASSERT(half.triangleCount() <= halfTri + 10, "pm_half: roughly half triangles");
    ASSERT(half.triangleCount() > 0, "pm_half: not empty");
    ASSERT(!half.positions.empty(), "pm_half: has positions");
}

TEST(progressive_mesh_ratio) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);

    auto pm = bromesh::buildProgressiveMesh(mesh);
    auto lod0 = bromesh::progressiveMeshAtRatio(pm, 1.0f);
    auto lod50 = bromesh::progressiveMeshAtRatio(pm, 0.5f);
    auto lodMin = bromesh::progressiveMeshAtRatio(pm, 0.0f);

    ASSERT(lod0.triangleCount() >= lod50.triangleCount(),
           "pm_ratio: 100% >= 50%");
    ASSERT(lod50.triangleCount() >= lodMin.triangleCount(),
           "pm_ratio: 50% >= 0%");
}

TEST(progressive_mesh_monotonic_decrease) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);

    auto pm = bromesh::buildProgressiveMesh(mesh);

    // Extract at several levels and verify monotonically decreasing triangle counts
    size_t prev = pm.maxTriangles();
    bool monotonic = true;
    for (float r = 0.9f; r >= 0.1f; r -= 0.1f) {
        auto lod = bromesh::progressiveMeshAtRatio(pm, r);
        if (lod.triangleCount() > prev) {
            monotonic = false; break;
        }
        prev = lod.triangleCount();
    }
    ASSERT(monotonic, "pm_monotonic: triangle count decreases with ratio");
}

TEST(progressive_mesh_serialize_roundtrip) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);

    auto pm = bromesh::buildProgressiveMesh(mesh);
    auto data = bromesh::serializeProgressiveMesh(pm);
    ASSERT(!data.empty(), "pm_serialize: data not empty");

    auto pm2 = bromesh::deserializeProgressiveMesh(data.data(), data.size());
    ASSERT(pm2.maxTriangles() == pm.maxTriangles(), "pm_roundtrip: same max triangles");
    ASSERT(pm2.collapses.size() == pm.collapses.size(), "pm_roundtrip: same collapse count");

    // Extract at half resolution from both and compare triangle counts
    size_t halfTri = pm.maxTriangles() / 2;
    auto lodA = bromesh::progressiveMeshAtTriangleCount(pm, halfTri);
    auto lodB = bromesh::progressiveMeshAtTriangleCount(pm2, halfTri);
    ASSERT(lodA.triangleCount() == lodB.triangleCount(),
           "pm_roundtrip: same LOD output");
}

TEST(progressive_mesh_bad_deserialize) {
    // Should not crash on garbage data
    uint8_t garbage[32] = {};
    auto pm = bromesh::deserializeProgressiveMesh(garbage, sizeof(garbage));
    ASSERT(pm.baseMesh.empty(), "pm_bad_data: returns empty on garbage");

    auto pm2 = bromesh::deserializeProgressiveMesh(nullptr, 0);
    ASSERT(pm2.baseMesh.empty(), "pm_null: returns empty on null");
}

TEST(progressive_mesh_preserves_attributes) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto pm = bromesh::buildProgressiveMesh(mesh);
    auto half = bromesh::progressiveMeshAtRatio(pm, 0.5f);

    ASSERT(half.hasNormals(), "pm_attrs: preserves normals");
    ASSERT(half.hasUVs(), "pm_attrs: preserves UVs");
}

