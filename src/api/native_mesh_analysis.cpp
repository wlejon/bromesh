#include "host_mesh_internal.h"
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace bromesh::api {

HostClass g_meshBvhClass;
HostClass g_progressiveMeshClass;

namespace {

Value makeRayHitObject(const bromesh::RayHit& h, const bromesh::MeshData& m) {
    if (!h.hit) return ev::null();
    ObjectBuilder obj;
    obj.set("hit", true);
    obj.set("distance", static_cast<double>(h.distance));
    obj.set("triangle", static_cast<double>(h.triangleIndex));
    obj.set("triangleIndex", static_cast<double>(h.triangleIndex));

    const float pos[3] = {h.position[0], h.position[1], h.position[2]};
    {
        ev::Persistent p(makeFloat32Array(pos, 3));
        obj.set("point", p.get());
        obj.set("position", p.get());
    }
    const float norm[3] = {h.normal[0], h.normal[1], h.normal[2]};
    {
        ev::Persistent n(makeFloat32Array(norm, 3));
        obj.set("normal", n.get());
    }
    const float bary[3] = {h.baryU, h.baryV, h.baryW};
    {
        ev::Persistent b(makeFloat32Array(bary, 3));
        obj.set("barycentric", b.get());
    }
    if (m.hasUVs() && h.triangleIndex * 3 + 2 < m.indices.size()) {
        const uint32_t i0 = m.indices[h.triangleIndex * 3];
        const uint32_t i1 = m.indices[h.triangleIndex * 3 + 1];
        const uint32_t i2 = m.indices[h.triangleIndex * 3 + 2];
        const float u = h.baryU * m.uvs[i0 * 2] + h.baryV * m.uvs[i1 * 2] + h.baryW * m.uvs[i2 * 2];
        const float v = h.baryU * m.uvs[i0 * 2 + 1] + h.baryV * m.uvs[i1 * 2 + 1] + h.baryW * m.uvs[i2 * 2 + 1];
        const float uvArr[2] = {u, v};
        ev::Persistent uv(makeFloat32Array(uvArr, 2));
        obj.set("uv", uv.get());
    }
    return obj.build();
}

inline void parseRayArgs(std::span<const Value> a, float o[3], float d[3], float& maxDist) {
    if (!a.empty() && (ev::isObject(a[0]) || a.size() <= 3)) {
        std::vector<float> ov = toFloatVector(a[0]);
        if (ov.size() >= 3) { o[0] = ov[0]; o[1] = ov[1]; o[2] = ov[2]; }
        if (a.size() > 1) {
            std::vector<float> dv = toFloatVector(a[1]);
            if (dv.size() >= 3) { d[0] = dv[0]; d[1] = dv[1]; d[2] = dv[2]; }
        }
        if (a.size() > 2) maxDist = static_cast<float>(numAt(a, 2));
    } else {
        ArgReader r(a);
        o[0] = static_cast<float>(r.getDouble(0, 0.0));
        o[1] = static_cast<float>(r.getDouble(1, 0.0));
        o[2] = static_cast<float>(r.getDouble(2, 0.0));
        d[0] = static_cast<float>(r.getDouble(3, 0.0));
        d[1] = static_cast<float>(r.getDouble(4, -1.0));
        d[2] = static_cast<float>(r.getDouble(5, 0.0));
        maxDist = static_cast<float>(r.getDouble(6, 0.0));
    }
}

inline void parsePointArg(std::span<const Value> a, float p[3]) {
    if (!a.empty() && ev::isObject(a[0])) {
        std::vector<float> pv = toFloatVector(a[0]);
        if (pv.size() >= 3) { p[0] = pv[0]; p[1] = pv[1]; p[2] = pv[2]; }
    } else {
        ArgReader r(a);
        p[0] = static_cast<float>(r.getDouble(0, 0.0));
        p[1] = static_cast<float>(r.getDouble(1, 0.0));
        p[2] = static_cast<float>(r.getDouble(2, 0.0));
    }
}

int computeGenus(const bromesh::MeshData& mesh) {
    if (mesh.indices.empty()) return 0;
    auto qi = [](float f) -> int32_t { return static_cast<int32_t>(std::round(f * 10000.0f)); };
    struct Vec3iHash {
        size_t operator()(const std::tuple<int32_t,int32_t,int32_t>& v) const {
            size_t h = std::hash<int32_t>{}(std::get<0>(v));
            h ^= std::hash<int32_t>{}(std::get<1>(v)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int32_t>{}(std::get<2>(v)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_map<std::tuple<int32_t,int32_t,int32_t>, uint32_t, Vec3iHash> posMap;
    std::vector<uint32_t> canonicalId(mesh.vertexCount());
    uint32_t nextId = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        auto key = std::make_tuple(qi(mesh.positions[v*3]), qi(mesh.positions[v*3+1]), qi(mesh.positions[v*3+2]));
        auto it = posMap.find(key);
        if (it == posMap.end()) {
            canonicalId[v] = nextId;
            posMap[key] = nextId++;
        } else {
            canonicalId[v] = it->second;
        }
    }
    std::unordered_set<uint64_t> edges;
    size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v0 = canonicalId[mesh.indices[t * 3]];
        uint32_t v1 = canonicalId[mesh.indices[t * 3 + 1]];
        uint32_t v2 = canonicalId[mesh.indices[t * 3 + 2]];
        auto addEdge = [&](uint32_t a, uint32_t b) {
            if (a > b) std::swap(a, b);
            edges.insert((static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b));
        };
        addEdge(v0, v1);
        addEdge(v1, v2);
        addEdge(v2, v0);
    }
    int64_t V = static_cast<int64_t>(nextId);
    int64_t E = static_cast<int64_t>(edges.size());
    int64_t F = static_cast<int64_t>(triCount);
    int64_t chi = V - E + F;
    int64_t g = (2 - chi) / 2;
    return g < 0 ? 0 : static_cast<int>(g);
}

int countNonManifoldEdges(const bromesh::MeshData& mesh) {
    if (mesh.indices.empty()) return 0;
    auto qi = [](float f) -> int32_t { return static_cast<int32_t>(std::round(f * 10000.0f)); };
    struct Vec3iHash {
        size_t operator()(const std::tuple<int32_t,int32_t,int32_t>& v) const {
            size_t h = std::hash<int32_t>{}(std::get<0>(v));
            h ^= std::hash<int32_t>{}(std::get<1>(v)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int32_t>{}(std::get<2>(v)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_map<std::tuple<int32_t,int32_t,int32_t>, uint32_t, Vec3iHash> posMap;
    std::vector<uint32_t> canonicalId(mesh.vertexCount());
    uint32_t nextId = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        auto key = std::make_tuple(qi(mesh.positions[v*3]), qi(mesh.positions[v*3+1]), qi(mesh.positions[v*3+2]));
        auto it = posMap.find(key);
        if (it == posMap.end()) {
            canonicalId[v] = nextId;
            posMap[key] = nextId++;
        } else {
            canonicalId[v] = it->second;
        }
    }
    std::unordered_map<uint64_t, int> edgeCount;
    size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v0 = canonicalId[mesh.indices[t * 3]];
        uint32_t v1 = canonicalId[mesh.indices[t * 3 + 1]];
        uint32_t v2 = canonicalId[mesh.indices[t * 3 + 2]];
        auto addEdge = [&](uint32_t a, uint32_t b) {
            if (a > b) std::swap(a, b);
            edgeCount[(static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b)]++;
        };
        addEdge(v0, v1);
        addEdge(v1, v2);
        addEdge(v2, v0);
    }
    int nonManifold = 0;
    for (const auto& [_, count] : edgeCount) {
        if (count != 2) nonManifold++;
    }
    return nonManifold;
}

} // namespace

void initMeshAnalysis(ObjectBuilder& proto, HostClass& cls) {
    proto.def("surfaceArea", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.surfaceArea: not a Mesh instance");
        return ev::fromDouble(static_cast<double>(bromesh::computeSurfaceArea(m->mesh)));
    });

    proto.def("volume", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.volume: not a Mesh instance");
        return ev::fromDouble(static_cast<double>(bromesh::computeVolume(m->mesh)));
    });

    proto.def("isManifold", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.isManifold: not a Mesh instance");
        return ev::fromBool(bromesh::isManifold(m->mesh));
    });

    proto.def("genus", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.genus: not a Mesh instance");
        return ev::fromDouble(static_cast<double>(computeGenus(m->mesh)));
    });

    proto.def("nonManifoldEdges", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.nonManifoldEdges: not a Mesh instance");
        return ev::fromDouble(static_cast<double>(countNonManifoldEdges(m->mesh)));
    });

    proto.def("curvature", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.curvature: not a Mesh instance");
        float scale = a.empty() ? 1.0f : static_cast<float>(numAt(a, 0));
        bromesh::MeshData copy = m->mesh;
        bromesh::bakeCurvature(copy, scale);
        return makeFloat32Array(copy.colors.data(), copy.colors.size());
    });

    proto.def("computeCurvature", 1, [](Value self, std::span<const Value> a) -> Value {
        return ev::call(ev::getProperty(self, "curvature"), self, a).value;
    });

    proto.def("triangleAreas", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.triangleAreas: not a Mesh instance");
        std::vector<float> areas = bromesh::computeTriangleAreas(m->mesh);
        return makeFloat32Array(areas.data(), areas.size());
    });

    proto.def("sampleSurface", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.sampleSurface: not a Mesh instance");
        ArgReader r(a);
        size_t count = static_cast<size_t>(r.getInt(0, 100));
        uint32_t seed = r.getUint(1, 0);
        return wrapMesh(bromesh::sampleSurface(m->mesh, count, seed));
    });

    // ---- Raycasting & Closest Point -----------------------------------------
    proto.def("raycast", 7, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.raycast: not a Mesh instance");
        float o[3] = {0, 0, 0}, d[3] = {0, -1, 0}, maxDist = 0.0f;
        parseRayArgs(a, o, d, maxDist);
        bromesh::RayHit hit = bromesh::raycast(m->mesh, o, d, maxDist);
        return makeRayHitObject(hit, m->mesh);
    });

    proto.def("raycastAll", 7, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.raycastAll: not a Mesh instance");
        float o[3] = {0, 0, 0}, d[3] = {0, -1, 0}, maxDist = 0.0f;
        parseRayArgs(a, o, d, maxDist);
        auto hits = bromesh::raycastAll(m->mesh, o, d, maxDist);
        return hostArrayOf(hits.size(), [&](size_t i) {
            return makeRayHitObject(hits[i], m->mesh);
        });
    });

    proto.def("raycastTest", 7, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.raycastTest: not a Mesh instance");
        float o[3] = {0, 0, 0}, d[3] = {0, -1, 0}, maxDist = 0.0f;
        parseRayArgs(a, o, d, maxDist);
        return ev::fromBool(bromesh::raycastTest(m->mesh, o, d, maxDist));
    });

    proto.def("closestPoint", 3, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.closestPoint: not a Mesh instance");
        float p[3] = {0, 0, 0};
        parsePointArg(a, p);
        bromesh::RayHit hit = bromesh::closestPoint(m->mesh, p);
        return makeRayHitObject(hit, m->mesh);
    });

    // ---- UV Operations -----------------------------------------------------
    proto.def("unwrapUVs", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.unwrapUVs: not a Mesh instance");
        bromesh::UnwrapParams opts;
        auto res = bromesh::unwrapUVs(m->mesh, opts);
        ObjectBuilder obj;
        obj.set("chartCount", static_cast<double>(res.chartCount));
        obj.set("atlasWidth", static_cast<double>(res.atlasWidth));
        obj.set("atlasHeight", static_cast<double>(res.atlasHeight));
        obj.set("success", ev::fromBool(res.success));
        return obj.build();
    });

    proto.def("projectUVs", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.projectUVs: not a Mesh instance");
        ArgReader r(a);
        int type = r.getInt(0, 0);
        float scale = static_cast<float>(r.getDouble(1, 1.0));
        bromesh::projectUVs(m->mesh, static_cast<bromesh::ProjectionType>(type), scale);
        return self;
    });

    proto.def("computeUVDistortion", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.computeUVDistortion: not a Mesh instance");
        auto d = bromesh::computeUVDistortion(m->mesh);
        std::vector<float> dist;
        dist.reserve(d.size() * 3);
        for (const auto& item : d) {
            dist.push_back(item.stretch);
            dist.push_back(item.areaDistortion);
            dist.push_back(item.angleDistortion);
        }
        return makeFloat32Array(dist.data(), dist.size());
    });

    proto.def("measureUVQuality", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.measureUVQuality: not a Mesh instance");
        auto q = bromesh::measureUVQuality(m->mesh);
        ObjectBuilder obj;
        obj.set("avgStretch", static_cast<double>(q.avgStretch));
        obj.set("maxStretch", static_cast<double>(q.maxStretch));
        obj.set("avgAreaDistortion", static_cast<double>(q.avgAreaDistortion));
        obj.set("maxAreaDistortion", static_cast<double>(q.maxAreaDistortion));
        obj.set("avgAngleDistortion", static_cast<double>(q.avgAngleDistortion));
        obj.set("maxAngleDistortion", static_cast<double>(q.maxAngleDistortion));
        obj.set("uvSpaceUsage", static_cast<double>(q.uvSpaceUsage));
        obj.set("triangleCount", static_cast<double>(q.triangleCount));
        return obj.build();
    });

    // ---- Optimization & Meshlets -------------------------------------------
    proto.def("buildMeshlets", 3, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.buildMeshlets: not a Mesh instance");
        ArgReader r(a);
        size_t maxV = static_cast<size_t>(r.getInt(0, 64));
        size_t maxT = static_cast<size_t>(r.getInt(1, 124));
        float coneW = static_cast<float>(r.getDouble(2, 0.5));
        bromesh::MeshletParams opts;
        opts.maxVertices = maxV;
        opts.maxTriangles = maxT;
        opts.coneWeight = coneW;
        auto meshlets = bromesh::buildMeshlets(m->mesh, opts);
        return hostArrayOf(meshlets.size(), [&](size_t i) {
            ObjectBuilder mo;
            mo.set("vertexCount", static_cast<double>(meshlets[i].vertexCount()));
            mo.set("triangleCount", static_cast<double>(meshlets[i].triangleCount()));
            ev::Persistent vb(makeUint32Array(meshlets[i].vertices.data(), meshlets[i].vertices.size()));
            mo.set("vertices", vb.get());
            ev::Persistent tb(makeUint8Array(meshlets[i].triangles.data(), meshlets[i].triangles.size()));
            mo.set("triangles", tb.get());
            return mo.build();
        });
    });

    proto.def("optimize", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.optimize: not a Mesh instance");
        bromesh::optimizeVertexCache(m->mesh);
        bromesh::optimizeVertexFetch(m->mesh);
        return self;
    });

    proto.def("optimizeVertexCache", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.optimizeVertexCache: not a Mesh instance");
        bromesh::optimizeVertexCache(m->mesh);
        return self;
    });

    proto.def("optimizeVertexFetch", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.optimizeVertexFetch: not a Mesh instance");
        bromesh::optimizeVertexFetch(m->mesh);
        return self;
    });

    proto.def("optimizeOverdraw", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.optimizeOverdraw: not a Mesh instance");
        float threshold = a.empty() ? 1.05f : static_cast<float>(numAt(a, 0));
        bromesh::optimizeOverdraw(m->mesh, threshold);
        return self;
    });

    proto.def("spatialSortTriangles", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.spatialSortTriangles: not a Mesh instance");
        bromesh::spatialSortTriangles(m->mesh);
        return self;
    });

    proto.def("spatialSortVertices", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.spatialSortVertices: not a Mesh instance");
        bromesh::spatialSortVertices(m->mesh);
        return self;
    });

    proto.def("generateShadowIndexBuffer", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.generateShadowIndexBuffer: not a Mesh instance");
        auto s = bromesh::generateShadowIndexBuffer(m->mesh);
        return makeUint32Array(s.data(), s.size());
    });

    // ---- Static Isosurfaces ------------------------------------------------
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        cls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindStatic("marchingCubes", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.marchingCubes: field required");
        std::vector<float> field = toFloatVector(a[0]);
        ArgReader r(a);
        int gx = r.getInt(1, 0), gy = r.getInt(2, 0), gz = r.getInt(3, 0);
        float iso = static_cast<float>(r.getDouble(4, 0.0));
        if (gx <= 0 || gy <= 0 || gz <= 0 || field.size() < static_cast<size_t>(gx * gy * gz)) {
            return ev::throwTypeError("Mesh.marchingCubes: invalid dimensions");
        }
        return wrapMesh(bromesh::marchingCubes(field.data(), gx, gy, gz, iso));
    });

    bindStatic("surfaceNets", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.surfaceNets: field required");
        std::vector<float> field = toFloatVector(a[0]);
        ArgReader r(a);
        int gx = r.getInt(1, 0), gy = r.getInt(2, 0), gz = r.getInt(3, 0);
        float iso = static_cast<float>(r.getDouble(4, 0.0));
        if (gx <= 0 || gy <= 0 || gz <= 0 || field.size() < static_cast<size_t>(gx * gy * gz)) {
            return ev::throwTypeError("Mesh.surfaceNets: invalid dimensions");
        }
        return wrapMesh(bromesh::surfaceNets(field.data(), gx, gy, gz, iso));
    });

    bindStatic("dualContouring", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.dualContouring: field required");
        std::vector<float> field = toFloatVector(a[0]);
        ArgReader r(a);
        int gx = r.getInt(1, 0), gy = r.getInt(2, 0), gz = r.getInt(3, 0);
        float iso = static_cast<float>(r.getDouble(4, 0.0));
        if (gx <= 0 || gy <= 0 || gz <= 0 || field.size() < static_cast<size_t>(gx * gy * gz)) {
            return ev::throwTypeError("Mesh.dualContouring: invalid dimensions");
        }
        return wrapMesh(bromesh::dualContour(field.data(), gx, gy, gz, iso));
    });

    bindStatic("transvoxel", 6, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.transvoxel: field required");
        std::vector<float> field = toFloatVector(a[0]);
        ArgReader r(a);
        int gSize = r.getInt(1, 16);
        int lod = r.getInt(2, 0);
        int nlods[6] = {0, 0, 0, 0, 0, 0};
        float iso = static_cast<float>(r.getDouble(4, 0.0));
        float cell = static_cast<float>(r.getDouble(5, 1.0));
        return wrapMesh(bromesh::transvoxel(field.data(), gSize, lod, nlods, iso, cell));
    });

    bindStatic("greedyMesh", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.greedyMesh: voxels required");
        std::vector<uint8_t> voxels = toUint8Vector(a[0]);
        ArgReader r(a);
        int sx = r.getInt(1, 16), sy = r.getInt(2, 16), sz = r.getInt(3, 16);
        double scale = r.getDouble(4, 1.0);
        if (sx <= 0 || sy <= 0 || sz <= 0 || voxels.size() < static_cast<size_t>(sx * sy * sz)) {
            return ev::throwTypeError("Mesh.greedyMesh: invalid dimensions");
        }
        return wrapMesh(bromesh::greedyMesh(voxels.data(), sx, sy, sz, static_cast<float>(scale)));
    });
}

// ---------------------------------------------------------------------------
// MeshBVH Class
// ---------------------------------------------------------------------------
void initMeshBvh(HostClass& cls) {
    cls.install("MeshBVH", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("MeshBVH constructor: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("MeshBVH: argument must be a Mesh");
        auto h = std::make_unique<HostMeshBVH>();
        h->meshCopy = m->mesh;
        h->bvh = std::make_unique<bromesh::MeshBVH>(bromesh::MeshBVH::build(h->meshCopy));
        return g_meshBvhClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.def("raycast", 7, [](Value self, std::span<const Value> a) -> Value {
            auto* b = unwrapBVH(self);
            if (!b || !b->bvh) return ev::throwTypeError("MeshBVH.raycast: not an instance");
            float o[3] = {0, 0, 0}, d[3] = {0, -1, 0}, maxDist = 0.0f;
            parseRayArgs(a, o, d, maxDist);
            bromesh::RayHit hit = b->bvh->raycast(b->meshCopy, o, d, maxDist);
            return makeRayHitObject(hit, b->meshCopy);
        });

        proto.def("closestPoint", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* b = unwrapBVH(self);
            if (!b || !b->bvh) return ev::throwTypeError("MeshBVH.closestPoint: not an instance");
            float p[3] = {0, 0, 0};
            parsePointArg(a, p);
            bromesh::RayHit hit = b->bvh->closestPoint(b->meshCopy, p);
            return makeRayHitObject(hit, b->meshCopy);
        });
    });
}

// ---------------------------------------------------------------------------
// ProgressiveMesh Class
// ---------------------------------------------------------------------------
void initProgressiveMesh(HostClass& cls) {
    cls.install("ProgressiveMesh", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("ProgressiveMesh constructor: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("ProgressiveMesh: argument must be a Mesh");
        auto h = std::make_unique<HostProgressiveMesh>();
        h->pm = std::make_unique<bromesh::ProgressiveMesh>(bromesh::buildProgressiveMesh(m->mesh));
        return g_progressiveMeshClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("maxTriangles", [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPM(self);
            return ev::fromDouble(p && p->pm ? static_cast<double>(p->pm->maxTriangles()) : 0.0);
        });
        proto.accessor("minTriangles", [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPM(self);
            return ev::fromDouble(p && p->pm ? static_cast<double>(p->pm->minTriangles()) : 0.0);
        });
        proto.def("atRatio", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPM(self);
            if (!p || !p->pm) return ev::throwTypeError("ProgressiveMesh.atRatio: not an instance");
            float r = a.empty() ? 1.0f : static_cast<float>(numAt(a, 0));
            return wrapMesh(bromesh::progressiveMeshAtRatio(*p->pm, r));
        });
        proto.def("atTriangleCount", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPM(self);
            if (!p || !p->pm) return ev::throwTypeError("ProgressiveMesh.atTriangleCount: not an instance");
            int c = a.empty() ? 0 : i32At(a, 0);
            return wrapMesh(bromesh::progressiveMeshAtTriangleCount(*p->pm, c > 0 ? c : 0));
        });
        proto.def("serialize", 0, [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPM(self);
            if (!p || !p->pm) return ev::throwTypeError("ProgressiveMesh.serialize: not an instance");
            std::vector<uint8_t> bytes = bromesh::serializeProgressiveMesh(*p->pm);
            return makeUint8Array(bytes.data(), bytes.size());
        });
    });

    cls.setStatic("deserialize", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("ProgressiveMesh.deserialize: bytes required");
        std::vector<uint8_t> bytes = toUint8Vector(a[0]);
        auto pm = bromesh::deserializeProgressiveMesh(bytes.data(), bytes.size());
        auto h = std::make_unique<HostProgressiveMesh>();
        h->pm = std::make_unique<bromesh::ProgressiveMesh>(std::move(pm));
        return g_progressiveMeshClass.createInstance(std::move(h));
    }, 1, "deserialize"));
}

} // namespace bromesh::api
