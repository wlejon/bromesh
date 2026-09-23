#include "host_mesh_internal.h"
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace bromesh::api {

HostClass g_meshBvhClass;
HostClass g_progressiveMeshClass;

namespace {

// A contour list is an array of flat coordinate arrays; a missing or
// non-array value is no contours.
std::vector<std::vector<float>> readContours(Value v) {
    std::vector<std::vector<float>> contours;
    if (!ev::isObject(v)) return contours;
    Value lenV = ev::getProperty(v, "length");
    if (!ev::isNumber(lenV)) return contours;
    const size_t len = static_cast<size_t>(ev::toDouble(lenV));
    contours.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        contours.push_back(toFloatVector(ev::getElement(v, static_cast<uint32_t>(i))));
    }
    return contours;
}

inline std::string triple(const float* v) {
    return "[" + std::to_string(v[0]) + "," + std::to_string(v[1]) + "," + std::to_string(v[2]) + "]";
}

Value makeRayHitObject(const bromesh::RayHit& h, const bromesh::MeshData& m) {
    if (!h.hit) return ev::null();
    std::string s = "{\"hit\":true,\"distance\":" + std::to_string(h.distance) +
                    ",\"triangle\":" + std::to_string(h.triangleIndex) +
                    ",\"triangleIndex\":" + std::to_string(h.triangleIndex) +
                    ",\"point\":" + triple(h.position) +
                    ",\"position\":" + triple(h.position) +
                    ",\"normal\":" + triple(h.normal) +
                    ",\"barycentric\":[" + std::to_string(h.baryU) + "," + std::to_string(h.baryV) + "," + std::to_string(h.baryW) + "]";
    const size_t t = static_cast<size_t>(h.triangleIndex);
    if (m.hasUVs() && t * 3 + 2 < m.indices.size()) {
        const uint32_t i0 = m.indices[t * 3], i1 = m.indices[t * 3 + 1], i2 = m.indices[t * 3 + 2];
        const float u = h.baryU * m.uvs[i0 * 2] + h.baryV * m.uvs[i1 * 2] + h.baryW * m.uvs[i2 * 2];
        const float v = h.baryU * m.uvs[i0 * 2 + 1] + h.baryV * m.uvs[i1 * 2 + 1] + h.baryW * m.uvs[i2 * 2 + 1];
        s += ",\"uv\":[" + std::to_string(u) + "," + std::to_string(v) + "]";
    }
    s += "}";
    return ev::parseJson(s).value;
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

    // ---- Self-intersection -------------------------------------------------
    // Dropped by the bronze port (bro docs/transition-drift.md H7); the three
    // together are how a caller validates a mesh before a boolean or a
    // physics bake.
    proto.def("hasSelfIntersections", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.hasSelfIntersections: not a Mesh instance");
        return ev::fromBool(bromesh::hasSelfIntersections(m->mesh));
    });

    // findSelfIntersections() -> [{ triA, triB }, ...]
    proto.def("findSelfIntersections", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.findSelfIntersections: not a Mesh instance");
        std::vector<bromesh::TrianglePair> pairs = bromesh::findSelfIntersections(m->mesh);
        return hostArrayOf(pairs.size(), [&pairs](size_t i) -> Value {
            ObjectBuilder o;
            o.set("triA", static_cast<double>(pairs[i].triA));
            o.set("triB", static_cast<double>(pairs[i].triB));
            return o.build();
        });
    });

    // intersectsMesh(other) -> bool — do the two surfaces overlap at all.
    proto.def("intersectsMesh", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.intersectsMesh: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.intersectsMesh: other mesh required");
        auto* other = unwrapMesh(a[0]);
        if (!other) return ev::throwTypeError("Mesh.intersectsMesh: argument must be a Mesh");
        return ev::fromBool(bromesh::meshesIntersect(m->mesh, other->mesh));
    });

    proto.def("buildBVH", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.buildBVH: not a Mesh instance");
        Value bvhCtor = g_meshBvhClass.constructor();
        std::array<Value, 2> args = {self, a.empty() ? ev::fromDouble(8.0) : a[0]};
        return ev::call(bvhCtor, ev::undefined(), args).value;
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
        int type = 0;
        if (!a.empty()) {
            if (ev::isString(a[0])) {
                std::string s = ev::toUtf8(a[0]);
                if (s == "box") type = 0;
                else if (s == "planarXY") type = 1;
                else if (s == "planarXZ") type = 2;
                else if (s == "planarYZ") type = 3;
                else if (s == "cylindrical") type = 4;
                else if (s == "spherical") type = 5;
                else return ev::throwTypeError(("Mesh.projectUVs: unknown projection '" + s + "'").c_str());
            } else if (ev::isNumber(a[0])) {
                type = static_cast<int>(ev::toDouble(a[0]));
            }
        }
        float scale = a.size() > 1 && ev::isNumber(a[1]) ? static_cast<float>(ev::toDouble(a[1])) : 1.0f;
        bromesh::projectUVs(m->mesh, static_cast<bromesh::ProjectionType>(type), scale);
        return self;
    });

    proto.def("generateUVs", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.generateUVs: not a Mesh instance");
        std::string method = "unwrap";
        if (!a.empty() && ev::isString(a[0])) method = ev::toUtf8(a[0]);
        if (method == "unwrap" || method == "xatlas") {
            return ev::call(ev::getProperty(self, "unwrapUVs"), self, {}).value;
        }
        return ev::call(ev::getProperty(self, "projectUVs"), self, std::array<Value, 2>{a[0], ev::fromDouble(1.0)}).value;
    });

    proto.def("computeUVDistortion", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.computeUVDistortion: not a Mesh instance");
        auto d = bromesh::computeUVDistortion(m->mesh);
        return hostArrayOf(d.size(), [&](size_t i) -> Value {
            ObjectBuilder item;
            item.set("stretch", static_cast<double>(d[i].stretch));
            item.set("areaDistortion", static_cast<double>(d[i].areaDistortion));
            item.set("angleDistortion", static_cast<double>(d[i].angleDistortion));
            return item.build();
        });
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
        size_t maxV = 64;
        size_t maxT = 124;
        float coneW = 0.5f;
        if (!a.empty() && ev::isObject(a[0])) {
            Value opts = a[0];
            Value mv = ev::getProperty(opts, "maxVertices");
            Value mt = ev::getProperty(opts, "maxTriangles");
            Value cw = ev::getProperty(opts, "coneWeight");
            if (!ev::isUndefined(mv) && !ev::isNull(mv)) maxV = static_cast<size_t>(ev::toDouble(mv));
            if (!ev::isUndefined(mt) && !ev::isNull(mt)) maxT = static_cast<size_t>(ev::toDouble(mt));
            if (!ev::isUndefined(cw) && !ev::isNull(cw)) coneW = static_cast<float>(ev::toDouble(cw));
        } else {
            ArgReader r(a);
            maxV = static_cast<size_t>(r.getInt(0, 64));
            maxT = static_cast<size_t>(r.getInt(1, 124));
            coneW = static_cast<float>(r.getDouble(2, 0.5));
        }
        bromesh::MeshletParams opts;
        opts.maxVertices = maxV > 0 ? maxV : 64;
        opts.maxTriangles = maxT > 0 ? maxT : 124;
        opts.coneWeight = coneW >= 0.0f ? coneW : 0.5f;
        auto meshlets = bromesh::buildMeshlets(m->mesh, opts);

        std::vector<uint32_t> allVertices;
        std::vector<uint8_t> allTriangles;
        std::vector<uint8_t> records;
        records.reserve(meshlets.size() * 64);

        auto putU32 = [&](uint32_t v) {
            uint8_t b[4];
            std::memcpy(b, &v, 4);
            records.insert(records.end(), b, b + 4);
        };
        auto putF32 = [&](float v) {
            uint8_t b[4];
            std::memcpy(b, &v, 4);
            records.insert(records.end(), b, b + 4);
        };

        for (const auto& item : meshlets) {
            putU32(static_cast<uint32_t>(allVertices.size()));
            putU32(static_cast<uint32_t>(item.vertices.size()));
            putU32(static_cast<uint32_t>(allTriangles.size()));
            putU32(static_cast<uint32_t>(item.triangles.size() / 3));
            for (float f : item.bounds.center) putF32(f);
            putF32(item.bounds.radius);
            for (float f : item.bounds.coneApex) putF32(f);
            for (float f : item.bounds.coneAxis) putF32(f);
            putF32(item.bounds.coneCutoff);
            putU32(0);

            allVertices.insert(allVertices.end(), item.vertices.begin(), item.vertices.end());
            allTriangles.insert(allTriangles.end(), item.triangles.begin(), item.triangles.end());
        }

        Value arr = hostArrayOf(meshlets.size(), [&](size_t i) -> Value {
            ObjectBuilder mo;
            mo.set("vertexCount", static_cast<double>(meshlets[i].vertices.size()));
            mo.set("triangleCount", static_cast<double>(meshlets[i].triangles.size() / 3));
            mo.set("vertices", makeUint32Array(meshlets[i].vertices.data(), meshlets[i].vertices.size()));
            mo.set("triangles", makeUint8Array(meshlets[i].triangles.data(), meshlets[i].triangles.size()));
            // center / coneApex / coneAxis are what a cluster-cull shader
            // actually needs; radius+cutoff alone cannot place the cone.
            const auto& bb = meshlets[i].bounds;
            ObjectBuilder bnd;
            {
                ev::Persistent c(hostArrayOf(3, [&bb](size_t k) {
                    return ev::fromDouble(static_cast<double>(bb.center[k]));
                }));
                bnd.set("center", c.get());
            }
            bnd.set("radius", static_cast<double>(bb.radius));
            {
                ev::Persistent ap(hostArrayOf(3, [&bb](size_t k) {
                    return ev::fromDouble(static_cast<double>(bb.coneApex[k]));
                }));
                bnd.set("coneApex", ap.get());
            }
            {
                ev::Persistent ax(hostArrayOf(3, [&bb](size_t k) {
                    return ev::fromDouble(static_cast<double>(bb.coneAxis[k]));
                }));
                bnd.set("coneAxis", ax.get());
            }
            bnd.set("coneCutoff", static_cast<double>(bb.coneCutoff));
            mo.set("bounds", bnd.build());
            return mo.build();
        });

        ev::setProperty(arr, "meshletCount", ev::fromDouble(static_cast<double>(meshlets.size())));
        ev::Persistent vBuf(makeUint32Array(allVertices.data(), allVertices.size()));
        ev::setProperty(arr, "vertices", vBuf.get());
        ev::Persistent tBuf(makeUint8Array(allTriangles.data(), allTriangles.size()));
        ev::setProperty(arr, "triangles", tBuf.get());
        ev::Persistent rBuf(makeUint8Array(records.data(), records.size()));
        ev::setProperty(arr, "meshlets", ev::getProperty(rBuf.get(), "buffer"));
        return arr;
    });


    proto.def("optimize", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.optimize: not a Mesh instance");
        bromesh::optimizeVertexCache(m->mesh);
        bromesh::optimizeVertexFetch(m->mesh);
        return self;
    });

    proto.def("analyzeVertexCache", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.analyzeVertexCache: not a Mesh instance");
        unsigned int cacheSize = !a.empty() && ev::isNumber(a[0]) ? static_cast<unsigned int>(ev::toDouble(a[0])) : 16u;
        auto st = bromesh::analyzeVertexCache(m->mesh, cacheSize > 0 ? cacheSize : 16u);
        ObjectBuilder out;
        out.set("verticesTransformed", static_cast<double>(st.verticesTransformed));
        out.set("warpsExecuted", static_cast<double>(st.warpsExecuted));
        out.set("acmr", static_cast<double>(st.acmr));
        out.set("atvr", static_cast<double>(st.atvr));
        return out.build();
    });

    proto.def("analyzeVertexFetch", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.analyzeVertexFetch: not a Mesh instance");
        size_t vertexSize = !a.empty() && ev::isNumber(a[0]) ? static_cast<size_t>(ev::toDouble(a[0])) : 32u;
        auto st = bromesh::analyzeVertexFetch(m->mesh, vertexSize > 0 ? vertexSize : 32u);
        ObjectBuilder out;
        out.set("bytesFetched", static_cast<double>(st.bytesFetched));
        out.set("overfetch", static_cast<double>(st.overfetch));
        return out.build();
    });

    proto.def("analyzeOverdraw", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.analyzeOverdraw: not a Mesh instance");
        auto st = bromesh::analyzeOverdraw(m->mesh);
        ObjectBuilder out;
        out.set("pixelsCovered", static_cast<double>(st.pixelsCovered));
        out.set("pixelsShaded", static_cast<double>(st.pixelsShaded));
        out.set("overdraw", static_cast<double>(st.overdraw));
        return out.build();
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

    proto.def("encode", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.encode: not a Mesh instance");
        bromesh::EncodedMesh enc = bromesh::encodeMesh(m->mesh);
        ObjectBuilder out;
        out.set("vertexData", makeUint8Array(enc.vertexData.data(), enc.vertexData.size()));
        out.set("indexData", makeUint8Array(enc.indexData.data(), enc.indexData.size()));
        out.set("vertexCount", static_cast<double>(enc.vertexCount));
        out.set("vertexSize", static_cast<double>(enc.vertexSize));
        out.set("indexCount", static_cast<double>(enc.indexCount));
        out.set("hasNormals", ev::fromBool(!m->mesh.normals.empty()));
        out.set("hasUVs", ev::fromBool(!m->mesh.uvs.empty()));
        out.set("hasColors", ev::fromBool(!m->mesh.colors.empty()));
        return out.build();
    });

    // ---- Static Isosurfaces & Compression ----------------------------------
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        cls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindStatic("stripify", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.stripify(indices, vertexCount): required arguments");
        std::vector<uint32_t> indices = toUint32Vector(a[0]);
        size_t vertexCount = static_cast<size_t>(ev::toDouble(a[1]));
        uint32_t restartIndex = a.size() > 2 && ev::isNumber(a[2]) ? static_cast<uint32_t>(ev::toDouble(a[2])) : 0xFFFFFFFF;
        auto strip = bromesh::stripify(indices, vertexCount, restartIndex);
        return makeUint32Array(strip.data(), strip.size());
    });

    bindStatic("unstripify", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.unstripify(strip): strip required");
        std::vector<uint32_t> strip = toUint32Vector(a[0]);
        uint32_t restartIndex = a.size() > 1 && ev::isNumber(a[1]) ? static_cast<uint32_t>(ev::toDouble(a[1])) : 0xFFFFFFFF;
        auto indices = bromesh::unstripify(strip, restartIndex);
        return makeUint32Array(indices.data(), indices.size());
    });

    bindStatic("decode", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isObject(a[0])) return ev::throwTypeError("Mesh.decode: encoded object required");
        Value encVal = a[0];
        bromesh::EncodedMesh enc;
        enc.vertexData = toUint8Vector(ev::getProperty(encVal, "vertexData"));
        enc.indexData = toUint8Vector(ev::getProperty(encVal, "indexData"));
        Value vc = ev::getProperty(encVal, "vertexCount");
        Value vs = ev::getProperty(encVal, "vertexSize");
        Value ic = ev::getProperty(encVal, "indexCount");
        enc.vertexCount = ev::isNumber(vc) ? static_cast<size_t>(ev::toDouble(vc)) : 0;
        enc.vertexSize = ev::isNumber(vs) ? static_cast<size_t>(ev::toDouble(vs)) : 0;
        enc.indexCount = ev::isNumber(ic) ? static_cast<size_t>(ev::toDouble(ic)) : 0;

        bool hasNormals = true;
        bool hasUVs = true;
        bool hasColors = false;
        Value hn = ev::getProperty(encVal, "hasNormals");
        Value hu = ev::getProperty(encVal, "hasUVs");
        Value hc = ev::getProperty(encVal, "hasColors");
        if (ev::isBool(hn)) hasNormals = ev::toBool(hn);
        if (ev::isBool(hu)) hasUVs = ev::toBool(hu);
        if (ev::isBool(hc)) hasColors = ev::toBool(hc);

        bromesh::MeshData mesh = bromesh::decodeMesh(enc, hasNormals, hasUVs, hasColors);
        return wrapMesh(std::move(mesh));
    });

    // ---- Polygon triangulation and point-cloud reconstruction -----------
    bindStatic("polygon2D", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.polygon2D: (outer[, holes[, z]]) required");
        std::vector<float> outer = toFloatVector(a[0]);
        std::vector<std::vector<float>> holes = readContours(argAt(a, 1));
        const float z = static_cast<float>(numAt(a, 2));
        return wrapMesh(bromesh::triangulatePolygon2D(outer, holes, z));
    });
    bindStatic("polygon3D", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.polygon3D: (outer, holes, normal) required");
        std::vector<float> outer = toFloatVector(a[0]);
        std::vector<std::vector<float>> holes = readContours(argAt(a, 1));
        float normal[3] = {0.0f, 0.0f, 1.0f};
        std::vector<float> n = toFloatVector(argAt(a, 2));
        if (n.size() >= 3) { normal[0] = n[0]; normal[1] = n[1]; normal[2] = n[2]; }
        return wrapMesh(bromesh::triangulatePolygon3D(outer, holes, normal));
    });
    bindStatic("reconstruct", 2, [](Value, std::span<const Value> a) -> Value {
        HostMesh* cloud = a.empty() ? nullptr : unwrapMesh(a[0]);
        if (!cloud) return ev::throwTypeError("Mesh.reconstruct: a Mesh point cloud (positions + normals) required");
        bromesh::ReconstructParams params;
        if (a.size() > 1 && ev::isObject(a[1])) {
            Value gr = ev::getProperty(a[1], "gridResolution");
            if (ev::isNumber(gr)) params.gridResolution = static_cast<int>(ev::toDouble(gr));
            Value sr = ev::getProperty(a[1], "supportRadius");
            if (ev::isNumber(sr)) params.supportRadius = static_cast<float>(ev::toDouble(sr));
            Value il = ev::getProperty(a[1], "isoLevel");
            if (ev::isNumber(il)) params.isoLevel = static_cast<float>(ev::toDouble(il));
        }
        return wrapMesh(bromesh::reconstructFromPointCloud(cloud->mesh, params));
    });

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
    bindStatic("dualContour", 5, [](Value, std::span<const Value> a) -> Value {
        return ev::call(ev::getProperty(g_meshClass.constructor(), "dualContouring"), ev::undefined(), a).value;
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
    cls.install("MeshBVH", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("MeshBVH constructor: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) {
            auto* b = unwrapBVH(a[0]);
            if (b) return ev::throwTypeError("expected a __bro_native.mesh.Mesh handle, got a __bro_native.mesh.MeshBVH handle");
            return ev::throwTypeError("expected a __bro_native.mesh.Mesh handle, got an ordinary object");
        }
        int leafSize = a.size() > 1 && ev::isNumber(a[1]) ? static_cast<int>(ev::toDouble(a[1])) : 8;
        auto h = std::make_unique<HostMeshBVH>();
        h->meshCopy = m->mesh;
        h->bvh = std::make_unique<bromesh::MeshBVH>(bromesh::MeshBVH::build(h->meshCopy, leafSize > 0 ? leafSize : 8));
        return g_meshBvhClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("empty", [](Value self, std::span<const Value>) -> Value {
            auto* b = unwrapBVH(self);
            return ev::fromBool(!b || !b->bvh || b->bvh->empty());
        });
        proto.accessor("triangleCount", [](Value self, std::span<const Value>) -> Value {
            auto* b = unwrapBVH(self);
            return ev::fromDouble(b && b->bvh ? static_cast<double>(b->bvh->triangleCount()) : 0.0);
        });
        proto.accessor("nodeCount", [](Value self, std::span<const Value>) -> Value {
            auto* b = unwrapBVH(self);
            return ev::fromDouble(b && b->bvh ? static_cast<double>(b->bvh->nodeCount()) : 0.0);
        });
        proto.def("bounds", 0, [](Value self, std::span<const Value>) -> Value {
            auto* b = unwrapBVH(self);
            if (!b || !b->bvh) return ev::throwTypeError("MeshBVH.bounds: not an instance");
            const bromath::AABB3 bb = b->bvh->bounds();
            ObjectBuilder out;
            out.set("minX", static_cast<double>(bb.min.x));
            out.set("minY", static_cast<double>(bb.min.y));
            out.set("minZ", static_cast<double>(bb.min.z));
            out.set("maxX", static_cast<double>(bb.max.x));
            out.set("maxY", static_cast<double>(bb.max.y));
            out.set("maxZ", static_cast<double>(bb.max.z));
            out.set("centerX", static_cast<double>((bb.min.x + bb.max.x) * 0.5f));
            out.set("centerY", static_cast<double>((bb.min.y + bb.max.y) * 0.5f));
            out.set("centerZ", static_cast<double>((bb.min.z + bb.max.z) * 0.5f));
            out.set("extentX", static_cast<double>(bb.max.x - bb.min.x));
            out.set("extentY", static_cast<double>(bb.max.y - bb.min.y));
            out.set("extentZ", static_cast<double>(bb.max.z - bb.min.z));
            out.set("min", hostArrayOf(3, [&](size_t i) {
                return ev::fromDouble(i == 0 ? bb.min.x : (i == 1 ? bb.min.y : bb.min.z));
            }));
            out.set("max", hostArrayOf(3, [&](size_t i) {
                return ev::fromDouble(i == 0 ? bb.max.x : (i == 1 ? bb.max.y : bb.max.z));
            }));
            return out.build();
        });
        proto.def("raycast", 7, [](Value self, std::span<const Value> a) -> Value {
            auto* b = unwrapBVH(self);
            if (!b || !b->bvh) return ev::throwTypeError("MeshBVH.raycast: not an instance");
            const bromesh::MeshData* targetMesh = &b->meshCopy;
            std::span<const Value> rayArgs = a;
            if (!a.empty()) {
                auto* m = unwrapMesh(a[0]);
                if (m) {
                    targetMesh = &m->mesh;
                    rayArgs = a.subspan(1);
                }
            }
            float o[3] = {0, 0, 0}, d[3] = {0, -1, 0}, maxDist = 0.0f;
            parseRayArgs(rayArgs, o, d, maxDist);
            bromesh::RayHit hit = b->bvh->raycast(*targetMesh, o, d, maxDist);
            return makeRayHitObject(hit, *targetMesh);
        });
        proto.def("raycastTest", 7, [](Value self, std::span<const Value> a) -> Value {
            auto* b = unwrapBVH(self);
            if (!b || !b->bvh) return ev::throwTypeError("MeshBVH.raycastTest: not an instance");
            const bromesh::MeshData* targetMesh = &b->meshCopy;
            std::span<const Value> rayArgs = a;
            if (!a.empty()) {
                auto* m = unwrapMesh(a[0]);
                if (m) {
                    targetMesh = &m->mesh;
                    rayArgs = a.subspan(1);
                }
            }
            float o[3] = {0, 0, 0}, d[3] = {0, -1, 0}, maxDist = 0.0f;
            parseRayArgs(rayArgs, o, d, maxDist);
            return ev::fromBool(b->bvh->raycastTest(*targetMesh, o, d, maxDist));
        });
        proto.def("closestPoint", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* b = unwrapBVH(self);
            if (!b || !b->bvh) return ev::throwTypeError("MeshBVH.closestPoint: not an instance");
            const bromesh::MeshData* targetMesh = &b->meshCopy;
            std::span<const Value> ptArgs = a;
            if (!a.empty()) {
                auto* m = unwrapMesh(a[0]);
                if (m) {
                    targetMesh = &m->mesh;
                    ptArgs = a.subspan(1);
                }
            }
            float p[3] = {0, 0, 0};
            parsePointArg(ptArgs, p);
            bromesh::RayHit hit = b->bvh->closestPoint(*targetMesh, p);
            return makeRayHitObject(hit, *targetMesh);
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
        proto.accessor("collapseCount", [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPM(self);
            return ev::fromDouble(p && p->pm ? static_cast<double>(p->pm->collapses.size()) : 0.0);
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
