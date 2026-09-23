#include "host_mesh_internal.h"
#if BROMESH_HAS_DRACO
#include <bromesh/io/draco.h>
#endif
#include <vector>
#include <utility>

namespace bromesh::api {

namespace {

Value makeTextureResult(const bromesh::TextureBuffer& tb) {
    ObjectBuilder res;
    res.set("width", static_cast<double>(tb.width));
    res.set("height", static_cast<double>(tb.height));
    res.set("channels", static_cast<double>(tb.channels));
    {
        ev::Persistent px(makeFloat32Array(tb.pixels.data(), tb.pixels.size()));
        res.set("pixels", px.get());
        res.set("data", px.get());
    }
    return res.build();
}

// A bake texture size (w at args[i], h at args[i + 1]); each axis is
// [1, kMaxAxis], and the bakers size their buffer as an int product.
bool texSizeArgs(std::span<const Value> a, size_t i, const char* fn, int& w, int& h) {
    const std::string name(fn);
    return countArg(a, i, name + ": width", 1, kMaxAxis, w) &&
           countArg(a, i + 1, name + ": height", 1, kMaxAxis, h);
}

// A hemisphere ray count: at least one (the bakers divide by it).
bool raysArg(std::span<const Value> a, size_t i, const char* fn, int& rays) {
    return countArg(a, i, std::string(fn) + ": rays", 1, kMaxIterations, rays);
}

} // namespace

void initMeshOps(ObjectBuilder& proto, HostClass& cls) {
    // ---- Normal Operations -------------------------------------------------
    // ---- Normal Operations -------------------------------------------------
    proto.def("computeNormals", 1, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) {
            auto* b = unwrapBVH(self);
            if (b) return ev::throwTypeError("expected a __bro_native.mesh.Mesh handle, got a __bro_native.mesh.MeshBVH handle");
            return ev::throwTypeError("Mesh.computeNormals: not a Mesh instance");
        }
        bromesh::computeNormals(m->mesh);
        return self;
    });


    proto.def("computeFlatNormals", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.computeFlatNormals: not a Mesh instance");
        return wrapMesh(bromesh::computeFlatNormals(m->mesh));
    });

    // Per-vertex tangents as a flat Float32Array (xyzw per vertex, w = the
    // bitangent sign). Needed for normal-mapped materials; dropped by the
    // bronze port (bro docs/transition-drift.md H7).
    proto.def("computeTangents", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.computeTangents: not a Mesh instance");
        std::vector<float> tangents = bromesh::computeTangents(m->mesh);
        return makeFloat32Array(tangents.data(), tangents.size());
    });

    proto.def("computeCreaseNormals", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.computeCreaseNormals: not a Mesh instance");
        double angle = a.empty() ? 60.0 : numAt(a, 0);
        m->mesh = bromesh::computeCreaseNormals(m->mesh, static_cast<float>(angle));
        return self;
    });

    proto.def("invertNormals", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.invertNormals: not a Mesh instance");
        for (float& val : m->mesh.normals) val = -val;
        return self;
    });

    proto.def("flipFaces", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.flipFaces: not a Mesh instance");
        for (size_t i = 0; i + 2 < m->mesh.indices.size(); i += 3) {
            std::swap(m->mesh.indices[i + 1], m->mesh.indices[i + 2]);
        }
        return self;
    });

    // ---- Clean / Repair / Weld ---------------------------------------------
    proto.def("weld", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.weld: not a Mesh instance");
        double eps = a.empty() ? 1e-4 : numAt(a, 0);
        m->mesh = bromesh::weldVertices(m->mesh, static_cast<float>(eps));
        return self;
    });

    proto.def("removeDegenerateTriangles", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.removeDegenerateTriangles: not a Mesh instance");
        double eps = a.empty() ? 1e-6 : numAt(a, 0);
        m->mesh = bromesh::removeDegenerateTriangles(m->mesh, static_cast<float>(eps));
        return self;
    });

    proto.def("removeDuplicateTriangles", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.removeDuplicateTriangles: not a Mesh instance");
        m->mesh = bromesh::removeDuplicateTriangles(m->mesh);
        return self;
    });

    proto.def("fillHoles", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.fillHoles: not a Mesh instance");
        int maxEdges = 32;
        if (!countArg(a, 0, "Mesh.fillHoles: maxEdges", 0, kMaxInt32, maxEdges)) return ev::undefined();
        m->mesh = bromesh::fillHoles(m->mesh, maxEdges);
        return self;
    });

    proto.def("repair", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.repair: not a Mesh instance");
        m->mesh = bromesh::removeDegenerateTriangles(m->mesh);
        m->mesh = bromesh::removeDuplicateTriangles(m->mesh);
        return self;
    });

    // ---- Simplification / LOD ----------------------------------------------
    proto.def("simplify", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.simplify: not a Mesh instance");
        ArgReader r(a);
        float ratio = static_cast<float>(r.getDouble(0, 0.5));
        float targetError = static_cast<float>(r.getDouble(1, 1e-3));
        m->mesh = bromesh::simplify(m->mesh, ratio, targetError);
        return self;
    });

    // simplifyWithAttributes(ratio, error=0.01, uvWeight=1, normalWeight=0.5) —
    // quadric simplification that folds UV and normal error into the metric, so
    // UV seams and hard edges survive. The plain simplify() does not.
    proto.def("simplifyWithAttributes", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.simplifyWithAttributes: not a Mesh instance");
        ArgReader r(a);
        float ratio = static_cast<float>(r.getDouble(0, 0.5));
        float targetError = static_cast<float>(r.getDouble(1, 0.01));
        float uvWeight = static_cast<float>(r.getDouble(2, 1.0));
        float normalWeight = static_cast<float>(r.getDouble(3, 0.5));
        m->mesh = bromesh::simplifyWithAttributes(m->mesh, ratio, targetError, uvWeight, normalWeight);
        return self;
    });

    proto.def("simplifyToTriangleCount", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.simplifyToTriangleCount: not a Mesh instance");
        ArgReader r(a);
        size_t count = 100;
        if (!countArg(a, 0, "Mesh.simplifyToTriangleCount: count", 0, kMaxUint32, count)) return ev::undefined();
        float targetError = static_cast<float>(r.getDouble(1, 1e-3));
        m->mesh = bromesh::simplifyToTriangleCount(m->mesh, count, targetError);
        return self;
    });

    proto.def("generateLODChain", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.generateLODChain: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.generateLODChain: ratios array required");
        std::vector<float> ratios = toFloatVector(a[0]);
        auto chain = bromesh::generateLODChain(m->mesh, ratios.data(), static_cast<int>(ratios.size()));
        return hostArrayOf(chain.size(), [&](size_t i) {
            return wrapMesh(std::move(chain[i]));
        });
    });

    // ---- Subdivision & Smoothing -------------------------------------------
    proto.def("subdivideLoop", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.subdivideLoop: not a Mesh instance");
        int iters = 1;  // 0 is a no-op
        if (!countArg(a, 0, "Mesh.subdivideLoop: iterations", 0, kMaxSubdivisions, iters)) return ev::undefined();
        m->mesh = bromesh::subdivideLoop(m->mesh, iters);
        return self;
    });

    proto.def("subdivideCatmullClark", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.subdivideCatmullClark: not a Mesh instance");
        int iters = 1;  // 0 is a no-op
        if (!countArg(a, 0, "Mesh.subdivideCatmullClark: iterations", 0, kMaxSubdivisions, iters)) {
            return ev::undefined();
        }
        m->mesh = bromesh::subdivideCatmullClark(m->mesh, iters);
        return self;
    });

    // Plain 1-to-4 midpoint split: no smoothing, so the surface is unchanged
    // and only density rises. The subdivision a displacement pass wants.
    proto.def("subdivideMidpoint", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.subdivideMidpoint: not a Mesh instance");
        int iters = 1;  // 0 is a no-op
        if (!countArg(a, 0, "Mesh.subdivideMidpoint: iterations", 0, kMaxSubdivisions, iters)) return ev::undefined();
        m->mesh = bromesh::subdivideMidpoint(m->mesh, iters);
        return self;
    });

    proto.def("smooth", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.smooth: not a Mesh instance");
        ArgReader r(a);
        float lambda = static_cast<float>(r.getDouble(0, 0.5));
        int iters = 1;  // 0 is a no-op
        if (!countArg(a, 1, "Mesh.smooth: iterations", 0, kMaxIterations, iters)) return ev::undefined();
        bromesh::smoothLaplacian(m->mesh, lambda, iters);
        return self;
    });

    proto.def("smoothLaplacian", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.smoothLaplacian: not a Mesh instance");
        ArgReader r(a);
        float lambda = static_cast<float>(r.getDouble(0, 0.5));
        int iters = 1;
        if (!countArg(a, 1, "Mesh.smoothLaplacian: iterations", 0, kMaxIterations, iters)) return ev::undefined();
        bromesh::smoothLaplacian(m->mesh, lambda, iters);
        return self;
    });

    proto.def("smoothTaubin", 3, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.smoothTaubin: not a Mesh instance");
        ArgReader r(a);
        float lambda = static_cast<float>(r.getDouble(0, 0.5));
        float mu = static_cast<float>(r.getDouble(1, -0.53));
        int iters = 1;
        if (!countArg(a, 2, "Mesh.smoothTaubin: iterations", 0, kMaxIterations, iters)) return ev::undefined();
        bromesh::smoothTaubin(m->mesh, lambda, mu, iters);
        return self;
    });

    proto.def("remesh", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.remesh: not a Mesh instance");
        double len = a.empty() ? 0.1 : numAt(a, 0);
        m->mesh = bromesh::remeshIsotropic(m->mesh, static_cast<float>(len), 3);
        return self;
    });

    proto.def("remeshIsotropic", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.remeshIsotropic: not a Mesh instance");
        ArgReader r(a);
        float len = static_cast<float>(r.getDouble(0, 0.1));
        int iters = 3;
        if (!countArg(a, 1, "Mesh.remeshIsotropic: iterations", 0, kMaxIterations, iters)) return ev::undefined();
        m->mesh = bromesh::remeshIsotropic(m->mesh, len, iters);
        return self;
    });

    proto.def("shrinkwrap", 5, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.shrinkwrap: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.shrinkwrap: target mesh required");
        auto* target = unwrapMesh(a[0]);
        if (!target) return ev::throwTypeError("Mesh.shrinkwrap: target must be a Mesh");
        ev::Persistent selfP(self);  // the axis read below may allocate
        int mode = 0;
        if (a.size() > 1) {
            if (ev::isString(a[1])) {
                std::string sm = ev::toUtf8(a[1]);
                if (sm == "normal" || sm == "projectAlongNormal") mode = 1;
                else if (sm == "axis" || sm == "projectAlongAxis") mode = 2;
                else mode = 0;
            } else if (!countArg(a, 1, "Mesh.shrinkwrap: mode", 0, 2, mode)) {
                return ev::undefined();
            }
        }
        ArgReader r(a);
        float maxDist = static_cast<float>(r.getDouble(2, 0.0));
        float offset = static_cast<float>(r.getDouble(3, 0.0));
        bromath::Vec3 axis{0, 1, 0};
        const float* axisPtr = nullptr;
        if (a.size() > 4 && !ev::isNull(a[4]) && !ev::isUndefined(a[4])) {
            std::vector<float> v = toFloatVector(a[4]);
            if (v.size() >= 3) {
                axis = {v[0], v[1], v[2]};
                axisPtr = &axis.x;
            }
        }
        bromesh::shrinkwrap(m->mesh, target->mesh, static_cast<bromesh::ShrinkwrapMode>(mode), maxDist, offset, axisPtr);
        return selfP.get();
    });


    // ---- Splitting / CSG / Decomposition -----------------------------------
    proto.def("splitByPlane", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.splitByPlane: not a Mesh instance");
        ArgReader r(a);
        float nx = static_cast<float>(r.getDouble(0, 0.0));
        float ny = static_cast<float>(r.getDouble(1, 1.0));
        float nz = static_cast<float>(r.getDouble(2, 0.0));
        float d = static_cast<float>(r.getDouble(3, 0.0));
        auto pair = bromesh::splitByPlane(m->mesh, nx, ny, nz, d);
        return hostArrayOf(2, [&](size_t i) {
            return wrapMesh(i == 0 ? std::move(pair.first) : std::move(pair.second));
        });
    });

    proto.def("splitComponents", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.splitComponents: not a Mesh instance");
        auto comps = bromesh::splitConnectedComponents(m->mesh);
        return hostArrayOf(comps.size(), [&](size_t i) {
            return wrapMesh(std::move(comps[i]));
        });
    });

    proto.def("convexDecomposition", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.convexDecomposition: not a Mesh instance");
        // Two call forms: the pre-transition options object
        // ({maxHulls, maxVerticesPerHull, resolution, minVolumePerHull}) and
        // the bronze-era positional one. The object form was the documented
        // shape and the only one the old binding read.
        int maxHulls = 16;
        int maxVerts = 32;
        double res = 100000.0;
        double minVol = 0.0001;
        if (!a.empty() && ev::isObject(a[0])) {
            ev::Persistent o(a[0]);
            if (!countField(o.get(), "maxHulls", "Mesh.convexDecomposition: maxHulls", 1, kMaxAxis, maxHulls) ||
                !countField(o.get(), "maxVerticesPerHull", "Mesh.convexDecomposition: maxVerticesPerHull", 4,
                            kMaxIterations, maxVerts)) {
                return ev::undefined();
            }
            Value v = ev::getProperty(o.get(), "resolution");
            if (ev::isNumber(v)) res = ev::toDouble(v);
            v = ev::getProperty(o.get(), "minVolumePerHull");
            if (ev::isNumber(v)) minVol = ev::toDouble(v);
        } else {
            ArgReader r(a);
            if (!countArg(a, 0, "Mesh.convexDecomposition: maxHulls", 1, kMaxAxis, maxHulls) ||
                !countArg(a, 1, "Mesh.convexDecomposition: maxVerticesPerHull", 4, kMaxIterations, maxVerts)) {
                return ev::undefined();
            }
            res = r.getDouble(2, 100000.0);
            minVol = r.getDouble(3, 0.0001);
        }
        bromesh::ConvexDecompParams opts;
        opts.maxHulls = maxHulls;
        opts.maxVerticesPerHull = maxVerts;
        opts.resolution = static_cast<float>(res > 0 ? res : 100000.0f);
        opts.minVolumePerHull = static_cast<float>(minVol);
        auto hulls = bromesh::convexDecomposition(m->mesh, opts);
        return hostArrayOf(hulls.size(), [&](size_t i) {
            return wrapMesh(std::move(hulls[i]));
        });
    });

    proto.def("convexHull", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.convexHull: not a Mesh instance");
        return wrapMesh(bromesh::convexHull(m->mesh));
    });

    proto.def("booleanUnion", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* aMesh = unwrapMesh(self);
        if (!aMesh) return ev::throwTypeError("Mesh.booleanUnion: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.booleanUnion: other mesh required");
        auto* bMesh = unwrapMesh(a[0]);
        if (!bMesh) return ev::throwTypeError("Mesh.booleanUnion: argument must be a Mesh");
        return wrapMesh(bromesh::booleanUnion(aMesh->mesh, bMesh->mesh));
    });

    proto.def("booleanDifference", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* aMesh = unwrapMesh(self);
        if (!aMesh) return ev::throwTypeError("Mesh.booleanDifference: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.booleanDifference: other mesh required");
        auto* bMesh = unwrapMesh(a[0]);
        if (!bMesh) return ev::throwTypeError("Mesh.booleanDifference: argument must be a Mesh");
        return wrapMesh(bromesh::booleanDifference(aMesh->mesh, bMesh->mesh));
    });

    proto.def("booleanIntersection", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* aMesh = unwrapMesh(self);
        if (!aMesh) return ev::throwTypeError("Mesh.booleanIntersection: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.booleanIntersection: other mesh required");
        auto* bMesh = unwrapMesh(a[0]);
        if (!bMesh) return ev::throwTypeError("Mesh.booleanIntersection: argument must be a Mesh");
        return wrapMesh(bromesh::booleanIntersection(aMesh->mesh, bMesh->mesh));
    });

    // Aliases
    proto.def("union", 1, [](Value self, std::span<const Value> a) -> Value {
        return callMethod(self, "booleanUnion", a);
    });
    proto.def("subtract", 1, [](Value self, std::span<const Value> a) -> Value {
        return callMethod(self, "booleanDifference", a);
    });
    proto.def("intersect", 1, [](Value self, std::span<const Value> a) -> Value {
        return callMethod(self, "booleanIntersection", a);
    });

    // ---- Baking ------------------------------------------------------------
    proto.def("bakeAmbientOcclusion", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeAmbientOcclusion: not a Mesh instance");
        ArgReader r(a);
        int rays = 64;
        if (!raysArg(a, 0, "Mesh.bakeAmbientOcclusion", rays)) return ev::undefined();
        float dist = static_cast<float>(r.getDouble(1, 0.0));
        bromesh::bakeAmbientOcclusion(m->mesh, rays, dist);
        return self;
    });

    proto.def("bakeCurvature", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeCurvature: not a Mesh instance");
        float scale = static_cast<float>(numAt(a, 0));
        if (scale <= 0.0f) scale = 1.0f;
        bromesh::bakeCurvature(m->mesh, scale);
        return self;
    });

    proto.def("bakeThickness", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeThickness: not a Mesh instance");
        ArgReader r(a);
        int rays = 32;
        if (!raysArg(a, 0, "Mesh.bakeThickness", rays)) return ev::undefined();
        float dist = static_cast<float>(r.getDouble(1, 0.0));
        bromesh::bakeThickness(m->mesh, rays, dist);
        return self;
    });

    proto.def("bakeAOToTexture", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeAOToTexture: not a Mesh instance");
        ArgReader r(a);
        int w = 512, h = 512, rays = 64;
        if (!texSizeArgs(a, 0, "Mesh.bakeAOToTexture", w, h) || !raysArg(a, 2, "Mesh.bakeAOToTexture", rays)) {
            return ev::undefined();
        }
        float dist = static_cast<float>(r.getDouble(3, 0.0));
        auto tb = bromesh::bakeAmbientOcclusionToTexture(m->mesh, w, h, rays, dist);
        return makeTextureResult(tb);
    });

    proto.def("bakeCurvatureToTexture", 3, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeCurvatureToTexture: not a Mesh instance");
        ArgReader r(a);
        int w = 512, h = 512;
        if (!texSizeArgs(a, 0, "Mesh.bakeCurvatureToTexture", w, h)) return ev::undefined();
        float scale = static_cast<float>(r.getDouble(2, 1.0));
        auto tb = bromesh::bakeCurvatureToTexture(m->mesh, w, h, scale);
        return makeTextureResult(tb);
    });

    proto.def("bakeThicknessToTexture", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeThicknessToTexture: not a Mesh instance");
        ArgReader r(a);
        int w = 512, h = 512, rays = 32;
        if (!texSizeArgs(a, 0, "Mesh.bakeThicknessToTexture", w, h) ||
            !raysArg(a, 2, "Mesh.bakeThicknessToTexture", rays)) {
            return ev::undefined();
        }
        float dist = static_cast<float>(r.getDouble(3, 0.0));
        auto tb = bromesh::bakeThicknessToTexture(m->mesh, w, h, rays, dist);
        return makeTextureResult(tb);
    });

    proto.def("bakeNormalsToTexture", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakeNormalsToTexture: not a Mesh instance");
        int w = 512, h = 512;
        if (!texSizeArgs(a, 0, "Mesh.bakeNormalsToTexture", w, h)) return ev::undefined();
        auto tb = bromesh::bakeNormalsToTexture(m->mesh, w, h);
        return makeTextureResult(tb);
    });

    proto.def("bakePositionToTexture", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bakePositionToTexture: not a Mesh instance");
        int w = 512, h = 512;
        if (!texSizeArgs(a, 0, "Mesh.bakePositionToTexture", w, h)) return ev::undefined();
        auto tb = bromesh::bakePositionToTexture(m->mesh, w, h);
        return makeTextureResult(tb);
    });

    proto.def("bakeNormalsFromReference", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* low = unwrapMesh(self);
        if (!low) return ev::throwTypeError("Mesh.bakeNormalsFromReference: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.bakeNormalsFromReference: high mesh required");
        auto* high = unwrapMesh(a[0]);
        if (!high) return ev::throwTypeError("Mesh.bakeNormalsFromReference: high must be a Mesh");
        ArgReader r(a);
        int w = 512, h = 512;
        if (!texSizeArgs(a, 1, "Mesh.bakeNormalsFromReference", w, h)) return ev::undefined();
        float dist = static_cast<float>(r.getDouble(3, 0.0));
        auto tb = bromesh::bakeNormalsFromReference(low->mesh, high->mesh, w, h, dist);
        return makeTextureResult(tb);
    });

    proto.def("bakeAOFromReference", 5, [](Value self, std::span<const Value> a) -> Value {
        auto* low = unwrapMesh(self);
        if (!low) return ev::throwTypeError("Mesh.bakeAOFromReference: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.bakeAOFromReference: high mesh required");
        auto* high = unwrapMesh(a[0]);
        if (!high) return ev::throwTypeError("Mesh.bakeAOFromReference: high must be a Mesh");
        ArgReader r(a);
        int w = 512, h = 512, rays = 64;
        if (!texSizeArgs(a, 1, "Mesh.bakeAOFromReference", w, h) ||
            !raysArg(a, 3, "Mesh.bakeAOFromReference", rays)) {
            return ev::undefined();
        }
        float dist = static_cast<float>(r.getDouble(4, 0.0));
        auto tb = bromesh::bakeAOFromReference(low->mesh, high->mesh, w, h, rays, dist);
        return makeTextureResult(tb);
    });

    // ---- Static helpers ----------------------------------------------------
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        cls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindStatic("merge", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapMesh(bromesh::MeshData{});
        std::vector<bromesh::MeshData> meshes;
        if (ev::isObject(a[0])) {
            Value lenVal = ev::getProperty(a[0], "length");
            if (ev::isNumber(lenVal)) {
                size_t n = lengthValue(lenVal);
                for (size_t i = 0; i < n; ++i) {
                    Value elem = ev::getElement(a[0], static_cast<uint32_t>(i));
                    auto* m = unwrapMesh(elem);
                    if (!m) {
                        return ev::throwTypeError("expected a __bro_native.mesh.Mesh handle, got a " + std::string(ev::isObject(elem) ? "wrong object" : "non-object"));
                    }
                    meshes.push_back(m->mesh);
                }
            } else {
                for (size_t i = 0; i < a.size(); ++i) {
                    auto* m = unwrapMesh(a[i]);
                    if (!m) {
                        return ev::throwTypeError("expected a __bro_native.mesh.Mesh handle, got a " + std::string(ev::isObject(a[i]) ? "wrong object" : "non-object"));
                    }
                    meshes.push_back(m->mesh);
                }
            }
        }
        return wrapMesh(bromesh::mergeMeshes(meshes));
    });

    bindStatic("splitByPlane", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.splitByPlane: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Mesh.splitByPlane: not a Mesh instance");
        ArgReader r(a.subspan(1));
        float nx = static_cast<float>(r.getDouble(0, 0.0));
        float ny = static_cast<float>(r.getDouble(1, 1.0));
        float nz = static_cast<float>(r.getDouble(2, 0.0));
        float d = static_cast<float>(r.getDouble(3, 0.0));
        auto pair = bromesh::splitByPlane(m->mesh, nx, ny, nz, d);
        return hostArrayOf(2, [&](size_t i) {
            return wrapMesh(i == 0 ? std::move(pair.first) : std::move(pair.second));
        });
    });


    bindStatic("booleanUnion", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.booleanUnion(a, b): two meshes required");
        auto* ma = unwrapMesh(a[0]);
        auto* mb = unwrapMesh(a[1]);
        if (!ma || !mb) return ev::throwTypeError("Mesh.booleanUnion: arguments must be Meshes");
        return wrapMesh(bromesh::booleanUnion(ma->mesh, mb->mesh));
    });
    bindStatic("union", 2, [](Value, std::span<const Value> a) -> Value {
        return callMethod(g_meshClass.constructor(), "booleanUnion", a);
    });

    bindStatic("booleanDifference", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.booleanDifference(a, b): two meshes required");
        auto* ma = unwrapMesh(a[0]);
        auto* mb = unwrapMesh(a[1]);
        if (!ma || !mb) return ev::throwTypeError("Mesh.booleanDifference: arguments must be Meshes");
        return wrapMesh(bromesh::booleanDifference(ma->mesh, mb->mesh));
    });
    bindStatic("subtract", 2, [](Value, std::span<const Value> a) -> Value {
        return callMethod(g_meshClass.constructor(), "booleanDifference", a);
    });

    bindStatic("booleanIntersection", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.booleanIntersection(a, b): two meshes required");
        auto* ma = unwrapMesh(a[0]);
        auto* mb = unwrapMesh(a[1]);
        if (!ma || !mb) return ev::throwTypeError("Mesh.booleanIntersection: arguments must be Meshes");
        return wrapMesh(bromesh::booleanIntersection(ma->mesh, mb->mesh));
    });
    bindStatic("intersect", 2, [](Value, std::span<const Value> a) -> Value {
        return callMethod(g_meshClass.constructor(), "booleanIntersection", a);
    });

    bindStatic("convexHull", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.convexHull: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Mesh.convexHull: argument must be a Mesh");
        return wrapMesh(bromesh::convexHull(m->mesh));
    });

#if BROMESH_HAS_DRACO
    bindStatic("decodeDraco", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.decodeDraco: bytes required");
        const uint8_t* data = nullptr;
        size_t size = 0;
        if (!readUint8Array(a[0], data, size)) {
            ev::TypedArrayInfo info = ev::typedArrayInfo(a[0]);
            if (info && info.data) {
                data = reinterpret_cast<const uint8_t*>(info.data);
                size = info.byteLength;
            } else {
                return ev::throwTypeError("Mesh.decodeDraco: expected Uint8Array or TypedArray view");
            }
        }
        bromesh::DracoDecoded dec = bromesh::decodeDraco(data, size);
        if (!dec.ok()) {
            return ev::throwError("Mesh.decodeDraco: " + (dec.error.empty() ? "decode failed" : dec.error));
        }
        ObjectBuilder res;
        res.set("positions", makeFloat32Array(dec.mesh.positions.data(), dec.mesh.positions.size()));
        res.set("indices", makeUint32Array(dec.mesh.indices.data(), dec.mesh.indices.size()));
        res.set("normals", makeFloat32Array(dec.mesh.normals.data(), dec.mesh.normals.size()));
        res.set("uvs", makeFloat32Array(dec.mesh.uvs.data(), dec.mesh.uvs.size()));
        res.set("colors", makeFloat32Array(dec.mesh.colors.data(), dec.mesh.colors.size()));
        res.set("mesh", wrapMesh(std::move(dec.mesh)));

        Value attrs = hostArrayOf(dec.attributes.size(), [&](size_t i) -> Value {
            const auto& attr = dec.attributes[i];
            ObjectBuilder ab;
            ab.set("type", ev::fromUtf8(attr.type.c_str()));
            ab.set("uniqueId", ev::fromDouble(attr.uniqueId));
            ab.set("components", ev::fromDouble(attr.components));
            ab.set("count", ev::fromDouble(attr.count));
            const char* kindName = "float32";
            switch (attr.kind) {
                case bromesh::DracoAttribute::Kind::Float32: kindName = "float32"; break;
                case bromesh::DracoAttribute::Kind::Int8:    kindName = "int8"; break;
                case bromesh::DracoAttribute::Kind::Uint8:   kindName = "uint8"; break;
                case bromesh::DracoAttribute::Kind::Int16:   kindName = "int16"; break;
                case bromesh::DracoAttribute::Kind::Uint16:  kindName = "uint16"; break;
                case bromesh::DracoAttribute::Kind::Int32:   kindName = "int32"; break;
                case bromesh::DracoAttribute::Kind::Uint32:  kindName = "uint32"; break;
            }
            ab.set("kind", ev::fromUtf8(kindName));
            ab.set("bytes", makeUint8Array(attr.bytes.data(), attr.bytes.size()));
            return ab.build();
        });
        res.set("attributes", attrs);
        return res.build();
    });

    bindStatic("encodeDraco", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.encodeDraco: mesh or {positions, indices} required");
        bromesh::MeshData mesh;
        auto* m = unwrapMesh(a[0]);
        if (m) {
            mesh = m->mesh;
        } else if (ev::isObject(a[0])) {
            // Each field is read right before it is converted: a read (and a
            // conversion) may allocate, which would stale an earlier read.
            mesh.positions = toFloatVector(ev::getProperty(a[0], "positions"));
            mesh.indices = toUint32Vector(ev::getProperty(a[0], "indices"));
            if (Value v = ev::getProperty(a[0], "normals"); !ev::isUndefined(v)) mesh.normals = toFloatVector(v);
            if (Value v = ev::getProperty(a[0], "uvs"); !ev::isUndefined(v)) mesh.uvs = toFloatVector(v);
            if (Value v = ev::getProperty(a[0], "colors"); !ev::isUndefined(v)) mesh.colors = toFloatVector(v);
        } else {
            return ev::throwTypeError("Mesh.encodeDraco: first argument must be a Mesh or object");
        }

        bromesh::DracoEncodeOptions opts;
        if (a.size() > 1 && ev::isObject(a[1])) {
            // Quantization is 1..30 bits and the speed dial 0..10, Draco's
            // own ranges; `speed` wins over the `compressionLevel` spelling.
            if (!countField(a[1], "positionBits", "Mesh.encodeDraco: positionBits", 1, 30, opts.positionBits) ||
                !countField(a[1], "normalBits", "Mesh.encodeDraco: normalBits", 1, 30, opts.normalBits) ||
                !countField(a[1], "uvBits", "Mesh.encodeDraco: uvBits", 1, 30, opts.uvBits) ||
                !countField(a[1], "colorBits", "Mesh.encodeDraco: colorBits", 1, 30, opts.colorBits) ||
                !countField(a[1], "compressionLevel", "Mesh.encodeDraco: compressionLevel", 0, 10, opts.speed) ||
                !countField(a[1], "speed", "Mesh.encodeDraco: speed", 0, 10, opts.speed)) {
                return ev::undefined();
            }
        }

        std::string err;
        std::vector<uint8_t> bytes = bromesh::encodeDraco(mesh, opts, &err);
        if (bytes.empty() && !err.empty()) {
            return ev::throwError("Mesh.encodeDraco: " + err);
        }
        return makeUint8Array(bytes.data(), bytes.size());
    });
#endif
}

} // namespace bromesh::api
