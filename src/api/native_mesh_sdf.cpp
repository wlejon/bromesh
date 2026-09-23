#include "host_mesh_internal.h"
#include "arg_reader.h"
#include "bromesh/isosurface/jit/sdf_node.h"
#include "bromesh/isosurface/jit/sdf_compiler.h"
#include "bromesh/isosurface/jit/jit_mesher.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace bromesh::api {

HostClass g_sdfGraphClass;

static bromath::Vec3 readVec3Param(Value v) {
    bromath::Vec3 r{0.0f, 0.0f, 0.0f};
    if (!ev::isObject(v)) return r;
    // Rooted: each read may allocate. Numbers are immediates, so a component
    // read is converted at once and nothing but the root is held across.
    ev::Persistent obj(v);
    auto comp = [](Value c, float& out) {
        if (ev::isNumber(c)) out = static_cast<float>(ev::toDouble(c));
    };
    Value x = ev::getProperty(obj.get(), "x");
    if (ev::isNumber(x)) {
        r.x = static_cast<float>(ev::toDouble(x));
        comp(ev::getProperty(obj.get(), "y"), r.y);
        comp(ev::getProperty(obj.get(), "z"), r.z);
        return r;
    }
    comp(ev::getElement(obj.get(), 0), r.x);
    comp(ev::getElement(obj.get(), 1), r.y);
    comp(ev::getElement(obj.get(), 2), r.z);
    return r;
}

static bromath::Vec3 parseVec3(std::span<const Value> a, size_t& idx) {
    if (idx < a.size() && ev::isObject(a[idx])) {
        return readVec3Param(a[idx++]);
    }
    float x = static_cast<float>(numAt(a, idx++));
    float y = static_cast<float>(numAt(a, idx++));
    float z = static_cast<float>(numAt(a, idx++));
    return {x, y, z};
}

static bromath::AABB3 parseBounds(Value v, bromath::AABB3 def = {{-2.0f, -2.0f, -2.0f}, {2.0f, 2.0f, 2.0f}}) {
    if (ev::isObject(v)) {
        ev::Persistent obj(v);  // rooted across the allocating reads
        ev::Persistent minVal(ev::getProperty(obj.get(), "min"));
        ev::Persistent maxVal(ev::getProperty(obj.get(), "max"));
        if (ev::isObject(minVal.get()) && ev::isObject(maxVal.get())) {
            const bromath::Vec3 lo = readVec3Param(minVal.get());
            const bromath::Vec3 hi = readVec3Param(maxVal.get());
            return {lo, hi};
        }
        Value sizeVal = ev::getProperty(obj.get(), "size");
        if (ev::isNumber(sizeVal)) {
            float s = static_cast<float>(ev::toDouble(sizeVal));
            return {{-s, -s, -s}, {s, s, s}};
        }
    }
    return def;
}

struct SdfMeshOptions {
    int dimX = 64;
    int dimY = 64;
    int dimZ = 64;
    bromath::AABB3 bounds = {{-2.0f, -2.0f, -2.0f}, {2.0f, 2.0f, 2.0f}};
    float isoLevel = 0.0f;
    bool closeBoundary = true;
    bool computeGradients = true;
};

static SdfMeshOptions parseMeshOptions(Value optsVal) {
    SdfMeshOptions opt;
    if (!ev::isObject(optsVal)) {
        return opt;
    }
    ev::Persistent opts(optsVal);  // rooted: every read below may allocate

    {
        ev::Persistent dims(ev::getProperty(opts.get(), "dims"));
        if (ev::isObject(dims.get())) {
            int* out[3] = {&opt.dimX, &opt.dimY, &opt.dimZ};
            for (uint32_t i = 0; i < 3; ++i) {
                Value d = ev::getElement(dims.get(), i);
                if (ev::isNumber(d)) *out[i] = static_cast<int>(ev::toDouble(d));
            }
        } else if (ev::isNumber(dims.get())) {
            int d = static_cast<int>(ev::toDouble(dims.get()));
            opt.dimX = opt.dimY = opt.dimZ = d;
        }
    }

    Value resVal = ev::getProperty(opts.get(), "resolution");
    if (ev::isNumber(resVal)) {
        int d = static_cast<int>(ev::toDouble(resVal));
        opt.dimX = opt.dimY = opt.dimZ = d;
    }

    Value dx = ev::getProperty(opts.get(), "dimX");
    if (ev::isNumber(dx)) opt.dimX = static_cast<int>(ev::toDouble(dx));
    Value dy = ev::getProperty(opts.get(), "dimY");
    if (ev::isNumber(dy)) opt.dimY = static_cast<int>(ev::toDouble(dy));
    Value dz = ev::getProperty(opts.get(), "dimZ");
    if (ev::isNumber(dz)) opt.dimZ = static_cast<int>(ev::toDouble(dz));

    Value boundsVal = ev::getProperty(opts.get(), "bounds");
    if (!ev::isUndefined(boundsVal)) {
        opt.bounds = parseBounds(boundsVal, opt.bounds);
    }

    Value isoVal = ev::getProperty(opts.get(), "isoLevel");
    if (ev::isNumber(isoVal)) opt.isoLevel = static_cast<float>(ev::toDouble(isoVal));

    Value cbVal = ev::getProperty(opts.get(), "closeBoundary");
    if (ev::isBool(cbVal)) opt.closeBoundary = ev::toBool(cbVal);

    Value cgVal = ev::getProperty(opts.get(), "computeGradients");
    if (ev::isBool(cgVal)) opt.computeGradients = ev::toBool(cgVal);

    return opt;
}

void initSdfGraph(HostClass& cls) {
    cls.install("SDFGraph", 0, [](Value, std::span<const Value>) -> Value {
        return g_sdfGraphClass.createInstance(std::make_unique<HostSdfGraph>());
    }, [](ObjectBuilder& proto) {
        // Primitives
        proto.def("sphere", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.sphere: invalid instance");
            float r = a.empty() ? 1.0f : static_cast<float>(numAt(a, 0));
            return ev::fromDouble(h->graph.sphere(r));
        });

        proto.def("box", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.box: invalid instance");
            size_t idx = 0;
            bromath::Vec3 ext = parseVec3(a, idx);
            return ev::fromDouble(h->graph.box(ext));
        });

        proto.def("roundedBox", 4, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.roundedBox: invalid instance");
            size_t idx = 0;
            bromath::Vec3 ext = parseVec3(a, idx);
            float r = static_cast<float>(numAt(a, idx));
            return ev::fromDouble(h->graph.roundedBox(ext, r));
        });

        proto.def("cylinder", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.cylinder: invalid instance");
            float r = static_cast<float>(numAt(a, 0));
            float hh = static_cast<float>(numAt(a, 1));
            return ev::fromDouble(h->graph.cylinder(r, hh));
        });

        proto.def("capsule", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.capsule: invalid instance");
            size_t idx = 0;
            bromath::Vec3 p0 = parseVec3(a, idx);
            bromath::Vec3 p1 = parseVec3(a, idx);
            float r = static_cast<float>(numAt(a, idx));
            return ev::fromDouble(h->graph.capsule(p0, p1, r));
        });

        proto.def("torus", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.torus: invalid instance");
            float majorR = static_cast<float>(numAt(a, 0));
            float minorR = static_cast<float>(numAt(a, 1));
            return ev::fromDouble(h->graph.torus(majorR, minorR));
        });

        proto.def("plane", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.plane: invalid instance");
            size_t idx = 0;
            bromath::Vec3 n = parseVec3(a, idx);
            float d = static_cast<float>(numAt(a, idx));
            return ev::fromDouble(h->graph.plane(n, d));
        });

        // CSG Combinators
        auto bindCsg = [&](const char* name, auto fn) {
            proto.def(name, 2, [fn](Value self, std::span<const Value> a) -> Value {
                auto* h = unwrapSdfGraph(self);
                if (!h) return ev::throwTypeError("SDFGraph method: invalid instance");
                int n1 = i32At(a, 0);
                int n2 = i32At(a, 1);
                return ev::fromDouble((h->graph.*fn)(n1, n2));
            });
        };
        bindCsg("opUnion", &SdfGraph::opUnion);
        bindCsg("union", &SdfGraph::opUnion);
        bindCsg("opIntersection", &SdfGraph::opIntersection);
        bindCsg("intersection", &SdfGraph::opIntersection);
        bindCsg("opSubtraction", &SdfGraph::opSubtraction);
        bindCsg("subtraction", &SdfGraph::opSubtraction);

        auto bindSmoothCsg = [&](const char* name, auto fn) {
            proto.def(name, 3, [fn](Value self, std::span<const Value> a) -> Value {
                auto* h = unwrapSdfGraph(self);
                if (!h) return ev::throwTypeError("SDFGraph method: invalid instance");
                int n1 = i32At(a, 0);
                int n2 = i32At(a, 1);
                float k = static_cast<float>(numAt(a, 2));
                return ev::fromDouble((h->graph.*fn)(n1, n2, k));
            });
        };
        bindSmoothCsg("opSmoothUnion", &SdfGraph::opSmoothUnion);
        bindSmoothCsg("smoothUnion", &SdfGraph::opSmoothUnion);
        bindSmoothCsg("opSmoothIntersection", &SdfGraph::opSmoothIntersection);
        bindSmoothCsg("smoothIntersection", &SdfGraph::opSmoothIntersection);
        bindSmoothCsg("opSmoothSubtraction", &SdfGraph::opSmoothSubtraction);
        bindSmoothCsg("smoothSubtraction", &SdfGraph::opSmoothSubtraction);

        // Transforms & Modifiers
        proto.def("translate", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.translate: invalid instance");
            int child = i32At(a, 0);
            size_t idx = 1;
            bromath::Vec3 offset = parseVec3(a, idx);
            return ev::fromDouble(h->graph.translate(child, offset));
        });

        proto.def("rotateY", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.rotateY: invalid instance");
            int child = i32At(a, 0);
            float angle = static_cast<float>(numAt(a, 1));
            return ev::fromDouble(h->graph.rotateY(child, angle));
        });

        proto.def("scaleUniform", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.scaleUniform: invalid instance");
            int child = i32At(a, 0);
            float s = static_cast<float>(numAt(a, 1));
            return ev::fromDouble(h->graph.scaleUniform(child, s));
        });

        proto.def("noise3D", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.noise3D: invalid instance");
            size_t idx = 0;
            bromath::Vec3 freq = parseVec3(a, idx);
            float amp = static_cast<float>(numAt(a, idx));
            return ev::fromDouble(h->graph.noise3D(freq, amp));
        });

        proto.def("displace", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.displace: invalid instance");
            int child = i32At(a, 0);
            int noise = i32At(a, 1);
            return ev::fromDouble(h->graph.displace(child, noise));
        });

        proto.def("setRoot", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.setRoot: invalid instance");
            h->graph.setRoot(i32At(a, 0));
            return self;
        });

        proto.accessor("root", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapSdfGraph(self);
            return ev::fromDouble(h ? h->graph.root() : -1);
        });

        proto.accessor("size", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapSdfGraph(self);
            return ev::fromDouble(h ? static_cast<double>(h->graph.size()) : 0.0);
        });

        proto.def("computeHash", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.computeHash: invalid instance");
            return ev::fromDouble(static_cast<double>(h->graph.computeHash()));
        });

        // Meshing methods on SDFGraph instance
        proto.def("marchingCubes", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.marchingCubes: invalid instance");
            SdfMeshOptions opt = parseMeshOptions(!a.empty() ? a[0] : ev::undefined());
            MeshData m = marchingCubesFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                              opt.bounds, opt.isoLevel,
                                              opt.closeBoundary, opt.computeGradients);
            return wrapMesh(std::move(m));
        });

        proto.def("surfaceNets", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.surfaceNets: invalid instance");
            SdfMeshOptions opt = parseMeshOptions(!a.empty() ? a[0] : ev::undefined());
            MeshData m = surfaceNetsFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                           opt.bounds, opt.isoLevel);
            return wrapMesh(std::move(m));
        });
    });
}

void initMeshSdf(ObjectBuilder&, HostClass& meshCls) {
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        meshCls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindStatic("createSDF", 0, [](Value, std::span<const Value>) -> Value {
        return g_sdfGraphClass.createInstance(std::make_unique<HostSdfGraph>());
    });

    auto extractGraphAndOptions = [](std::span<const Value> a, SdfGraph& outGraph, SdfMeshOptions& outOpt) -> bool {
        if (a.empty()) return false;
        if (auto* h = unwrapSdfGraph(a[0])) {
            outGraph = h->graph;
            outOpt = parseMeshOptions(a.size() > 1 ? a[1] : ev::undefined());
            return true;
        }
        if (ev::isObject(a[0])) {
            Value gVal = ev::getProperty(a[0], "graph");
            if (auto* h = unwrapSdfGraph(gVal)) {
                outGraph = h->graph;
                outOpt = parseMeshOptions(a[0]);
                return true;
            }
            // Declarative shorthand: e.g. { type: "sphere", radius: 1.0 }
            Value typeVal = ev::getProperty(a[0], "type");
            if (!ev::isString(typeVal)) typeVal = ev::getProperty(a[0], "op");
            if (ev::isString(typeVal)) {
                std::string t = ev::toUtf8(typeVal);
                if (t == "sphere") {
                    Value rVal = ev::getProperty(a[0], "radius");
                    float r = ev::isNumber(rVal) ? static_cast<float>(ev::toDouble(rVal)) : 1.0f;
                    outGraph.sphere(r);
                } else if (t == "box") {
                    Value extVal = ev::getProperty(a[0], "halfExtents");
                    bromath::Vec3 ext = ev::isObject(extVal) ? readVec3Param(extVal) : bromath::Vec3{1.0f, 1.0f, 1.0f};
                    outGraph.box(ext);
                } else if (t == "torus") {
                    float majR = static_cast<float>(ev::toDouble(ev::getProperty(a[0], "majorRadius")));
                    float minR = static_cast<float>(ev::toDouble(ev::getProperty(a[0], "minorRadius")));
                    outGraph.torus(majR > 0 ? majR : 1.0f, minR > 0 ? minR : 0.3f);
                } else if (t == "cylinder") {
                    float r = static_cast<float>(ev::toDouble(ev::getProperty(a[0], "radius")));
                    float hh = static_cast<float>(ev::toDouble(ev::getProperty(a[0], "halfHeight")));
                    outGraph.cylinder(r > 0 ? r : 0.5f, hh > 0 ? hh : 1.0f);
                }
                outOpt = parseMeshOptions(a[0]);
                return true;
            }
        }
        return false;
    };

    bindStatic("marchingCubesSDF", 2, [extractGraphAndOptions](Value, std::span<const Value> a) -> Value {
        SdfGraph graph;
        SdfMeshOptions opt;
        if (!extractGraphAndOptions(a, graph, opt)) {
            return ev::throwTypeError("Mesh.marchingCubesSDF: valid SDFGraph or options required");
        }
        MeshData m = marchingCubesFromSDF(graph, opt.dimX, opt.dimY, opt.dimZ,
                                          opt.bounds, opt.isoLevel,
                                          opt.closeBoundary, opt.computeGradients);
        return wrapMesh(std::move(m));
    });

    bindStatic("surfaceNetsSDF", 2, [extractGraphAndOptions](Value, std::span<const Value> a) -> Value {
        SdfGraph graph;
        SdfMeshOptions opt;
        if (!extractGraphAndOptions(a, graph, opt)) {
            return ev::throwTypeError("Mesh.surfaceNetsSDF: valid SDFGraph or options required");
        }
        MeshData m = surfaceNetsFromSDF(graph, opt.dimX, opt.dimY, opt.dimZ,
                                       opt.bounds, opt.isoLevel);
        return wrapMesh(std::move(m));
    });
}

} // namespace bromesh::api
