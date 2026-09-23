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

// A node id argument: an existing node of `g`, [0, size). A node can only
// name nodes built before it, which is also what keeps the graph acyclic —
// a combinator given its own (not yet allocated) id would recurse forever
// when evaluated.
static bool nodeArg(const SdfGraph& g, std::span<const Value> a, size_t i, const char* what, int& out) {
    if (g.size() == 0) {
        ev::throwRangeError(std::string(what) + ": the graph has no nodes yet");
        return false;
    }
    if (i >= a.size()) {
        ev::throwTypeError(std::string(what) + " must be a number");
        return false;
    }
    int64_t v = 0;
    if (!intValue(a[i], what, 0, static_cast<double>(g.size() - 1), v)) return false;
    out = static_cast<int>(v);
    return true;
}

// Returns false once it has thrown (a bad grid dimension).
static bool parseMeshOptions(Value optsVal, SdfMeshOptions& opt) {
    if (!ev::isObject(optsVal)) {
        return true;
    }
    ev::Persistent opts(optsVal);  // rooted: every read below may allocate

    // Each axis is [2, kMaxVolumeAxis] samples (a grid needs two to hold a
    // cell) and the grid at most kMaxVolumeCells.
    {
        ev::Persistent dims(ev::getProperty(opts.get(), "dims"));
        if (ev::isObject(dims.get())) {
            int* out[3] = {&opt.dimX, &opt.dimY, &opt.dimZ};
            for (uint32_t i = 0; i < 3; ++i) {
                Value d = ev::getElement(dims.get(), i);
                if (ev::isUndefined(d)) continue;
                int64_t n = 0;
                if (!intValue(d, "SDF mesh options: dims", 2, kMaxVolumeAxis, n)) return false;
                *out[i] = static_cast<int>(n);
            }
        } else if (!ev::isUndefined(dims.get())) {
            int64_t n = 0;
            if (!intValue(dims.get(), "SDF mesh options: dims", 2, kMaxVolumeAxis, n)) return false;
            opt.dimX = opt.dimY = opt.dimZ = static_cast<int>(n);
        }
    }

    {
        int res = 0;
        if (!countField(opts.get(), "resolution", "SDF mesh options: resolution", 2, kMaxVolumeAxis, res)) {
            return false;
        }
        if (res > 0) opt.dimX = opt.dimY = opt.dimZ = res;
    }

    if (!countField(opts.get(), "dimX", "SDF mesh options: dimX", 2, kMaxVolumeAxis, opt.dimX) ||
        !countField(opts.get(), "dimY", "SDF mesh options: dimY", 2, kMaxVolumeAxis, opt.dimY) ||
        !countField(opts.get(), "dimZ", "SDF mesh options: dimZ", 2, kMaxVolumeAxis, opt.dimZ) ||
        !volumeCellsOk("SDF mesh options", opt.dimX, opt.dimY, opt.dimZ)) {
        return false;
    }

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

    return true;
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
                int n1 = 0, n2 = 0;
                if (!nodeArg(h->graph, a, 0, "SDFGraph: node a", n1) ||
                    !nodeArg(h->graph, a, 1, "SDFGraph: node b", n2)) {
                    return ev::undefined();
                }
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
                int n1 = 0, n2 = 0;
                if (!nodeArg(h->graph, a, 0, "SDFGraph: node a", n1) ||
                    !nodeArg(h->graph, a, 1, "SDFGraph: node b", n2)) {
                    return ev::undefined();
                }
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
            int child = 0;
            if (!nodeArg(h->graph, a, 0, "SDFGraph.translate: child", child)) return ev::undefined();
            size_t idx = 1;
            bromath::Vec3 offset = parseVec3(a, idx);
            return ev::fromDouble(h->graph.translate(child, offset));
        });

        proto.def("rotateY", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.rotateY: invalid instance");
            int child = 0;
            if (!nodeArg(h->graph, a, 0, "SDFGraph.rotateY: child", child)) return ev::undefined();
            float angle = static_cast<float>(numAt(a, 1));
            return ev::fromDouble(h->graph.rotateY(child, angle));
        });

        proto.def("scaleUniform", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.scaleUniform: invalid instance");
            int child = 0;
            if (!nodeArg(h->graph, a, 0, "SDFGraph.scaleUniform: child", child)) return ev::undefined();
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
            int child = 0, noise = 0;
            if (!nodeArg(h->graph, a, 0, "SDFGraph.displace: child", child) ||
                !nodeArg(h->graph, a, 1, "SDFGraph.displace: noise", noise)) {
                return ev::undefined();
            }
            return ev::fromDouble(h->graph.displace(child, noise));
        });

        proto.def("setRoot", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.setRoot: invalid instance");
            int root = 0;
            if (!nodeArg(h->graph, a, 0, "SDFGraph.setRoot: node", root)) return ev::undefined();
            h->graph.setRoot(root);
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
            SdfMeshOptions opt;
            if (!parseMeshOptions(!a.empty() ? a[0] : ev::undefined(), opt)) return ev::undefined();
            MeshData m = marchingCubesFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                              opt.bounds, opt.isoLevel,
                                              opt.closeBoundary, opt.computeGradients);
            return wrapMesh(std::move(m));
        });

        proto.def("surfaceNets", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.surfaceNets: invalid instance");
            SdfMeshOptions opt;
            if (!parseMeshOptions(!a.empty() ? a[0] : ev::undefined(), opt)) return ev::undefined();
            MeshData m = surfaceNetsFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                           opt.bounds, opt.isoLevel);
            return wrapMesh(std::move(m));
        });
    });
}

void initMeshSdf(ObjectBuilder&, HostClass& meshCls) {
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        meshCls.setStatic(name, hostFunction(std::move(fn), arity, name));
    };

    bindStatic("createSDF", 0, [](Value, std::span<const Value>) -> Value {
        return g_sdfGraphClass.createInstance(std::make_unique<HostSdfGraph>());
    });

    // NoGraph: nothing usable was passed. Threw: an option was bad and the
    // exception is pending.
    enum class Extract { Ok, NoGraph, Threw };
    auto extractGraphAndOptions = [](std::span<const Value> a, SdfGraph& outGraph, SdfMeshOptions& outOpt) -> Extract {
        if (a.empty()) return Extract::NoGraph;
        if (auto* h = unwrapSdfGraph(a[0])) {
            outGraph = h->graph;
            return parseMeshOptions(a.size() > 1 ? a[1] : ev::undefined(), outOpt) ? Extract::Ok : Extract::Threw;
        }
        if (ev::isObject(a[0])) {
            Value gVal = ev::getProperty(a[0], "graph");
            if (auto* h = unwrapSdfGraph(gVal)) {
                outGraph = h->graph;
                return parseMeshOptions(a[0], outOpt) ? Extract::Ok : Extract::Threw;
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
                return parseMeshOptions(a[0], outOpt) ? Extract::Ok : Extract::Threw;
            }
        }
        return Extract::NoGraph;
    };

    bindStatic("marchingCubesSDF", 2, [extractGraphAndOptions](Value, std::span<const Value> a) -> Value {
        SdfGraph graph;
        SdfMeshOptions opt;
        switch (extractGraphAndOptions(a, graph, opt)) {
            case Extract::Ok: break;
            case Extract::Threw: return ev::undefined();
            case Extract::NoGraph:
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
        switch (extractGraphAndOptions(a, graph, opt)) {
            case Extract::Ok: break;
            case Extract::Threw: return ev::undefined();
            case Extract::NoGraph:
                return ev::throwTypeError("Mesh.surfaceNetsSDF: valid SDFGraph or options required");
        }
        MeshData m = surfaceNetsFromSDF(graph, opt.dimX, opt.dimY, opt.dimZ,
                                       opt.bounds, opt.isoLevel);
        return wrapMesh(std::move(m));
    });
}

} // namespace bromesh::api
