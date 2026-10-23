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
        if (ev::isNumber(c)) out = finiteFloat(c, out);
    };
    Value x = ev::getProperty(obj.get(), "x");
    if (ev::isNumber(x)) {
        r.x = finiteFloat(x);
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
            float s = finiteFloat(sizeVal, 2.0f);
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
    bool enableMacroBlockPruning = true;
    std::string algorithm = "marchingCubes";
};

static bool extractCustomExprSource(Value v, std::string& outExpr) {
    if (ev::isString(v)) {
        outExpr = ev::toUtf8(v);
        return !outExpr.empty();
    }
    if (ev::isFunction(v)) {
        ev::Persistent fn(v);
        ev::Persistent toStrFn(ev::getProperty(fn.get(), "toString"));
        if (ev::isFunction(toStrFn.get())) {
            ev::CallResult r = ev::call(toStrFn.get(), fn.get(), {});
            if (!r.thrown && ev::isString(r.value)) {
                outExpr = ev::toUtf8(r.value);
                return !outExpr.empty();
            }
        }
    }
    return false;
}

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
    if (ev::isNumber(isoVal)) opt.isoLevel = finiteFloat(isoVal, opt.isoLevel);

    Value cbVal = ev::getProperty(opts.get(), "closeBoundary");
    if (ev::isBool(cbVal)) opt.closeBoundary = ev::toBool(cbVal);

    Value cgVal = ev::getProperty(opts.get(), "computeGradients");
    if (ev::isBool(cgVal)) opt.computeGradients = ev::toBool(cgVal);

    Value mbpVal = ev::getProperty(opts.get(), "macroBlockPruning");
    if (ev::isBool(mbpVal)) opt.enableMacroBlockPruning = ev::toBool(mbpVal);

    Value embpVal = ev::getProperty(opts.get(), "enableMacroBlockPruning");
    if (ev::isBool(embpVal)) opt.enableMacroBlockPruning = ev::toBool(embpVal);

    Value algVal = ev::getProperty(opts.get(), "algorithm");
    if (ev::isString(algVal)) {
        std::string alg = ev::toUtf8(algVal);
        if (alg == "marchingCubes" || alg == "marching_cubes" || alg == "mc") {
            opt.algorithm = "marchingCubes";
        } else if (alg == "surfaceNets" || alg == "surface_nets" || alg == "sn") {
            opt.algorithm = "surfaceNets";
        } else {
            ev::throwRangeError("SDF mesh options: unknown algorithm '" + alg + "'");
            return false;
        }
    }

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

        auto customExprFn = [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.customExpr: invalid instance");
            std::string expr;
            if (a.empty() || !extractCustomExprSource(a[0], expr)) {
                return ev::throwTypeError("SDFGraph.customExpr: formula string or function required");
            }
            float lip = (a.size() > 1 && ev::isNumber(a[1])) ? finiteFloat(a[1], 1.0f) : 1.0f;
            try {
                return ev::fromDouble(h->graph.customExpr(expr, lip));
            } catch (const std::exception& e) {
                return ev::throwRangeError(e.what());
            }
        };
        proto.def("customExpr", 2, customExprFn);
        proto.def("customPrimitive", 2, customExprFn);

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

        proto.def("evalPoint", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.evalPoint: invalid instance");
            size_t idx = 0;
            bromath::Vec3 p = parseVec3(a, idx);
            CompiledSdfKernel k = JitSdfCompiler::instance().compile(h->graph);
            float d = k.pointFn ? k.pointFn(p.x, p.y, p.z) : 0.0f;
            return ev::fromDouble(d);
        });

        proto.def("queryContact", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.queryContact: invalid instance");
            size_t idx = 0;
            bromath::Vec3 p = parseVec3(a, idx);
            CompiledSdfKernel k = JitSdfCompiler::instance().compile(h->graph);
            SdfContact c = k.queryContact(p.x, p.y, p.z);

            ObjectBuilder normObj;
            normObj.set("x", c.normal.x);
            normObj.set("y", c.normal.y);
            normObj.set("z", c.normal.z);

            ObjectBuilder ptObj;
            ptObj.set("x", c.closestSurfacePoint.x);
            ptObj.set("y", c.closestSurfacePoint.y);
            ptObj.set("z", c.closestSurfacePoint.z);

            ObjectBuilder gradObj;
            gradObj.set("x", c.gradient.x);
            gradObj.set("y", c.gradient.y);
            gradObj.set("z", c.gradient.z);

            ObjectBuilder res;
            res.set("distance", c.distance);
            res.set("normal", normObj.build());
            res.set("closestSurfacePoint", ptObj.build());
            res.set("gradient", gradObj.build());
            return res.build();
        });

        proto.def("evalBlockBounds", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.evalBlockBounds: invalid instance");
            bromath::AABB3 b = parseBounds(!a.empty() ? a[0] : ev::undefined());
            float dMin = 0.0f, dMax = 0.0f;
            evalBlockBounds(h->graph, -1, b, dMin, dMax);
            ObjectBuilder res;
            res.set("min", dMin);
            res.set("max", dMax);
            return res.build();
        });

        proto.def("evaluateVolumeStats", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.evaluateVolumeStats: invalid instance");
            SdfMeshOptions opt;
            if (!parseMeshOptions(!a.empty() ? a[0] : ev::undefined(), opt)) return ev::undefined();
            SDFVolume vol = JitSdfCompiler::instance().evaluateVolume(
                h->graph, opt.dimX, opt.dimY, opt.dimZ, opt.bounds,
                -1, true, opt.isoLevel, opt.enableMacroBlockPruning);
            size_t activeBlocks = vol.activeBlockCount;
            size_t totalBlocks = vol.totalBlockCount;
            size_t prunedBlocks = (totalBlocks >= activeBlocks) ? (totalBlocks - activeBlocks) : 0;
            double pruningRatio = (totalBlocks > 0)
                                      ? (static_cast<double>(prunedBlocks) / static_cast<double>(totalBlocks))
                                      : 0.0;
            ObjectBuilder res;
            res.set("activeBlocks", static_cast<double>(activeBlocks));
            res.set("totalBlocks", static_cast<double>(totalBlocks));
            res.set("prunedBlocks", static_cast<double>(prunedBlocks));
            res.set("pruningRatio", pruningRatio);
            res.set("dimX", static_cast<double>(vol.dimX));
            res.set("dimY", static_cast<double>(vol.dimY));
            res.set("dimZ", static_cast<double>(vol.dimZ));
            return res.build();
        });

        proto.def("compileSpirv", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.compileSpirv: invalid instance");
            std::vector<uint32_t> spv = JitSdfCompiler::instance().compileSpirv(h->graph);
            return makeUint32Array(spv.data(), spv.size());
        });

        proto.def("compilePtx", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.compilePtx: invalid instance");
            std::string ptx = JitSdfCompiler::instance().compilePtx(h->graph);
            return ev::fromUtf8(ptx);
        });

        // Meshing methods on SDFGraph instance
        proto.def("marchingCubes", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.marchingCubes: invalid instance");
            SdfMeshOptions opt;
            if (!parseMeshOptions(!a.empty() ? a[0] : ev::undefined(), opt)) return ev::undefined();
            MeshData m = marchingCubesFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                              opt.bounds, opt.isoLevel,
                                              opt.closeBoundary, opt.computeGradients,
                                              opt.enableMacroBlockPruning);
            return wrapMesh(std::move(m));
        });

        proto.def("surfaceNets", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.surfaceNets: invalid instance");
            SdfMeshOptions opt;
            if (!parseMeshOptions(!a.empty() ? a[0] : ev::undefined(), opt)) return ev::undefined();
            MeshData m = surfaceNetsFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                           opt.bounds, opt.isoLevel,
                                           opt.enableMacroBlockPruning);
            return wrapMesh(std::move(m));
        });

        proto.def("toMesh", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapSdfGraph(self);
            if (!h) return ev::throwTypeError("SDFGraph.toMesh: invalid instance");
            SdfMeshOptions opt;
            if (!parseMeshOptions(!a.empty() ? a[0] : ev::undefined(), opt)) return ev::undefined();
            MeshData m = (opt.algorithm == "surfaceNets")
                             ? surfaceNetsFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                                  opt.bounds, opt.isoLevel,
                                                  opt.enableMacroBlockPruning)
                             : marchingCubesFromSDF(h->graph, opt.dimX, opt.dimY, opt.dimZ,
                                                    opt.bounds, opt.isoLevel,
                                                    opt.closeBoundary, opt.computeGradients,
                                                    opt.enableMacroBlockPruning);
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
        if (ev::isFunction(a[0])) {
            std::string expr;
            if (!extractCustomExprSource(a[0], expr)) return Extract::NoGraph;
            try {
                outGraph.customExpr(expr, 1.0f);
            } catch (const std::exception& e) {
                ev::throwRangeError(e.what());
                return Extract::Threw;
            }
            return parseMeshOptions(a.size() > 1 ? a[1] : ev::undefined(), outOpt) ? Extract::Ok : Extract::Threw;
        }
        if (ev::isObject(a[0])) {
            ev::Persistent obj(a[0]);
            Value gVal = ev::getProperty(obj.get(), "graph");
            if (auto* h = unwrapSdfGraph(gVal)) {
                outGraph = h->graph;
                return parseMeshOptions(obj.get(), outOpt) ? Extract::Ok : Extract::Threw;
            }
            // Declarative shorthand: e.g. { type: "sphere", radius: 1.0 }
            Value typeVal = ev::getProperty(obj.get(), "type");
            if (!ev::isString(typeVal)) typeVal = ev::getProperty(obj.get(), "op");
            if (ev::isString(typeVal)) {
                std::string t = ev::toUtf8(typeVal);
                if (t == "sphere") {
                    Value rVal = ev::getProperty(obj.get(), "radius");
                    float r = finiteFloat(rVal, 1.0f);
                    outGraph.sphere(r);
                } else if (t == "box") {
                    Value extVal = ev::getProperty(obj.get(), "halfExtents");
                    bromath::Vec3 ext = ev::isObject(extVal) ? readVec3Param(extVal) : bromath::Vec3{1.0f, 1.0f, 1.0f};
                    outGraph.box(ext);
                } else if (t == "torus") {
                    float majR = finiteFloat(ev::getProperty(obj.get(), "majorRadius"), 1.0f);
                    float minR = finiteFloat(ev::getProperty(obj.get(), "minorRadius"), 0.3f);
                    outGraph.torus(majR > 0 ? majR : 1.0f, minR > 0 ? minR : 0.3f);
                } else if (t == "cylinder") {
                    float r = finiteFloat(ev::getProperty(obj.get(), "radius"), 0.5f);
                    float hh = finiteFloat(ev::getProperty(obj.get(), "halfHeight"), 1.0f);
                    outGraph.cylinder(r > 0 ? r : 0.5f, hh > 0 ? hh : 1.0f);
                } else if (t == "custom" || t == "customExpr" || t == "customPrimitive") {
                    Value exprVal = ev::getProperty(obj.get(), "expr");
                    if (ev::isUndefined(exprVal)) exprVal = ev::getProperty(obj.get(), "formula");
                    if (ev::isUndefined(exprVal)) exprVal = ev::getProperty(obj.get(), "fn");
                    if (ev::isUndefined(exprVal)) exprVal = ev::getProperty(obj.get(), "callback");
                    std::string expr;
                    if (!extractCustomExprSource(exprVal, expr)) return Extract::NoGraph;
                    float lip = finiteFloat(ev::getProperty(obj.get(), "lipschitzBound"), 1.0f);
                    try {
                        outGraph.customExpr(expr, lip);
                    } catch (const std::exception& e) {
                        ev::throwRangeError(e.what());
                        return Extract::Threw;
                    }
                }
                return parseMeshOptions(obj.get(), outOpt) ? Extract::Ok : Extract::Threw;
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
                                          opt.closeBoundary, opt.computeGradients,
                                          opt.enableMacroBlockPruning);
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
                                       opt.bounds, opt.isoLevel,
                                       opt.enableMacroBlockPruning);
        return wrapMesh(std::move(m));
    });

    bindStatic("fromSDF", 2, [extractGraphAndOptions](Value, std::span<const Value> a) -> Value {
        SdfGraph graph;
        SdfMeshOptions opt;
        switch (extractGraphAndOptions(a, graph, opt)) {
            case Extract::Ok: break;
            case Extract::Threw: return ev::undefined();
            case Extract::NoGraph:
                return ev::throwTypeError("Mesh.fromSDF: valid SDFGraph or options required");
        }
        MeshData m = (opt.algorithm == "surfaceNets")
                         ? surfaceNetsFromSDF(graph, opt.dimX, opt.dimY, opt.dimZ,
                                              opt.bounds, opt.isoLevel,
                                              opt.enableMacroBlockPruning)
                         : marchingCubesFromSDF(graph, opt.dimX, opt.dimY, opt.dimZ,
                                                opt.bounds, opt.isoLevel,
                                                opt.closeBoundary, opt.computeGradients,
                                                opt.enableMacroBlockPruning);
        return wrapMesh(std::move(m));
    });
}

} // namespace bromesh::api
