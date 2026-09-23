#include "host_mesh_internal.h"
#include <string>

namespace bromesh::api {

HostClass g_meshClass;

namespace {

Value meshConstructor(Value, std::span<const Value> args) {
    auto hm = std::make_unique<HostMesh>();
    if (args.empty()) {
        return g_meshClass.createInstance(std::move(hm));
    }

    // Check if args[0] is options object
    if (ev::isObject(args[0]) && !ev::isTypedArray(args[0])) {
        // args[0] is a rooted slot; each field is read right before it is
        // consumed, since every read (and every conversion) may allocate.
        auto field = [&](const char* name) { return ev::getProperty(args[0], name); };
        auto present = [](Value v) { return !ev::isUndefined(v) && !ev::isNull(v); };

        Value pVal = field("positions");
        if (present(pVal)) {
            if (!ev::isTypedArray(pVal)) {
                return ev::throwTypeError("expected a Float32Array, got a non-typed-array object");
            }
            ev::TypedArrayInfo info = ev::typedArrayInfo(pVal);
            if (info.elementKind != ev::elements::Float32) {
                return ev::throwTypeError("expected a Float32Array, got a Float64Array");
            }
            hm->mesh.positions = toFloatVector(pVal);
        }
        if (Value v = field("normals"); present(v)) hm->mesh.normals = toFloatVector(v);
        if (Value v = field("uvs"); present(v)) hm->mesh.uvs = toFloatVector(v);
        if (Value v = field("colors"); present(v)) hm->mesh.colors = toFloatVector(v);
        if (Value v = field("indices"); present(v)) hm->mesh.indices = toUint32Vector(v);
    } else {
        // Positional arguments: positions, normals, uvs, colors, indices
        if (args.size() > 0 && !ev::isUndefined(args[0])) hm->mesh.positions = toFloatVector(args[0]);
        if (args.size() > 1 && !ev::isUndefined(args[1])) hm->mesh.normals = toFloatVector(args[1]);
        if (args.size() > 2 && !ev::isUndefined(args[2])) hm->mesh.uvs = toFloatVector(args[2]);
        if (args.size() > 3 && !ev::isUndefined(args[3])) hm->mesh.colors = toFloatVector(args[3]);
        if (args.size() > 4 && !ev::isUndefined(args[4])) hm->mesh.indices = toUint32Vector(args[4]);
    }

    if (hm->mesh.positions.size() % 3 != 0) {
        return ev::throwTypeError("Mesh: positions length must be a multiple of 3");
    }
    const size_t verts = hm->mesh.vertexCount();
    if (!hm->mesh.normals.empty() && hm->mesh.normals.size() != verts * 3) {
        return ev::throwTypeError("Mesh: normals must hold one xyz per vertex");
    }
    if (!hm->mesh.uvs.empty() && hm->mesh.uvs.size() != verts * 2) {
        return ev::throwTypeError("Mesh: uvs must hold one uv per vertex");
    }
    if (!hm->mesh.colors.empty() && hm->mesh.colors.size() != verts * 4) {
        return ev::throwTypeError("Mesh: colors must hold one rgba per vertex");
    }
    if (hm->mesh.indices.size() % 3 != 0) {
        return ev::throwTypeError("Mesh: indices length must be a multiple of 3");
    }

    return g_meshClass.createInstance(std::move(hm));
}

} // namespace

void initMeshCore(ObjectBuilder& proto, HostClass& cls) {
    // ---- Accessors ---------------------------------------------------------
    proto.accessor("positions",
        [](Value self, std::span<const Value>) -> Value {
            auto* m = unwrapMesh(self);
            if (!m) return ev::undefined();
            return makeFloat32Array(m->mesh.positions.data(), m->mesh.positions.size());
        },
        [](Value self, std::span<const Value> a) -> Value {
            auto* m = unwrapMesh(self);
            if (!m || a.empty()) return ev::undefined();
            std::vector<float> v = toFloatVector(a[0]);
            if (v.size() % 3 != 0) return ev::throwTypeError("positions length must be a multiple of 3");
            m->mesh.positions = std::move(v);
            return ev::undefined();
        });

    proto.accessor("normals",
        [](Value self, std::span<const Value>) -> Value {
            auto* m = unwrapMesh(self);
            if (!m) return ev::undefined();
            return makeFloat32Array(m->mesh.normals.data(), m->mesh.normals.size());
        },
        [](Value self, std::span<const Value> a) -> Value {
            auto* m = unwrapMesh(self);
            if (!m || a.empty()) return ev::undefined();
            std::vector<float> v = toFloatVector(a[0]);
            if (!v.empty() && v.size() != m->mesh.vertexCount() * 3) {
                return ev::throwTypeError("normals must hold one xyz per vertex");
            }
            m->mesh.normals = std::move(v);
            return ev::undefined();
        });

    proto.accessor("uvs",
        [](Value self, std::span<const Value>) -> Value {
            auto* m = unwrapMesh(self);
            if (!m) return ev::undefined();
            return makeFloat32Array(m->mesh.uvs.data(), m->mesh.uvs.size());
        },
        [](Value self, std::span<const Value> a) -> Value {
            auto* m = unwrapMesh(self);
            if (!m || a.empty()) return ev::undefined();
            std::vector<float> v = toFloatVector(a[0]);
            if (!v.empty() && v.size() != m->mesh.vertexCount() * 2) {
                return ev::throwTypeError("uvs must hold one uv per vertex");
            }
            m->mesh.uvs = std::move(v);
            return ev::undefined();
        });

    proto.accessor("colors",
        [](Value self, std::span<const Value>) -> Value {
            auto* m = unwrapMesh(self);
            if (!m) return ev::undefined();
            return makeFloat32Array(m->mesh.colors.data(), m->mesh.colors.size());
        },
        [](Value self, std::span<const Value> a) -> Value {
            auto* m = unwrapMesh(self);
            if (!m || a.empty()) return ev::undefined();
            std::vector<float> v = toFloatVector(a[0]);
            if (!v.empty() && v.size() != m->mesh.vertexCount() * 4) {
                return ev::throwTypeError("colors must hold one rgba per vertex");
            }
            m->mesh.colors = std::move(v);
            return ev::undefined();
        });

    proto.accessor("indices",
        [](Value self, std::span<const Value>) -> Value {
            auto* m = unwrapMesh(self);
            if (!m) return ev::undefined();
            return makeUint32Array(m->mesh.indices.data(), m->mesh.indices.size());
        },
        [](Value self, std::span<const Value> a) -> Value {
            auto* m = unwrapMesh(self);
            if (!m || a.empty()) return ev::undefined();
            std::vector<uint32_t> v = toUint32Vector(a[0]);
            if (v.size() % 3 != 0) return ev::throwTypeError("indices length must be a multiple of 3");
            const size_t verts = m->mesh.vertexCount();
            for (uint32_t idx : v) {
                if (idx >= verts) return ev::throwRangeError("Mesh.indices: index " + std::to_string(idx) + " is out of range");
            }
            m->mesh.indices = std::move(v);
            return ev::undefined();
        });

    proto.accessor("vertexCount", [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        return ev::fromDouble(m ? static_cast<double>(m->mesh.vertexCount()) : 0.0);
    });

    proto.accessor("triangleCount", [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        return ev::fromDouble(m ? static_cast<double>(m->mesh.triangleCount()) : 0.0);
    });

    proto.accessor("hasNormals", [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        return ev::fromBool(m ? m->mesh.hasNormals() : false);
    });

    proto.accessor("hasUVs", [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        return ev::fromBool(m ? m->mesh.hasUVs() : false);
    });

    proto.accessor("hasColors", [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        return ev::fromBool(m ? m->mesh.hasColors() : false);
    });

    proto.accessor("empty", [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        return ev::fromBool(m ? m->mesh.empty() : true);
    });

    // ---- Clone & Bounds ----------------------------------------------------
    proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.clone: not a Mesh instance");
        return wrapMesh(m->mesh);
    });

    proto.def("bounds", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.bounds: not a Mesh instance");
        const bromath::AABB3 bb = bromesh::computeBBox(m->mesh);
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
    proto.def("computeBBox", 0, [](Value self, std::span<const Value> a) -> Value {
        return callMethod(self, "bounds", a);
    });
    proto.def("computeVolume", 0, [](Value self, std::span<const Value> a) -> Value {
        return callMethod(self, "volume", a);
    });
    proto.def("computeSurfaceArea", 0, [](Value self, std::span<const Value> a) -> Value {
        return callMethod(self, "surfaceArea", a);
    });

    // ---- In-place transforms -----------------------------------------------
    proto.def("translate", 3, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.translate: not a Mesh instance");
        ArgReader r(a);
        bromesh::translateMesh(m->mesh,
            static_cast<float>(r.getDouble(0, 0.0)),
            static_cast<float>(r.getDouble(1, 0.0)),
            static_cast<float>(r.getDouble(2, 0.0)));
        return self;
    });

    proto.def("scale", 3, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.scale: not a Mesh instance");
        ArgReader r(a);
        double sx = r.getDouble(0, 1.0);
        double sy = r.has(1) ? r.getDouble(1, sx) : sx;
        double sz = r.has(2) ? r.getDouble(2, sx) : sx;
        bromesh::scaleMesh(m->mesh, static_cast<float>(sx), static_cast<float>(sy), static_cast<float>(sz));
        return self;
    });

    proto.def("rotate", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.rotate: not a Mesh instance");
        ArgReader r(a);
        bromesh::rotateMesh(m->mesh,
            static_cast<float>(r.getDouble(0, 0.0)),
            static_cast<float>(r.getDouble(1, 1.0)),
            static_cast<float>(r.getDouble(2, 0.0)),
            static_cast<float>(r.getDouble(3, 0.0)));
        return self;
    });

    proto.def("center", 0, [](Value self, std::span<const Value>) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.center: not a Mesh instance");
        bromesh::centerMesh(m->mesh);
        return self;
    });

    proto.def("fitToBox", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.fitToBox: not a Mesh instance");
        if (m->mesh.empty()) return self;
        double size = numAt(a, 0);
        if (size <= 0.0) size = 1.0;
        bromesh::centerMesh(m->mesh);
        const bromath::AABB3 bb = bromesh::computeBBox(m->mesh);
        const float ex = bb.max.x - bb.min.x, ey = bb.max.y - bb.min.y, ez = bb.max.z - bb.min.z;
        const float longest = std::fmax(ex, std::fmax(ey, ez));
        if (longest > 0.0f) bromesh::scaleMesh(m->mesh, static_cast<float>(size) / longest);
        return self;
    });

    proto.def("transform", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.transform: not a Mesh instance");
        if (a.empty()) return self;
        ev::Persistent selfP(self);  // toFloatVector may allocate
        std::vector<float> mat = toFloatVector(a[0]);
        if (mat.size() != 16) {
            return ev::throwTypeError("Mesh.transform: matrix must have 16 elements, got " + std::to_string(mat.size()));
        }
        bromesh::transformMesh(m->mesh, mat.data());
        return selfP.get();
    });


    proto.def("mirror", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.mirror: not a Mesh instance");
        int axis = 0;
        if (!countArg(a, 0, "Mesh.mirror: axis (0 X, 1 Y, 2 Z)", 0, 2, axis)) return ev::undefined();
        bromesh::mirrorMesh(m->mesh, axis);
        return self;
    });

    // ---- Static Primitive Factories ----------------------------------------
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        cls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindStatic("box", 3, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float w = static_cast<float>(r.getDouble(0, 0.5));
        float h = static_cast<float>(r.has(1) ? r.getDouble(1, w) : w);
        float d = static_cast<float>(r.has(2) ? r.getDouble(2, w) : w);
        return wrapMesh(bromesh::box(w, h, d));
    });

    bindStatic("sphere", 3, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 1.0));
        int segs = 16, rings = 12;
        if (!countArg(a, 1, "Mesh.sphere: segments", 3, kMaxAxis, segs) ||
            !countArg(a, 2, "Mesh.sphere: rings", 2, kMaxAxis, rings)) {
            return ev::undefined();
        }
        return wrapMesh(bromesh::sphere(radius, segs, rings));
    });

    bindStatic("cylinder", 3, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 0.5));
        float halfH = static_cast<float>(r.getDouble(1, 1.0));
        int segs = 16;
        if (!countArg(a, 2, "Mesh.cylinder: segments", 3, kMaxAxis, segs)) return ev::undefined();
        return wrapMesh(bromesh::cylinder(radius, halfH, segs));
    });

    bindStatic("capsule", 4, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 0.5));
        float halfH = static_cast<float>(r.getDouble(1, 1.0));
        int segs = 16, rings = 8;
        if (!countArg(a, 2, "Mesh.capsule: segments", 3, kMaxAxis, segs) ||
            !countArg(a, 3, "Mesh.capsule: rings", 1, kMaxAxis, rings)) {
            return ev::undefined();
        }
        return wrapMesh(bromesh::capsule(radius, halfH, segs, rings));
    });

    bindStatic("cone", 5, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 0.5));
        float height = static_cast<float>(r.getDouble(1, 1.0));
        int segs = 16, stacks = 4;
        if (!countArg(a, 2, "Mesh.cone: segments", 3, kMaxAxis, segs) ||
            !countArg(a, 3, "Mesh.cone: stacks", 1, kMaxAxis, stacks)) {
            return ev::undefined();
        }
        bool capped = r.has(4) ? r.getBool(4, true) : true;
        return wrapMesh(bromesh::cone(radius, height, segs, stacks, capped));
    });

    bindStatic("plane", 4, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float w = static_cast<float>(r.getDouble(0, 1.0));
        float h = static_cast<float>(r.has(1) ? r.getDouble(1, w) : w);
        int segW = 1, segH = 1;
        if (!countArg(a, 2, "Mesh.plane: segmentsW", 1, kMaxAxis, segW) ||
            !countArg(a, 3, "Mesh.plane: segmentsH", 1, kMaxAxis, segH)) {
            return ev::undefined();
        }
        return wrapMesh(bromesh::plane(w, h, segW, segH));
    });

    bindStatic("torus", 4, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 1.0));
        float tubeRadius = static_cast<float>(r.getDouble(1, 0.3));
        int segs = 24, tubeSegs = 12;
        if (!countArg(a, 2, "Mesh.torus: segments", 3, kMaxAxis, segs) ||
            !countArg(a, 3, "Mesh.torus: tubeSegments", 3, kMaxAxis, tubeSegs)) {
            return ev::undefined();
        }
        return wrapMesh(bromesh::torus(radius, tubeRadius, segs, tubeSegs));
    });

    bindStatic("icosahedron", 1, [](Value, std::span<const Value> a) -> Value {
        float r = static_cast<float>(numAt(a, 0));
        bromesh::MeshData m = bromesh::icosahedron();
        if (r > 0.0f && r != 1.0f) bromesh::scaleMesh(m, r);
        return wrapMesh(std::move(m));
    });

    bindStatic("dodecahedron", 1, [](Value, std::span<const Value> a) -> Value {
        float r = static_cast<float>(numAt(a, 0));
        bromesh::MeshData m = bromesh::dodecahedron();
        if (r > 0.0f && r != 1.0f) bromesh::scaleMesh(m, r);
        return wrapMesh(std::move(m));
    });

    bindStatic("octahedron", 1, [](Value, std::span<const Value> a) -> Value {
        float r = static_cast<float>(numAt(a, 0));
        bromesh::MeshData m = bromesh::octahedron();
        if (r > 0.0f && r != 1.0f) bromesh::scaleMesh(m, r);
        return wrapMesh(std::move(m));
    });

    bindStatic("tetrahedron", 1, [](Value, std::span<const Value> a) -> Value {
        float r = static_cast<float>(numAt(a, 0));
        bromesh::MeshData m = bromesh::tetrahedron();
        if (r > 0.0f && r != 1.0f) bromesh::scaleMesh(m, r);
        return wrapMesh(std::move(m));
    });

    auto disc = [](std::string_view fn, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 1.0));
        int segs = 16;
        if (!countArg(a, 1, std::string(fn) + ": segments", 3, kMaxAxis, segs)) return ev::undefined();
        return wrapMesh(bromesh::disc(radius > 0.0f ? radius : 1.0f, segs));
    };
    bindStatic("disk", 2, [disc](Value, std::span<const Value> a) -> Value { return disc("Mesh.disk", a); });
    bindStatic("disc", 2, [disc](Value, std::span<const Value> a) -> Value { return disc("Mesh.disc", a); });

#if BROMESH_HAS_PAR_SHAPES
    bindStatic("geodesicSphere", 2, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 1.0));
        int subdiv = 2;  // 20 * 4^subdivisions faces
        if (!countArg(a, 1, "Mesh.geodesicSphere: subdivisions", 0, kMaxSubdivisions, subdiv)) return ev::undefined();
        return wrapMesh(bromesh::geodesicSphere(radius > 0.0f ? radius : 1.0f, subdiv));
    });

    bindStatic("rock", 3, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        float radius = static_cast<float>(r.getDouble(0, 1.0));
        int seed = r.getInt(1, 1);
        int subdiv = 2;
        if (!countArg(a, 2, "Mesh.rock: subdivisions", 0, kMaxSubdivisions, subdiv)) return ev::undefined();
        return wrapMesh(bromesh::rock(radius > 0.0f ? radius : 1.0f, seed, subdiv));
    });
#endif

    // blob, tube, sweep and the plant/branch/L-system statics live in
    // native_mesh_plants.cpp (initMeshPlants).

    bindStatic("heightmapGrid", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.heightmapGrid: heights required");
        std::vector<float> heights = toFloatVector(a[0]);
        ArgReader r(a);
        int gw = 0, gh = 0, b = 0;
        if (!countArg(a, 1, "Mesh.heightmapGrid: gridW", 0, kMaxInt32, gw) ||
            !countArg(a, 2, "Mesh.heightmapGrid: gridH", 0, kMaxInt32, gh) ||
            !countArg(a, 4, "Mesh.heightmapGrid: border", 0, kMaxAxis, b)) {
            return ev::undefined();
        }
        double cellSize = r.getDouble(3, 1.0);

        if (gw <= 0 || gh <= 0) return ev::throwTypeError("Mesh.heightmapGrid: width and height must be positive");
        // In double: (gw + 2b) * (gh + 2b) overflows an int long before the
        // heights array could hold it.
        const double expected = (static_cast<double>(gw) + 2.0 * b) * (static_cast<double>(gh) + 2.0 * b);
        if (static_cast<double>(heights.size()) < expected) {
            return ev::throwTypeError("Mesh.heightmapGrid: heights array too short");
        }
        return wrapMesh(bromesh::heightmapGrid(heights.data(), gw, gh, static_cast<float>(cellSize > 0.0 ? cellSize : 1.0f), b));
    });
}

void ensureMeshClassesInstalled() {
    // Once per THREAD, not per process: a class's constructor and prototype
    // are the installing thread's (host_class.h), so a Worker realm installs
    // its own.
    static thread_local bool installed = false;
    if (installed) return;
    installed = true;

    g_meshClass.install("Mesh", 1, meshConstructor, [](ObjectBuilder& proto) {
        initMeshCore(proto, g_meshClass);
        initMeshOps(proto, g_meshClass);
        initMeshAnalysis(proto, g_meshClass);
        initMeshPlants(g_meshClass);
        initMeshIo(proto, g_meshClass);
        initMeshSdf(proto, g_meshClass);
    });

    auto objVal = ev::globalValue("Object");
    if (objVal.found && ev::isObject(objVal.value)) {
        Rooted createFn(ev::getProperty(objVal.value, "create"));
        if (ev::isFunction(createFn)) {
            Value baseProto = g_meshClass.prototype();
            Value subProto = ev::call(createFn.get(), ev::undefined(), std::span<const Value>(&baseProto, 1)).value;
            if (ev::isObject(subProto)) {
                g_meshClass.setInstancePrototype(subProto);
            }
        }
    }

    initMeshBvh(g_meshBvhClass);
    initProgressiveMesh(g_progressiveMeshClass);
    initCapsuleField(g_capsuleFieldClass);
    initLSystem(g_lsystemClass);
    initPolyMesh(g_polyMeshClass);
    initSdfGraph(g_sdfGraphClass);
}

} // namespace bromesh::api
