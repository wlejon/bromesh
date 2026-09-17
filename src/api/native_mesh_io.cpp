// File I/O statics and methods of the Mesh class — OBJ, PLY, STL, VOX and
// FBX by path, the way `Mesh.loadGLTF` / `mesh.saveGLTF` already work — plus
// the PolyMesh host class: half-edge adjacency over N-gon faces, the edit
// topology a mesh editor keeps beside the triangle mesh it renders
// (bromesh/manipulation/poly_mesh.h).
//
// Paths go through the host's resolver when one is set (setPathResolver):
// a relative path then means what it means to the host's `fs`, anchored at
// the app, rather than at the process's working directory.

#include "host_mesh_internal.h"

#include <bromesh/io/obj.h>
#include <bromesh/io/ply.h>
#include <bromesh/io/stl.h>
#include <bromesh/io/vox.h>
#include <bromesh/io/fbx.h>
#include <bromesh/io/splat_ply.h>
#include <bromesh/manipulation/poly_mesh.h>

#include <filesystem>
#include <functional>

namespace bromesh::api {

HostClass g_polyMeshClass;

namespace {

std::function<std::string(const std::string&)>& pathResolver() {
    static std::function<std::string(const std::string&)> r;
    return r;
}

std::string resolvePath(const std::string& path) {
    auto& r = pathResolver();
    return r ? r(path) : path;
}

// A file to be WRITTEN does not exist yet, so the host's resolver — which
// answers by finding the file — would hand a relative path back unchanged
// and the save would land in the process's working directory rather than
// beside the directory an `fs.mkdirSync` just made under the app. Resolve
// the parent directory instead, which does exist, and put the file in it.
std::string resolveWritePath(const std::string& path) {
    namespace fs = std::filesystem;
    fs::path p(path);
    if (p.is_absolute() || !p.has_parent_path()) return resolvePath(path);
    std::string dir = resolvePath(p.parent_path().generic_string());
    return (fs::path(dir) / p.filename()).generic_string();
}

constexpr uint32_t kHostPolyMeshTag = 0x504F4C59u; // 'POLY'

struct HostPolyMesh {
    std::unique_ptr<bromesh::PolyMesh> pm = std::make_unique<bromesh::PolyMesh>();
    uint32_t tag = kHostPolyMeshTag;
};

HostPolyMesh* unwrapPolyMesh(Value v) {
    void* ptr = g_polyMeshClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostPolyMesh*>(ptr);
    return (h && h->tag == kHostPolyMeshTag) ? h : nullptr;
}

Value wrapPolyMesh(bromesh::PolyMesh pm) {
    auto h = std::make_unique<HostPolyMesh>();
    *h->pm = std::move(pm);
    return g_polyMeshClass.createInstance(std::move(h));
}

Value makeInt32Array(const int32_t* data, size_t count) {
    Value arr = ev::createTypedArray(ev::elements::Int32, static_cast<uint32_t>(count));
    if (data && count > 0) {
        ev::fillTypedArray(arr, std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(data),
                                                         count * sizeof(int32_t)));
    }
    return arr;
}

Value makeIntList(const std::vector<int32_t>& v) {
    return hostArrayOf(v.size(), [&](size_t i) { return ev::fromDouble(v[i]); });
}

Value makeFloat3(const float p[3]) {
    return hostArrayOf(3, [&](size_t i) { return ev::fromDouble(p[i]); });
}

// A group list arrives as Int32Array, Uint32Array or a plain array.
std::vector<int32_t> toInt32Vector(Value val) {
    std::vector<int32_t> out;
    ev::TypedArrayInfo info = ev::typedArrayInfo(val);
    if (info && (info.elementKind == ev::elements::Int32 || info.elementKind == ev::elements::Uint32)) {
        const auto* p = reinterpret_cast<const int32_t*>(info.data);
        out.assign(p, p + info.elementCount);
        return out;
    }
    for (uint32_t u : toUint32Vector(val)) out.push_back(static_cast<int32_t>(u));
    return out;
}

// [x,y,z] or {x,y,z}; anything else reads as the zero vector.
void readVec3Into(Value v, float out[3]) {
    out[0] = out[1] = out[2] = 0.0f;
    if (!ev::isObject(v)) return;
    Value x = ev::getProperty(v, "x");
    if (ev::isNumber(x)) {
        Value y = ev::getProperty(v, "y");
        Value z = ev::getProperty(v, "z");
        out[0] = static_cast<float>(ev::toDouble(x));
        if (ev::isNumber(y)) out[1] = static_cast<float>(ev::toDouble(y));
        if (ev::isNumber(z)) out[2] = static_cast<float>(ev::toDouble(z));
        return;
    }
    for (uint32_t i = 0; i < 3; ++i) {
        Value e = ev::getElement(v, i);
        if (ev::isNumber(e)) out[i] = static_cast<float>(ev::toDouble(e));
    }
}

void readOffset(std::span<const Value> a, size_t i, float out[3]) {
    readVec3Into(a.size() > i ? a[i] : ev::undefined(), out);
}

Value makeSplatCloud(const bromesh::GaussianSplatCloud& c) {
    ObjectBuilder o;
    o.set("positions", makeFloat32Array(c.positions.data(), c.positions.size()));
    o.set("scales", makeFloat32Array(c.scales.data(), c.scales.size()));
    o.set("rotations", makeFloat32Array(c.rotations.data(), c.rotations.size()));
    o.set("opacities", makeFloat32Array(c.opacities.data(), c.opacities.size()));
    o.set("sh", makeFloat32Array(c.sh.data(), c.sh.size()));
    o.set("shDegree", static_cast<double>(c.shDegree));
    o.set("count", static_cast<double>(c.count()));
    return o.build();
}

bool readSplatCloud(Value obj, bromesh::GaussianSplatCloud& cloud, std::string& err) {
    if (!ev::isObject(obj)) { err = "cloud must be an object"; return false; }
    cloud.positions = toFloatVector(ev::getProperty(obj, "positions"));
    cloud.scales = toFloatVector(ev::getProperty(obj, "scales"));
    cloud.rotations = toFloatVector(ev::getProperty(obj, "rotations"));
    cloud.opacities = toFloatVector(ev::getProperty(obj, "opacities"));
    cloud.sh = toFloatVector(ev::getProperty(obj, "sh"));
    Value shd = ev::getProperty(obj, "shDegree");
    cloud.shDegree = ev::isNumber(shd) ? static_cast<int>(ev::toDouble(shd)) : 0;
    if (cloud.positions.empty()) { err = "cloud has no positions"; return false; }
    if (!cloud.validate()) { err = "attribute array lengths disagree with the point count"; return false; }
    return true;
}

int32_t i32OrAt(std::span<const Value> a, size_t i, int32_t def) {
    return (a.size() > i && ev::isNumber(a[i])) ? i32At(a, i) : def;
}

bool boolOrAt(std::span<const Value> a, size_t i, bool def) {
    return (a.size() > i && !ev::isUndefined(a[i])) ? ev::toBool(a[i]) : def;
}

}  // namespace

void setPathResolver(std::function<std::string(const std::string&)> resolver) {
    pathResolver() = std::move(resolver);
}

std::string resolveMeshPath(const std::string& path) {
    return resolvePath(path);
}

std::string resolveMeshWritePath(const std::string& path) {
    return resolveWritePath(path);
}

void initMeshIo(ObjectBuilder& proto, HostClass& cls) {
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        cls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    // ---- Loaders: one Mesh per file for the single-mesh formats ---------
    bindStatic("loadOBJ", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.loadOBJ: path required");
        return wrapMesh(bromesh::loadOBJ(resolvePath(ev::toUtf8(a[0]))));
    });
    bindStatic("loadPLY", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.loadPLY: path required");
        return wrapMesh(bromesh::loadPLY(resolvePath(ev::toUtf8(a[0]))));
    });
    bindStatic("loadSTL", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.loadSTL: path required");
        return wrapMesh(bromesh::loadSTL(resolvePath(ev::toUtf8(a[0]))));
    });
    // An FBX is a scene: every mesh it holds, in file order.
    bindStatic("loadFBX", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.loadFBX: path required");
        auto meshes = bromesh::loadFBX(resolvePath(ev::toUtf8(a[0])));
        return hostArrayOf(meshes.size(), [&](size_t i) { return wrapMesh(std::move(meshes[i])); });
    });
    // A .vox is a voxel grid, not a mesh: {sizeX, sizeY, sizeZ, voxels:
    // Uint8Array (0 = empty, else palette index), palette: Float32Array RGBA
    // x 256}. Mesh it with Mesh.greedyMesh.
    bindStatic("loadVOX", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.loadVOX: path required");
        bromesh::VoxData vox = bromesh::loadVOX(resolvePath(ev::toUtf8(a[0])));
        ObjectBuilder o;
        o.set("sizeX", static_cast<double>(vox.sizeX));
        o.set("sizeY", static_cast<double>(vox.sizeY));
        o.set("sizeZ", static_cast<double>(vox.sizeZ));
        o.set("voxels", makeUint8Array(vox.voxels.data(), vox.voxels.size()));
        o.set("palette", makeFloat32Array(vox.palette, 256 * 4));
        return o.build();
    });

    // A Gaussian splat .ply is a cloud, not a mesh: the plain object
    // scene.createGaussianSplat takes ({positions, scales, rotations,
    // opacities, sh, shDegree, count}), and the same shape saves back.
    bindStatic("loadSplatPLY", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.loadSplatPLY: path required");
        return makeSplatCloud(bromesh::loadSplatPLY(resolvePath(ev::toUtf8(a[0]))));
    });
    bindStatic("saveSplatPLY", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2 || !ev::isString(a[0])) return ev::throwTypeError("Mesh.saveSplatPLY: (path, cloud) required");
        bromesh::GaussianSplatCloud cloud;
        std::string err;
        if (!readSplatCloud(a[1], cloud, err)) return ev::throwTypeError(("Mesh.saveSplatPLY: " + err).c_str());
        return ev::fromBool(bromesh::saveSplatPLY(cloud, resolveWritePath(ev::toUtf8(a[0]))));
    });

    // ---- Savers: true when the file was written --------------------------
    proto.def("saveOBJ", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.saveOBJ: not a Mesh instance");
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.saveOBJ: path required");
        return ev::fromBool(bromesh::saveOBJ(m->mesh, resolveWritePath(ev::toUtf8(a[0]))));
    });
    proto.def("savePLY", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.savePLY: not a Mesh instance");
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.savePLY: path required");
        return ev::fromBool(bromesh::savePLY(m->mesh, resolveWritePath(ev::toUtf8(a[0]))));
    });
    proto.def("saveSTL", 1, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.saveSTL: not a Mesh instance");
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.saveSTL: path required");
        return ev::fromBool(bromesh::saveSTL(m->mesh, resolveWritePath(ev::toUtf8(a[0]))));
    });
}

// ---------------------------------------------------------------------------
// PolyMesh class
// ---------------------------------------------------------------------------
void initPolyMesh(HostClass& cls) {
    cls.install("PolyMesh", 0, [](Value, std::span<const Value>) -> Value {
        // Empty by default; the static factories build populated ones.
        return g_polyMeshClass.createInstance(std::make_unique<HostPolyMesh>());
    }, [](ObjectBuilder& proto) {
        // ---- Inspection ----------------------------------------------------
        proto.accessor("vertexCount", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            return ev::fromDouble(h ? static_cast<double>(h->pm->vertexCount()) : 0.0);
        });
        proto.accessor("halfEdgeCount", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            return ev::fromDouble(h ? static_cast<double>(h->pm->halfEdgeCount()) : 0.0);
        });
        proto.accessor("faceCount", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            return ev::fromDouble(h ? static_cast<double>(h->pm->faceCount()) : 0.0);
        });
        proto.def("faceVertexCount", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.faceVertexCount: not a PolyMesh instance");
            return ev::fromDouble(h->pm->faceVertexCount(i32At(a, 0)));
        });
        proto.def("faceVertices", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.faceVertices: not a PolyMesh instance");
            return makeIntList(h->pm->faceVertices(i32At(a, 0)));
        });
        proto.def("faceHalfEdges", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.faceHalfEdges: not a PolyMesh instance");
            return makeIntList(h->pm->faceHalfEdges(i32At(a, 0)));
        });
        proto.def("getVertex", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.getVertex: not a PolyMesh instance");
            float p[3] = {0, 0, 0};
            h->pm->getVertex(i32At(a, 0), p);
            return makeFloat3(p);
        });
        proto.def("computeFaceNormal", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.computeFaceNormal: not a PolyMesh instance");
            float n[3] = {0, 0, 0};
            h->pm->computeFaceNormal(i32At(a, 0), n);
            return makeFloat3(n);
        });
        proto.def("faceGroup", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.faceGroup: not a PolyMesh instance");
            int f = i32At(a, 0);
            if (f < 0 || f >= static_cast<int>(h->pm->faceCount())) return ev::fromDouble(-1);
            return ev::fromDouble(h->pm->faces()[static_cast<size_t>(f)].group);
        });
        proto.def("setFaceGroup", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.setFaceGroup: not a PolyMesh instance");
            int f = i32At(a, 0);
            if (f >= 0 && f < static_cast<int>(h->pm->faceCount())) {
                h->pm->faces()[static_cast<size_t>(f)].group = i32At(a, 1);
            }
            return ev::undefined();
        });
        proto.def("facesInGroup", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.facesInGroup: not a PolyMesh instance");
            return makeIntList(h->pm->facesInGroup(i32At(a, 0)));
        });
        proto.def("isBoundaryVertex", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.isBoundaryVertex: not a PolyMesh instance");
            return ev::fromBool(h->pm->isBoundaryVertex(i32At(a, 0)));
        });
        proto.def("isBoundaryHalfEdge", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.isBoundaryHalfEdge: not a PolyMesh instance");
            return ev::fromBool(h->pm->isBoundaryHalfEdge(i32At(a, 0)));
        });
        auto loops = [](const std::vector<std::vector<int32_t>>& ls) {
            return hostArrayOf(ls.size(), [&](size_t i) { return makeIntList(ls[i]); });
        };
        proto.def("findFaceBoundary", 1, [loops](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.findFaceBoundary: not a PolyMesh instance");
            return loops(h->pm->findFaceBoundary(i32At(a, 0)));
        });
        proto.def("findGroupBoundary", 1, [loops](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.findGroupBoundary: not a PolyMesh instance");
            return loops(h->pm->findGroupBoundary(i32At(a, 0)));
        });

        // ---- Tessellation / validation -------------------------------------
        proto.def("tessellate", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.tessellate: not a PolyMesh instance");
            auto t = h->pm->tessellate();
            ObjectBuilder o;
            o.set("positions", makeFloat32Array(t.positions.data(), t.positions.size()));
            o.set("normals", makeFloat32Array(t.normals.data(), t.normals.size()));
            o.set("indices", makeUint32Array(t.indices.data(), t.indices.size()));
            o.set("triToFace", makeInt32Array(t.triToFace.data(), t.triToFace.size()));
            o.set("triToGroup", makeInt32Array(t.triToGroup.data(), t.triToGroup.size()));
            return o.build();
        });
        proto.def("toMesh", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.toMesh: not a PolyMesh instance");
            auto t = h->pm->tessellate();
            bromesh::MeshData md;
            md.positions = std::move(t.positions);
            md.normals = std::move(t.normals);
            md.indices = std::move(t.indices);
            return wrapMesh(std::move(md));
        });
        proto.def("validate", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.validate: not a PolyMesh instance");
            auto v = h->pm->validate();
            ObjectBuilder o;
            o.set("valid", v.valid);
            o.set("isClosed", v.isClosed);
            o.set("boundaryHalfEdges", static_cast<double>(v.boundaryHalfEdges));
            o.set("errors", hostArrayOf(v.errors.size(), [&](size_t i) { return ev::fromUtf8(v.errors[i]); }));
            return o.build();
        });

        // ---- Surgery -------------------------------------------------------
        proto.def("addVertex", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.addVertex: not a PolyMesh instance");
            return ev::fromDouble(h->pm->addVertex(static_cast<float>(numAt(a, 0)),
                                                   static_cast<float>(numAt(a, 1)),
                                                   static_cast<float>(numAt(a, 2))));
        });
        proto.def("addFace", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.addFace: not a PolyMesh instance");
            if (a.empty() || !ev::isObject(a[0])) {
                return ev::throwTypeError("PolyMesh.addFace: expects an array of vertex indices");
            }
            return ev::fromDouble(h->pm->addFace(toInt32Vector(a[0]), i32OrAt(a, 1, -1)));
        });
        proto.def("deleteFace", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.deleteFace: not a PolyMesh instance");
            h->pm->deleteFace(i32At(a, 0));
            return ev::undefined();
        });
        proto.def("translateVertex", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.translateVertex: not a PolyMesh instance");
            float o[3]; readOffset(a, 1, o);
            h->pm->translateVertex(i32At(a, 0), o);
            return ev::undefined();
        });
        proto.def("translateFace", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.translateFace: not a PolyMesh instance");
            float o[3]; readOffset(a, 1, o);
            h->pm->translateFace(i32At(a, 0), o);
            return ev::undefined();
        });
        proto.def("translateFaceWithRing", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.translateFaceWithRing: not a PolyMesh instance");
            float o[3]; readOffset(a, 1, o);
            h->pm->translateFaceWithRing(i32At(a, 0), o);
            return ev::undefined();
        });
        proto.def("extrudeFace", 5, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.extrudeFace: not a PolyMesh instance");
            float o[3]; readOffset(a, 1, o);
            auto r = h->pm->extrudeFace(i32At(a, 0), o, boolOrAt(a, 2, true),
                                        i32OrAt(a, 3, -1), i32OrAt(a, 4, -1));
            ObjectBuilder res;
            res.set("dupVerts", makeInt32Array(r.dupVerts.data(), r.dupVerts.size()));
            res.set("bridgeFaces", makeInt32Array(r.bridgeFaces.data(), r.bridgeFaces.size()));
            res.set("bridgeAdjGroup", makeInt32Array(r.bridgeAdjGroup.data(), r.bridgeAdjGroup.size()));
            res.set("backFace", static_cast<double>(r.backFace));
            return res.build();
        });
        proto.def("insetFace", 4, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.insetFace: not a PolyMesh instance");
            auto r = h->pm->insetFace(i32At(a, 0), static_cast<float>(numAt(a, 1)),
                                      boolOrAt(a, 2, false), i32OrAt(a, 3, -1));
            ObjectBuilder res;
            res.set("innerFace", static_cast<double>(r.innerFace));
            res.set("innerVerts", makeInt32Array(r.innerVerts.data(), r.innerVerts.size()));
            res.set("bridgeFaces", makeInt32Array(r.bridgeFaces.data(), r.bridgeFaces.size()));
            return res.build();
        });
        proto.def("splitEdge", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.splitEdge: not a PolyMesh instance");
            float p[3];
            const bool hasPos = a.size() > 1 && ev::isObject(a[1]);
            if (hasPos) readOffset(a, 1, p);
            return ev::fromDouble(h->pm->splitEdge(i32At(a, 0), hasPos ? p : nullptr));
        });
        proto.def("flipEdge", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.flipEdge: not a PolyMesh instance");
            return ev::fromBool(h->pm->flipEdge(i32At(a, 0)));
        });
        proto.def("collapseEdge", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.collapseEdge: not a PolyMesh instance");
            float p[3];
            const bool hasPos = a.size() > 1 && ev::isObject(a[1]);
            if (hasPos) readOffset(a, 1, p);
            return ev::fromBool(h->pm->collapseEdge(i32At(a, 0), hasPos ? p : nullptr));
        });
        proto.def("rematchTwins", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.rematchTwins: not a PolyMesh instance");
            h->pm->rematchTwins();
            return ev::undefined();
        });
        proto.def("mergeFacesByGroup", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.mergeFacesByGroup: not a PolyMesh instance");
            h->pm->mergeFacesByGroup();
            return ev::undefined();
        });
        proto.def("compact", 0, [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapPolyMesh(self);
            if (!h) return ev::throwTypeError("PolyMesh.compact: not a PolyMesh instance");
            h->pm->compact();
            return ev::undefined();
        });
    });

    // ---- Factories -----------------------------------------------------------
    cls.setStatic("fromMeshData", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("PolyMesh.fromMeshData(positions, indices[, triToGroup])");
        std::vector<float> positions = toFloatVector(a[0]);
        std::vector<uint32_t> indices = toUint32Vector(a[1]);
        std::vector<int32_t> triToGroup;
        if (a.size() > 2 && ev::isObject(a[2])) triToGroup = toInt32Vector(a[2]);
        return wrapPolyMesh(bromesh::PolyMesh::fromMeshData(positions, indices, triToGroup));
    }, 3, "fromMeshData"));
    cls.setStatic("fromMesh", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        auto* m = a.empty() ? nullptr : unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("PolyMesh.fromMesh: expects a Mesh");
        std::vector<int32_t> triToGroup;
        if (a.size() > 1 && ev::isObject(a[1])) triToGroup = toInt32Vector(a[1]);
        return wrapPolyMesh(bromesh::PolyMesh::fromMeshData(m->mesh.positions, m->mesh.indices, triToGroup));
    }, 2, "fromMesh"));
    cls.setStatic("fromPolygon", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("PolyMesh.fromPolygon(positionsXYZ, normal[, group])");
        std::vector<float> positions = toFloatVector(a[0]);
        float n[3];
        readVec3Into(a[1], n);
        return wrapPolyMesh(bromesh::PolyMesh::fromPolygon(positions, n, i32OrAt(a, 2, 0)));
    }, 3, "fromPolygon"));
    cls.setStatic("fromPolygons", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) {
            return ev::throwTypeError("PolyMesh.fromPolygons(positions, polyVerts, polyOffsets[, faceGroups])");
        }
        std::vector<float> positions = toFloatVector(a[0]);
        std::vector<uint32_t> polyVerts = toUint32Vector(a[1]);
        std::vector<uint32_t> polyOffsets = toUint32Vector(a[2]);
        std::vector<int32_t> groups;
        if (a.size() > 3 && ev::isObject(a[3])) groups = toInt32Vector(a[3]);
        return wrapPolyMesh(bromesh::PolyMesh::fromPolygons(positions, polyVerts, polyOffsets, groups));
    }, 4, "fromPolygons"));
}

}  // namespace bromesh::api
