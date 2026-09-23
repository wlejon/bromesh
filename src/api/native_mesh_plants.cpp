// Plant / procedural surface of the Mesh class: sweeps, leaf and flower cards,
// space-colonization branch trees, leaf scattering, obstacle fields, anchor
// packing and L-systems — plus the CapsuleField and LSystem host classes.
//
// Argument shapes follow docs/mesh-api.js (the JS surface the plant demos
// were written against):
//   Vec3 list:      Float32Array(3N) | Float64Array(3N) | [[x,y,z]|{x,y,z}|n...]
//   Vec2 list:      Float32Array(2N) | [[x,y] ...]
//   BranchSegment:  { parent, from:[x,y,z], to:[x,y,z], radius, depth }
//   Capsule:        { a:[x,y,z], b:[x,y,z], radius, tag? }
//   Sphere:         { center:[x,y,z], radius, tag? }
//   Module:         { symbol:'F', params:[...] }

#include "host_mesh_internal.h"

#include <bromesh/manipulation/bezier_sweep.h>
#include <bromesh/procedural/branches.h>
#include <bromesh/procedural/leaf_scatter.h>
#include <bromesh/procedural/lsystem.h>
#include <bromesh/procedural/lsystem_turtle.h>
#include <bromesh/procedural/obstacle_field.h>
#include <bromesh/procedural/plants.h>
#include <bromesh/procedural/space_colonization.h>

#include <string_view>

namespace bromesh::api {

HostClass g_capsuleFieldClass;
HostClass g_lsystemClass;

namespace {

constexpr uint32_t kHostCapsuleFieldTag = 0x43415046u; // 'CAPF'
constexpr uint32_t kHostLSystemTag      = 0x4C535953u; // 'LSYS'

struct HostCapsuleField {
    std::unique_ptr<bromesh::CapsuleField> field;
    uint32_t tag = kHostCapsuleFieldTag;
};

struct HostLSystem {
    std::unique_ptr<bromesh::LSystem> ls = std::make_unique<bromesh::LSystem>();
    std::vector<bromesh::Module> axiom;
    uint32_t tag = kHostLSystemTag;
};

HostCapsuleField* unwrapCapsuleField(Value v) {
    if (!ev::isObject(v)) return nullptr;
    void* ptr = g_capsuleFieldClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostCapsuleField*>(ptr);
    return (h && h->tag == kHostCapsuleFieldTag) ? h : nullptr;
}

HostLSystem* unwrapLSystem(Value v) {
    if (!ev::isObject(v)) return nullptr;
    void* ptr = g_lsystemClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostLSystem*>(ptr);
    return (h && h->tag == kHostLSystemTag) ? h : nullptr;
}

// ---- option-object readers -------------------------------------------------

bool isPresent(Value v) { return !ev::isUndefined(v) && !ev::isNull(v); }

// An array-like: any object carrying a numeric `length`. bronze's embed API
// has no Array.isArray; every list reader here accepts anything array-like.
bool arrayLength(Value v, size_t& n) {
    n = 0;
    if (!ev::isObject(v)) return false;
    Value lenVal = ev::getProperty(v, "length");
    if (!ev::isNumber(lenVal)) return false;
    double d = ev::toDouble(lenVal);
    if (!(d >= 0.0)) return false;
    n = static_cast<size_t>(d);
    return true;
}

double objNum(Value o, const char* key, double def) {
    if (!ev::isObject(o)) return def;
    Value v = ev::getProperty(o, key);
    if (!ev::isNumber(v)) return def;
    double d = ev::toDouble(v);
    return std::isnan(d) ? def : d;
}

int objInt(Value o, const char* key, int def) {
    return static_cast<int>(objNum(o, key, static_cast<double>(def)));
}

bool objBool(Value o, const char* key, bool def) {
    if (!ev::isObject(o)) return def;
    Value v = ev::getProperty(o, key);
    if (!isPresent(v)) return def;
    return ev::toBool(v);
}

// [x,y,z] or {x,y,z}; anything else reads as the zero vector.
bromath::Vec3 readVec3(Value in) {
    bromath::Vec3 r{};
    if (!ev::isObject(in)) return r;
    Rooted v(in);  // every read below may allocate
    Value x = ev::getProperty(v, "x");
    if (ev::isNumber(x)) {
        r.x = static_cast<float>(ev::toDouble(x));
        r.y = static_cast<float>(objNum(v, "y", 0.0));
        r.z = static_cast<float>(objNum(v, "z", 0.0));
        return r;
    }
    // Numbers are immediates: converting each read at once holds nothing.
    if (Value e = ev::getElement(v, 0); ev::isNumber(e)) r.x = static_cast<float>(ev::toDouble(e));
    if (Value e = ev::getElement(v, 1); ev::isNumber(e)) r.y = static_cast<float>(ev::toDouble(e));
    if (Value e = ev::getElement(v, 2); ev::isNumber(e)) r.z = static_cast<float>(ev::toDouble(e));
    return r;
}

bool objVec3(Value o, const char* key, bromath::Vec3& out) {
    if (!ev::isObject(o)) return false;
    Value v = ev::getProperty(o, key);
    if (!ev::isObject(v)) return false;
    out = readVec3(v);
    return true;
}

// Typed-array fast path shared by the Vec2/Vec3 list readers: a Float32Array
// or Float64Array whose length is a multiple of `stride`.
bool readTypedFloats(Value v, size_t stride, std::vector<float>& out) {
    ev::TypedArrayInfo info = ev::typedArrayInfo(v);
    if (!info) return false;
    if (info.elementCount % stride != 0) return false;
    out.resize(info.elementCount);
    if (info.elementKind == ev::elements::Float32) {
        std::memcpy(out.data(), info.data, info.elementCount * sizeof(float));
        return true;
    }
    if (info.elementKind == ev::elements::Float64) {
        const double* d = reinterpret_cast<const double*>(info.data);
        for (size_t i = 0; i < info.elementCount; ++i) out[i] = static_cast<float>(d[i]);
        return true;
    }
    return false;
}

bool readVec3List(Value in, std::vector<bromath::Vec3>& out) {
    out.clear();
    std::vector<float> flat;
    if (readTypedFloats(in, 3, flat)) {
        out.resize(flat.size() / 3);
        for (size_t i = 0; i < out.size(); ++i) out[i] = {flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]};
        return true;
    }
    Rooted v(in);  // the length and element reads allocate
    size_t n = 0;
    if (!arrayLength(v, n)) return false;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        Value e = ev::getElement(v, static_cast<uint32_t>(i));
        if (ev::isNumber(e)) {
            // Flat [x,y,z,x,y,z,...]
            if (n % 3 != 0) return false;
            out.resize(n / 3);
            for (size_t k = 0; k < out.size(); ++k) {
                out[k].x = static_cast<float>(ev::toDouble(ev::getElement(v, static_cast<uint32_t>(3 * k))));
                out[k].y = static_cast<float>(ev::toDouble(ev::getElement(v, static_cast<uint32_t>(3 * k + 1))));
                out[k].z = static_cast<float>(ev::toDouble(ev::getElement(v, static_cast<uint32_t>(3 * k + 2))));
            }
            return true;
        }
        if (!ev::isObject(e)) return false;
        out.push_back(readVec3(e));
    }
    return true;
}

bool readVec2List(Value in, std::vector<bromath::Vec2>& out) {
    out.clear();
    std::vector<float> flat;
    if (readTypedFloats(in, 2, flat)) {
        out.resize(flat.size() / 2);
        for (size_t i = 0; i < out.size(); ++i) out[i] = {flat[2 * i], flat[2 * i + 1]};
        return true;
    }
    Rooted v(in);  // the length and element reads allocate
    size_t n = 0;
    if (!arrayLength(v, n)) return false;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        Rooted e(ev::getElement(v, static_cast<uint32_t>(i)));
        if (ev::isNumber(e)) {
            if (n % 2 != 0) return false;
            out.resize(n / 2);
            for (size_t k = 0; k < out.size(); ++k) {
                out[k].x = static_cast<float>(ev::toDouble(ev::getElement(v, static_cast<uint32_t>(2 * k))));
                out[k].y = static_cast<float>(ev::toDouble(ev::getElement(v, static_cast<uint32_t>(2 * k + 1))));
            }
            return true;
        }
        if (!ev::isObject(e)) return false;
        bromath::Vec2 p{};
        Value x = ev::getProperty(e, "x");
        if (ev::isNumber(x)) {
            p.x = static_cast<float>(ev::toDouble(x));
            p.y = static_cast<float>(objNum(e, "y", 0.0));
        } else {
            if (Value c = ev::getElement(e, 0); ev::isNumber(c)) p.x = static_cast<float>(ev::toDouble(c));
            if (Value c = ev::getElement(e, 1); ev::isNumber(c)) p.y = static_cast<float>(ev::toDouble(c));
        }
        out.push_back(p);
    }
    return true;
}

// A number, a typed array, or an array of numbers → float list. Used for the
// per-ring profileScale / twist / radii options.
bool readFloatLike(Value in, std::vector<float>& out) {
    out.clear();
    if (ev::isNumber(in)) {
        out.push_back(static_cast<float>(ev::toDouble(in)));
        return true;
    }
    if (!ev::isObject(in)) return false;
    Rooted v(in);  // arrayLength allocates before toFloatVector reads v
    ev::TypedArrayInfo info = ev::typedArrayInfo(v);
    size_t n = 0;
    if (!info && !arrayLength(v, n)) return false;
    out = toFloatVector(v);
    return true;
}

void readFloatLikeOpt(Value o, const char* key, std::vector<float>& out) {
    if (!ev::isObject(o)) return;
    Value v = ev::getProperty(o, key);
    if (!isPresent(v)) return;
    std::vector<float> tmp;
    if (readFloatLike(v, tmp)) out = std::move(tmp);
}

bool readBranchSegments(Value in, std::vector<bromesh::BranchSegment>& out) {
    out.clear();
    Rooted v(in);  // the list and each element are rooted across the reads
    size_t n = 0;
    if (!arrayLength(v, n)) return false;
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        Rooted o(ev::getElement(v, static_cast<uint32_t>(i)));
        bromesh::BranchSegment s{};
        s.parent = objInt(o, "parent", -1);
        s.depth  = objInt(o, "depth", 0);
        s.radius = static_cast<float>(objNum(o, "radius", 0.0));
        objVec3(o, "from", s.from);
        objVec3(o, "to", s.to);
        out[i] = s;
    }
    return true;
}

Value makeVec3Array(const bromath::Vec3& p) {
    return hostArrayOf(3, [&](size_t i) {
        return ev::fromDouble(i == 0 ? p.x : (i == 1 ? p.y : p.z));
    });
}

Value makeBranchSegments(const std::vector<bromesh::BranchSegment>& segs) {
    return hostArrayOf(segs.size(), [&](size_t i) -> Value {
        const auto& s = segs[i];
        ObjectBuilder o;
        o.set("parent", static_cast<double>(s.parent));
        o.set("from", makeVec3Array(s.from));
        o.set("to", makeVec3Array(s.to));
        o.set("radius", static_cast<double>(s.radius));
        o.set("depth", static_cast<double>(s.depth));
        return o.build();
    });
}

bool readCapsules(Value in, std::vector<bromesh::Capsule>& out) {
    out.clear();
    Rooted v(in);  // the list and each element are rooted across the reads
    size_t n = 0;
    if (!arrayLength(v, n)) return false;
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        Rooted o(ev::getElement(v, static_cast<uint32_t>(i)));
        bromesh::Capsule c{};
        objVec3(o, "a", c.a);
        objVec3(o, "b", c.b);
        c.radius = static_cast<float>(objNum(o, "radius", 0.0));
        c.tag    = objInt(o, "tag", -1);
        out[i] = c;
    }
    return true;
}

bool readSpheres(Value in, std::vector<bromesh::Sphere>& out) {
    out.clear();
    Rooted v(in);  // the list and each element are rooted across the reads
    size_t n = 0;
    if (!arrayLength(v, n)) return false;
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        Rooted o(ev::getElement(v, static_cast<uint32_t>(i)));
        bromesh::Sphere s{};
        objVec3(o, "center", s.center);
        s.radius = static_cast<float>(objNum(o, "radius", 0.0));
        s.tag    = objInt(o, "tag", -1);
        out[i] = s;
    }
    return true;
}

void readSpheresOpt(Value o, const char* key, std::vector<bromesh::Sphere>& out) {
    if (!ev::isObject(o)) return;
    Value v = ev::getProperty(o, key);
    if (!isPresent(v)) return;
    readSpheres(v, out);
}

bool readModules(Value in, std::vector<bromesh::Module>& out) {
    out.clear();
    Rooted v(in);  // the list and each element are rooted across the reads
    size_t n = 0;
    if (!arrayLength(v, n)) return false;
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        Rooted o(ev::getElement(v, static_cast<uint32_t>(i)));
        bromesh::Module m{};
        if (ev::isObject(o)) {
            Value sv = ev::getProperty(o, "symbol");
            if (ev::isString(sv)) {
                std::string s = ev::toUtf8(sv);
                if (!s.empty()) m.symbol = s[0];
            }
            Value pv = ev::getProperty(o, "params");
            if (ev::isObject(pv)) m.params = toFloatVector(pv);
        }
        out[i] = std::move(m);
    }
    return true;
}

Value makeModules(const std::vector<bromesh::Module>& mods) {
    return hostArrayOf(mods.size(), [&](size_t i) -> Value {
        const auto& m = mods[i];
        ObjectBuilder o;
        const char sym[2] = {m.symbol, 0};
        o.set("symbol", sym);
        o.set("params", hostArrayOf(m.params.size(), [&](size_t j) {
            return ev::fromDouble(m.params[j]);
        }));
        return o.build();
    });
}

Value makeInt32Array(const int32_t* data, size_t count) {
    Value arr = ev::createTypedArray(ev::elements::Int32, static_cast<uint32_t>(count));
    if (data && count > 0) {
        ev::fillTypedArray(arr, std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(data), count * sizeof(int32_t)));
    }
    return arr;
}

// The CapsuleField behind an option property. Its lifetime is the JS
// wrapper's, which the caller's option object keeps alive across the call.
const bromesh::CapsuleField* readAvoidField(Value opts, const char* key) {
    if (!ev::isObject(opts)) return nullptr;
    Value v = ev::getProperty(opts, key);
    auto* h = unwrapCapsuleField(v);
    return (h && h->field) ? h->field.get() : nullptr;
}

bromesh::LeafShape parseLeafShape(Value v) {
    if (ev::isNumber(v)) {
        int i = static_cast<int>(ev::toDouble(v));
        if (i >= 0 && i <= 5) return static_cast<bromesh::LeafShape>(i);
        return bromesh::LeafShape::Oval;
    }
    if (!ev::isString(v)) return bromesh::LeafShape::Oval;
    std::string s = ev::toUtf8(v);
    if (s == "pointed") return bromesh::LeafShape::Pointed;
    if (s == "lobed")   return bromesh::LeafShape::Lobed;
    if (s == "needle")  return bromesh::LeafShape::Needle;
    if (s == "frond")   return bromesh::LeafShape::Frond;
    if (s == "petal")   return bromesh::LeafShape::Petal;
    return bromesh::LeafShape::Oval;
}

void readLeafPlacementOptions(Value in, bromesh::LeafPlacementOptions& opts) {
    if (!ev::isObject(in)) return;
    Rooted o(in);  // re-read at every field: each read allocates
    opts.maxRadius      = static_cast<float>(objNum(o, "maxRadius",      opts.maxRadius));
    opts.minDepth       = objInt(o, "minDepth", opts.minDepth);
    opts.terminalOnly   = objBool(o, "terminalOnly", opts.terminalOnly);
    opts.perUnitLength  = static_cast<float>(objNum(o, "perUnitLength",  opts.perUnitLength));
    opts.densityFalloff = static_cast<float>(objNum(o, "densityFalloff", opts.densityFalloff));
    opts.upBias         = static_cast<float>(objNum(o, "upBias",         opts.upBias));
    opts.tiltJitter     = static_cast<float>(objNum(o, "tiltJitter",     opts.tiltJitter));
    opts.rollJitter     = static_cast<float>(objNum(o, "rollJitter",     opts.rollJitter));
    opts.baseScale      = static_cast<float>(objNum(o, "baseScale",      opts.baseScale));
    opts.scaleJitter    = static_cast<float>(objNum(o, "scaleJitter",    opts.scaleJitter));
    opts.scaleByRadius  = static_cast<float>(objNum(o, "scaleByRadius",  opts.scaleByRadius));
    opts.dedupRadius    = static_cast<float>(objNum(o, "dedupRadius",    opts.dedupRadius));
    opts.seed           = static_cast<uint64_t>(objNum(o, "seed", static_cast<double>(opts.seed)));
    readFloatLikeOpt(o, "densityWeight", opts.densityWeight);
    opts.avoid             = readAvoidField(o, "avoid");
    opts.obstacleClearance = static_cast<float>(objNum(o, "obstacleClearance", opts.obstacleClearance));
    opts.obstaclePushout   = static_cast<float>(objNum(o, "obstaclePushout",   opts.obstaclePushout));
    readSpheresOpt(o, "keepOut", opts.keepOut);
}

void readColonizeOptions(Value in, bromesh::SpaceColonizationOptions& opts) {
    if (!ev::isObject(in)) return;
    Rooted o(in);  // re-read at every field: each read allocates
    opts.attractionRadius = static_cast<float>(objNum(o, "attractionRadius", opts.attractionRadius));
    opts.killRadius       = static_cast<float>(objNum(o, "killRadius",       opts.killRadius));
    opts.segmentLength    = static_cast<float>(objNum(o, "segmentLength",    opts.segmentLength));
    opts.maxIterations    = objInt(o, "maxIterations", opts.maxIterations);
    opts.tropismWeight    = static_cast<float>(objNum(o, "tropismWeight",    opts.tropismWeight));
    objVec3(o, "tropism", opts.tropism);
    opts.obstacles         = readAvoidField(o, "obstacles");
    opts.obstacleClearance = static_cast<float>(objNum(o, "obstacleClearance", opts.obstacleClearance));
    opts.obstacleSteer     = static_cast<float>(objNum(o, "obstacleSteer",     opts.obstacleSteer));
}

Value makeCapsuleField(std::vector<bromesh::Capsule> caps, std::vector<bromesh::Sphere> sphs, float cellSize) {
    auto h = std::make_unique<HostCapsuleField>();
    h->field = std::make_unique<bromesh::CapsuleField>(std::move(caps), std::move(sphs), cellSize);
    return g_capsuleFieldClass.createInstance(std::move(h));
}

Value capsuleFieldFromArgs(std::span<const Value> a) {
    std::vector<bromesh::Capsule> caps;
    std::vector<bromesh::Sphere>  sphs;
    if (a.size() >= 1 && ev::isObject(a[0])) readCapsules(a[0], caps);
    if (a.size() >= 2 && ev::isObject(a[1])) readSpheres(a[1], sphs);
    float cellSize = 0.0f;
    if (a.size() >= 3 && ev::isNumber(a[2])) cellSize = static_cast<float>(ev::toDouble(a[2]));
    return makeCapsuleField(std::move(caps), std::move(sphs), cellSize);
}

int excludeTagAt(std::span<const Value> a, size_t i) {
    return (i < a.size() && ev::isNumber(a[i])) ? static_cast<int>(ev::toDouble(a[i])) : -1;
}

} // namespace

// ---------------------------------------------------------------------------
// Mesh statics
// ---------------------------------------------------------------------------
void initMeshPlants(HostClass& cls) {
    auto bindStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        cls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    // ---- Sweeps ------------------------------------------------------------

    // Mesh.sweep(profile, path, opts?) — extrude a 2D profile along a 3D path.
    bindStatic("sweep", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.sweep requires (profile, path[, opts])");
        std::vector<bromath::Vec2> profile;
        std::vector<bromath::Vec3> path;
        if (!readVec2List(a[0], profile)) return ev::throwTypeError("Mesh.sweep: profile must be Float32Array(2N) or [[x,y],...]");
        if (!readVec3List(a[1], path)) return ev::throwTypeError("Mesh.sweep: path must be Float32Array(3N) or [[x,y,z],...]");
        bromesh::SweepOptions opts;
        if (a.size() > 2 && ev::isObject(a[2])) {
            Rooted o(a[2]);
            opts.closeProfile = objBool(o, "closeProfile", opts.closeProfile);
            opts.capStart     = objBool(o, "capStart",     opts.capStart);
            opts.capEnd       = objBool(o, "capEnd",       opts.capEnd);
            opts.miterJoints  = objBool(o, "miterJoints",  opts.miterJoints);
            readFloatLikeOpt(o, "profileScale", opts.profileScale);
            readFloatLikeOpt(o, "twist", opts.twist);
        }
        return wrapMesh(bromesh::sweep(profile, path, opts));
    });

    // Mesh.bezierSweep(controlPoints, profile, opts?) — cubic Bézier spine.
    bindStatic("bezierSweep", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.bezierSweep requires (controlPoints, profile[, opts])");
        std::vector<bromath::Vec3> ctrl;
        std::vector<bromath::Vec2> profile;
        if (!readVec3List(a[0], ctrl)) return ev::throwTypeError("Mesh.bezierSweep: controlPoints must be a Vec3 list");
        if (!readVec2List(a[1], profile)) return ev::throwTypeError("Mesh.bezierSweep: profile must be a Vec2 list");
        bromesh::BezierSweepOptions o;
        if (a.size() > 2 && ev::isObject(a[2])) {
            Rooted ov(a[2]);
            o.samples      = objInt(ov, "samples", o.samples);
            o.capStart     = objBool(ov, "capStart", o.capStart);
            o.capEnd       = objBool(ov, "capEnd", o.capEnd);
            o.closeProfile = objBool(ov, "closeProfile", o.closeProfile);
            o.miterJoints  = objBool(ov, "miterJoints", o.miterJoints);
            readFloatLikeOpt(ov, "profileScale", o.profileScale);
            readFloatLikeOpt(ov, "twist", o.twist);
        }
        return wrapMesh(bromesh::bezierSweep(ctrl, profile, o));
    });

    // Mesh.tube(path, radius | radii[], sides = 8, opts?) — circular sweep.
    bindStatic("tube", 4, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.tube requires (path[, radius[, sides[, opts]]])");
        std::vector<bromath::Vec3> path;
        if (!readVec3List(a[0], path) || path.size() < 2) return ev::throwTypeError("Mesh.tube: path must be a Vec3 list of at least 2 points");
        std::vector<float> radii;
        if (a.size() > 1 && isPresent(a[1])) {
            if (!readFloatLike(a[1], radii)) return ev::throwTypeError("Mesh.tube: radius must be a number or float list");
        }
        if (radii.empty()) radii.push_back(0.1f);
        bromesh::TubeOptions opts;
        if (a.size() > 2 && ev::isNumber(a[2])) {
            int s = static_cast<int>(ev::toDouble(a[2]));
            if (s >= 3) opts.sides = s;
        }
        if (a.size() > 3 && ev::isObject(a[3])) {
            opts.capStart    = objBool(a[3], "capStart",    opts.capStart);
            opts.capEnd      = objBool(a[3], "capEnd",      opts.capEnd);
            opts.miterJoints = objBool(a[3], "miterJoints", opts.miterJoints);
        }
        return wrapMesh(bromesh::tube(path, radii, opts));
    });

    // ---- Plant cards -------------------------------------------------------

    // Mesh.blob({ radius, seed, nsub, scale, center }) — noise-displaced
    // sphere with non-uniform scale and translation baked in. A leading
    // number reads as (radius, seed, nsub) positionally.
    bindStatic("blob", 1, [](Value, std::span<const Value> a) -> Value {
        double radius = 0.5;
        int seed = 42;
        int nsub = 2;
        bromath::Vec3 scale{1.0f, 1.0f, 1.0f};
        bromath::Vec3 center{0.0f, 0.0f, 0.0f};
        if (!a.empty() && ev::isObject(a[0])) {
            Rooted o(a[0]);
            radius = objNum(o, "radius", radius);
            seed   = objInt(o, "seed", seed);
            nsub   = objInt(o, "nsub", nsub);
            Value sv = ev::getProperty(o, "scale");
            if (ev::isNumber(sv)) {
                float s = static_cast<float>(ev::toDouble(sv));
                scale = {s, s, s};
            } else if (ev::isObject(sv)) {
                scale = readVec3(sv);
            }
            objVec3(o, "center", center);
        } else if (!a.empty() && ev::isNumber(a[0])) {
            ArgReader r(a);
            radius = r.getDouble(0, radius);
            seed   = r.getInt(1, seed);
            nsub   = r.getInt(2, nsub);
        }
        if (nsub < 0) nsub = 0;
        return wrapMesh(bromesh::blob(static_cast<float>(radius), seed, nsub,
                                      scale.x, scale.y, scale.z,
                                      center.x, center.y, center.z));
    });

    // Mesh.leafCard(shape, opts?) — a bent/curled/cupped leaf card.
    bindStatic("leafCard", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.leafCard requires (shape[, opts])");
        bromesh::LeafShape shape = parseLeafShape(a[0]);
        bromesh::LeafCardOptions o;
        if (a.size() > 1 && ev::isObject(a[1])) {
            Rooted ov(a[1]);
            o.width  = static_cast<float>(objNum(ov, "width",  o.width));
            o.length = static_cast<float>(objNum(ov, "length", o.length));
            o.bend   = static_cast<float>(objNum(ov, "bend",   o.bend));
            o.curl   = static_cast<float>(objNum(ov, "curl",   o.curl));
            o.cup    = static_cast<float>(objNum(ov, "cup",    o.cup));
            o.stemOffset       = objBool(ov, "stemOffset", o.stemOffset);
            o.widthSegments    = objInt(ov, "widthSegments",  o.widthSegments);
            o.lengthSegments   = objInt(ov, "lengthSegments", o.lengthSegments);
            o.fullUV           = objBool(ov, "fullUV", o.fullUV);
            o.shapedSilhouette = objBool(ov, "shapedSilhouette", o.shapedSilhouette);
        }
        return wrapMesh(bromesh::leafCard(shape, o));
    });

    // Mesh.flower(opts?) — petal rings around a center disc.
    bindStatic("flower", 1, [](Value, std::span<const Value> a) -> Value {
        bromesh::FlowerOptions o;
        if (!a.empty() && ev::isObject(a[0])) {
            Rooted ov(a[0]);
            o.petalCount = objInt(ov, "petalCount", o.petalCount);
            Value ps = ev::getProperty(ov, "petalShape");
            if (isPresent(ps)) o.petalShape = parseLeafShape(ps);
            o.petalLength  = static_cast<float>(objNum(ov, "petalLength",  o.petalLength));
            o.petalWidth   = static_cast<float>(objNum(ov, "petalWidth",   o.petalWidth));
            o.petalCurl    = static_cast<float>(objNum(ov, "petalCurl",    o.petalCurl));
            o.petalBend    = static_cast<float>(objNum(ov, "petalBend",    o.petalBend));
            o.layers       = objInt(ov, "layers", o.layers);
            o.layerTwist   = static_cast<float>(objNum(ov, "layerTwist",   o.layerTwist));
            o.centerRadius = static_cast<float>(objNum(ov, "centerRadius", o.centerRadius));
            o.centerHeight = static_cast<float>(objNum(ov, "centerHeight", o.centerHeight));
            o.outerTilt    = static_cast<float>(objNum(ov, "outerTilt",    o.outerTilt));
            o.innerTilt    = static_cast<float>(objNum(ov, "innerTilt",    o.innerTilt));
            o.layerScaleFalloff = static_cast<float>(objNum(ov, "layerScaleFalloff", o.layerScaleFalloff));
            o.outerYLift   = static_cast<float>(objNum(ov, "outerYLift",   o.outerYLift));
            o.innerYLift   = static_cast<float>(objNum(ov, "innerYLift",   o.innerYLift));
            o.petalCup     = static_cast<float>(objNum(ov, "petalCup",     o.petalCup));
            o.shapedPetals = objBool(ov, "shapedPetals", o.shapedPetals);
            Value cc = ev::getProperty(ov, "centerColor");
            if (ev::isObject(cc)) {
                std::vector<float> c = toFloatVector(cc);
                for (size_t i = 0; i < 3 && i < c.size(); ++i) o.centerColor[i] = c[i];
            }
        }
        return wrapMesh(bromesh::flower(o));
    });

    // Mesh.bladeStrip(path, opts?) — diamond-profile sweep for grass blades.
    bindStatic("bladeStrip", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.bladeStrip requires (path[, opts])");
        std::vector<bromath::Vec3> path;
        if (!readVec3List(a[0], path)) return ev::throwTypeError("Mesh.bladeStrip: path must be a Vec3 list");
        bromesh::BladeStripOptions o;
        if (a.size() > 1 && ev::isObject(a[1])) {
            Rooted ov(a[1]);
            o.width       = static_cast<float>(objNum(ov, "width",     o.width));
            o.thickness   = static_cast<float>(objNum(ov, "thickness", o.thickness));
            o.capStart    = objBool(ov, "capStart",    o.capStart);
            o.capEnd      = objBool(ov, "capEnd",      o.capEnd);
            o.miterJoints = objBool(ov, "miterJoints", o.miterJoints);
            readFloatLikeOpt(ov, "profileScale", o.profileScale);
            readFloatLikeOpt(ov, "twist", o.twist);
        }
        return wrapMesh(bromesh::bladeStrip(path, o));
    });

    // Mesh.bladePath(opts?) → [[x,y,z], ...] quadratic-Bézier blade spine.
    bindStatic("bladePath", 1, [](Value, std::span<const Value> a) -> Value {
        bromesh::BladePathOptions o;
        if (!a.empty() && ev::isObject(a[0])) {
            Rooted ov(a[0]);
            objVec3(ov, "base", o.base);
            objVec3(ov, "tipDir", o.tipDir);
            o.length   = static_cast<float>(objNum(ov, "length", o.length));
            o.bend     = static_cast<float>(objNum(ov, "bend",   o.bend));
            o.lift     = static_cast<float>(objNum(ov, "lift",   o.lift));
            o.segments = objInt(ov, "segments", o.segments);
        }
        std::vector<bromath::Vec3> pts = bromesh::bladePath(o);
        return hostArrayOf(pts.size(), [&](size_t i) { return makeVec3Array(pts[i]); });
    });

    // ---- Branch trees ------------------------------------------------------

    bindStatic("spaceColonize", 4, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("Mesh.spaceColonize requires (attractors, seedPoints, initialDirection[, opts])");
        std::vector<bromath::Vec3> attractors, seeds;
        if (!readVec3List(a[0], attractors)) return ev::throwTypeError("Mesh.spaceColonize: attractors must be a Vec3 list");
        if (!readVec3List(a[1], seeds)) return ev::throwTypeError("Mesh.spaceColonize: seedPoints must be a Vec3 list");
        bromath::Vec3 initDir = readVec3(a[2]);
        bromesh::SpaceColonizationOptions opts;
        if (a.size() > 3) readColonizeOptions(a[3], opts);
        return makeBranchSegments(bromesh::spaceColonize(attractors, seeds, initDir, opts));
    });

    bindStatic("thickenBranches", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.thickenBranches requires (segments[, leafRadius[, pipeExp]])");
        std::vector<bromesh::BranchSegment> segs;
        if (!readBranchSegments(a[0], segs)) return ev::throwTypeError("Mesh.thickenBranches: segments must be an array of branch segment objects");
        ArgReader r(a);
        float leafR   = static_cast<float>(r.getDouble(1, 0.02));
        float pipeExp = static_cast<float>(r.getDouble(2, 2.5));
        bromesh::thickenBranches(segs, leafR, pipeExp);
        return makeBranchSegments(segs);
    });

    bindStatic("meshBranches", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.meshBranches requires (segments[, sides])");
        std::vector<bromesh::BranchSegment> segs;
        if (!readBranchSegments(a[0], segs)) return ev::throwTypeError("Mesh.meshBranches: segments must be an array of branch segment objects");
        int sides = 8;
        if (a.size() > 1 && ev::isNumber(a[1])) {
            int s = static_cast<int>(ev::toDouble(a[1]));
            if (s >= 3) sides = s;
        }
        return wrapMesh(bromesh::meshBranches(segs, sides));
    });

    bindStatic("placeLeavesOnBranches", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.placeLeavesOnBranches requires (segments[, opts])");
        std::vector<bromesh::BranchSegment> segs;
        if (!readBranchSegments(a[0], segs)) return ev::throwTypeError("Mesh.placeLeavesOnBranches: segments must be an array of branch segment objects");
        bromesh::LeafPlacementOptions opts;
        if (a.size() > 1) readLeafPlacementOptions(a[1], opts);
        bromesh::LeafPlacements p = bromesh::placeLeavesOnBranches(segs, opts);
        ObjectBuilder obj;
        obj.set("count", static_cast<double>(p.count()));
        obj.set("transforms", makeFloat32Array(p.transforms.data(), p.transforms.size()));
        obj.set("branchRadius", makeFloat32Array(p.branchRadius.data(), p.branchRadius.size()));
        obj.set("branchDepth", makeInt32Array(p.branchDepth.data(), p.branchDepth.size()));
        return obj.build();
    });

    bindStatic("scatterLeaves", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Mesh.scatterLeaves requires (segments, leaf[, opts])");
        std::vector<bromesh::BranchSegment> segs;
        if (!readBranchSegments(a[0], segs)) return ev::throwTypeError("Mesh.scatterLeaves: segments must be an array of branch segment objects");
        auto* leaf = unwrapMesh(a[1]);
        if (!leaf) return ev::throwTypeError("Mesh.scatterLeaves: leaf must be a Mesh");
        bromesh::LeafPlacementOptions opts;
        if (a.size() > 2) readLeafPlacementOptions(a[2], opts);
        return wrapMesh(bromesh::scatterLeaves(segs, leaf->mesh, opts));
    });

    // Mesh.tree(opts?) → { segments, branches }: spaceColonize →
    // thickenBranches → meshBranches in one call.
    bindStatic("tree", 1, [](Value, std::span<const Value> a) -> Value {
        bromesh::TreeOptions o;
        if (!a.empty() && ev::isObject(a[0])) {
            Rooted ov(a[0]);
            objVec3(ov, "base", o.base);
            objVec3(ov, "canopyCenter", o.canopyCenter);
            o.canopyRadius   = static_cast<float>(objNum(ov, "canopyRadius", o.canopyRadius));
            o.attractorCount = objInt(ov, "attractorCount", o.attractorCount);
            o.sides          = objInt(ov, "sides", o.sides);
            o.leafRadius     = static_cast<float>(objNum(ov, "leafRadius", o.leafRadius));
            o.pipeExp        = static_cast<float>(objNum(ov, "pipeExp", o.pipeExp));
            o.seed           = objInt(ov, "seed", o.seed);
            Value colv = ev::getProperty(ov, "colonize");
            if (ev::isObject(colv)) readColonizeOptions(colv, o.colonize);
        }
        bromesh::TreeResult r = bromesh::tree(o);
        ObjectBuilder out;
        out.set("segments", makeBranchSegments(r.segments));
        out.set("branches", wrapMesh(std::move(r.branches)));
        return out.build();
    });

    // ---- Collision-aware placement ----------------------------------------

    bindStatic("capsuleField", 3, [](Value, std::span<const Value> a) -> Value {
        return capsuleFieldFromArgs(a);
    });

    bindStatic("capsuleFieldFromSegments", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.capsuleFieldFromSegments requires (segments[, radiusScale[, spheres]])");
        std::vector<bromesh::BranchSegment> segs;
        if (!readBranchSegments(a[0], segs)) return ev::throwTypeError("Mesh.capsuleFieldFromSegments: segments must be a branch-segment array");
        float radiusScale = 1.0f;
        if (a.size() >= 2 && ev::isNumber(a[1])) radiusScale = static_cast<float>(ev::toDouble(a[1]));
        std::vector<bromesh::Sphere> sphs;
        if (a.size() >= 3 && ev::isObject(a[2])) readSpheres(a[2], sphs);
        return makeCapsuleField(bromesh::CapsuleField::capsulesFromSegments(segs, radiusScale), std::move(sphs), 0.0f);
    });

    // Mesh.packAnchors(candidates, opts?) → Int32Array of accepted indices.
    bindStatic("packAnchors", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.packAnchors requires (candidates[, opts])");
        std::vector<bromath::Vec3> cand;
        if (!readVec3List(a[0], cand)) return ev::throwTypeError("Mesh.packAnchors: candidates must be a Vec3 list");
        bromesh::AnchorPackOptions opts;
        const bromesh::CapsuleField* avoid = nullptr;
        std::vector<bromesh::Sphere> keepOut;
        if (a.size() > 1 && ev::isObject(a[1])) {
            Rooted ov(a[1]);
            opts.minSpacing          = static_cast<float>(objNum(ov, "minSpacing",          opts.minSpacing));
            opts.minObstacleDistance = static_cast<float>(objNum(ov, "minObstacleDistance", opts.minObstacleDistance));
            opts.maxCount            = objInt(ov, "maxCount", opts.maxCount);
            opts.seed                = static_cast<uint64_t>(objNum(ov, "seed", static_cast<double>(opts.seed)));
            avoid = readAvoidField(ov, "avoid");
            readSpheresOpt(ov, "keepOut", keepOut);
        }
        std::vector<int> idx = bromesh::packAnchors(cand, avoid, keepOut, opts);
        static_assert(sizeof(int) == sizeof(int32_t));
        return makeInt32Array(reinterpret_cast<const int32_t*>(idx.data()), idx.size());
    });

    // ---- L-systems ---------------------------------------------------------

    bindStatic("parseLSystem", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isString(a[0])) return ev::throwTypeError("Mesh.parseLSystem requires (text)");
        std::string s = ev::toUtf8(a[0]);
        return makeModules(bromesh::parseModules(std::string_view(s)));
    });

    bindStatic("lsystemToBranches", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.lsystemToBranches requires (modules[, opts])");
        std::vector<bromesh::Module> mods;
        if (ev::isString(a[0])) {
            mods = bromesh::parseModules(ev::toUtf8(a[0]));
        } else if (!readModules(a[0], mods)) {
            return ev::throwTypeError("Mesh.lsystemToBranches: modules must be an array of {symbol, params}");
        }
        bromesh::TurtleOptions to;
        if (a.size() > 1 && ev::isObject(a[1])) {
            Rooted ov(a[1]);
            to.stepLength = static_cast<float>(objNum(ov, "stepLength", to.stepLength));
            to.angle      = static_cast<float>(objNum(ov, "angle",      to.angle));
            to.radius     = static_cast<float>(objNum(ov, "radius",     to.radius));
            objVec3(ov, "position", to.position);
            objVec3(ov, "heading", to.heading);
            objVec3(ov, "up", to.up);
        }
        return makeBranchSegments(bromesh::lsystemToBranches(mods, to));
    });
}

// ---------------------------------------------------------------------------
// CapsuleField class — capsule + sphere occupancy field for placement tests
// ---------------------------------------------------------------------------
void initCapsuleField(HostClass& cls) {
    cls.install("CapsuleField", 3, [](Value, std::span<const Value> a) -> Value {
        return capsuleFieldFromArgs(a);
    }, [](ObjectBuilder& proto) {
        proto.accessor("empty", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapCapsuleField(self);
            return ev::fromBool(!h || !h->field || h->field->empty());
        });
        proto.accessor("capsuleCount", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapCapsuleField(self);
            return ev::fromDouble(h && h->field ? static_cast<double>(h->field->capsuleCount()) : 0.0);
        });
        proto.accessor("sphereCount", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapCapsuleField(self);
            return ev::fromDouble(h && h->field ? static_cast<double>(h->field->sphereCount()) : 0.0);
        });
        proto.accessor("cellSize", [](Value self, std::span<const Value>) -> Value {
            auto* h = unwrapCapsuleField(self);
            return ev::fromDouble(h && h->field ? static_cast<double>(h->field->cellSize()) : 0.0);
        });
        proto.def("contains", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapCapsuleField(self);
            if (!h || !h->field || a.empty()) return ev::fromBool(false);
            float extra = (a.size() > 2 && ev::isNumber(a[2])) ? static_cast<float>(ev::toDouble(a[2])) : 0.0f;
            return ev::fromBool(h->field->contains(readVec3(a[0]), excludeTagAt(a, 1), extra));
        });
        proto.def("tooClose", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapCapsuleField(self);
            if (!h || !h->field || a.empty()) return ev::fromBool(false);
            float clearance = (a.size() > 1 && ev::isNumber(a[1])) ? static_cast<float>(ev::toDouble(a[1])) : 0.0f;
            return ev::fromBool(h->field->tooClose(readVec3(a[0]), clearance, excludeTagAt(a, 2)));
        });
        proto.def("distance", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapCapsuleField(self);
            if (!h || !h->field || a.empty()) return ev::fromDouble(0.0);
            return ev::fromDouble(h->field->distance(readVec3(a[0]), excludeTagAt(a, 1)));
        });
        proto.def("nearest", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapCapsuleField(self);
            if (!h || !h->field || a.empty()) return ev::null();
            auto n = h->field->nearest(readVec3(a[0]), excludeTagAt(a, 1));
            ObjectBuilder o;
            o.set("point", makeVec3Array(n.point));
            o.set("normal", makeVec3Array(n.normal));
            o.set("distance", static_cast<double>(n.distance));
            o.set("tag", static_cast<double>(n.tag));
            return o.build();
        });
        proto.def("intersectsSphere", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapCapsuleField(self);
            if (!h || !h->field || a.size() < 2) return ev::fromBool(false);
            float radius = ev::isNumber(a[1]) ? static_cast<float>(ev::toDouble(a[1])) : 0.0f;
            return ev::fromBool(h->field->intersectsSphere(readVec3(a[0]), radius, excludeTagAt(a, 2)));
        });
    });
}

// ---------------------------------------------------------------------------
// LSystem class — string-rule stochastic L-systems. Parametric rules with
// conditions stay native-only; parseLSystem / deriveModules expose the
// module stream for callers who interpret it themselves.
// ---------------------------------------------------------------------------
void initLSystem(HostClass& cls) {
    cls.install("LSystem", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostLSystem>();
        if (!a.empty() && ev::isString(a[0])) {
            h->axiom = bromesh::parseModules(ev::toUtf8(a[0]));
            h->ls->setAxiom(h->axiom);
        }
        return g_lsystemClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.def("setAxiom", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapLSystem(self);
            if (!h) return ev::throwTypeError("LSystem.setAxiom: not an LSystem instance");
            h->axiom = bromesh::parseModules(strAt(a, 0));
            h->ls->setAxiom(h->axiom);
            return self;
        });
        proto.def("addRule", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapLSystem(self);
            if (!h) return ev::throwTypeError("LSystem.addRule: not an LSystem instance");
            std::string predecessor = strAt(a, 0);
            if (predecessor.empty()) return ev::throwTypeError("LSystem.addRule: predecessor must be a one-character symbol");
            bromesh::ProductionRule rule;
            rule.predecessor = predecessor[0];
            rule.weight = static_cast<float>(hasArg(a, 2) ? numAt(a, 2) : 1.0);
            auto mods = bromesh::parseModules(strAt(a, 1));
            rule.successor = [mods](const std::vector<float>&) { return mods; };
            h->ls->addRule(std::move(rule));
            return self;
        });
        proto.def("derive", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapLSystem(self);
            if (!h) return ev::throwTypeError("LSystem.derive: not an LSystem instance");
            auto mods = h->ls->derive(i32At(a, 0), static_cast<uint64_t>(i64At(a, 1)));
            return ev::fromUtf8(bromesh::serializeModules(mods));
        });
        proto.def("deriveModules", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* h = unwrapLSystem(self);
            if (!h) return ev::throwTypeError("LSystem.deriveModules: not an LSystem instance");
            return makeModules(h->ls->derive(i32At(a, 0), static_cast<uint64_t>(i64At(a, 1))));
        });
    });
}

} // namespace bromesh::api
