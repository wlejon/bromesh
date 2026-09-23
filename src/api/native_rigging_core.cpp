#include "host_mesh_internal.h"
#include <fstream>
#include <sstream>

namespace bromesh::api {

HostClass g_skinDataClass;
HostClass g_skeletonClass;
HostClass g_jointClass;
HostClass g_skeletonRigClass;
HostClass g_voxelChunkClass;

namespace {

Value makeSkinValidationObject(const bromesh::SkinValidation& sv) {
    ObjectBuilder obj;
    obj.set("clean", sv.clean());
    obj.set("vertexCount", static_cast<double>(sv.vertexCount));
    obj.set("orphanCount", static_cast<double>(sv.orphanCount));
    obj.set("badSumCount", static_cast<double>(sv.badSumCount));
    obj.set("nanCount", static_cast<double>(sv.nanCount));
    obj.set("maxSumDeviation", static_cast<double>(sv.maxSumDeviation));
    obj.set("maxInfluencesObserved", static_cast<double>(sv.maxInfluencesObserved));
    return obj.build();
}

bromesh::Landmarks landmarksFromObject(Value v) {
    bromesh::Landmarks lm;
    if (!ev::isObject(v)) return lm;
    // Every step allocates (the reads, Object.keys, the conversions), so the
    // source object, the key list and Object.keys itself are all rooted.
    Rooted src(v);
    Rooted root(src.get());
    {
        Value pts = ev::getProperty(src, "points");
        if (ev::isObject(pts)) root.p.set(pts);
    }

    ev::GlobalValue objectCtor = ev::globalValue("Object");
    if (!objectCtor.found || !ev::isObject(objectCtor.value)) return lm;
    Rooted keysFn(ev::getProperty(objectCtor.value, "keys"));
    if (ev::isFunction(keysFn)) {
        const Value args[1] = {root.get()};
        auto res = ev::call(keysFn, ev::undefined(), std::span<const Value>(args, 1));
        if (!res.thrown && ev::isObject(res.value)) {
            Rooted keys(res.value);
            Value lenVal = ev::getProperty(keys, "length");
            size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
            for (size_t i = 0; i < n; ++i) {
                std::string k = ev::toUtf8(ev::getElement(keys, static_cast<uint32_t>(i)));
                std::vector<float> pt = toFloatVector(ev::getProperty(root, k));
                if (pt.size() >= 3) {
                    lm.set(k, pt[0], pt[1], pt[2]);
                }
            }
        }
    }
    return lm;
}

Value landmarksToObject(const bromesh::Landmarks& lm) {
    ObjectBuilder obj;
    ObjectBuilder pts;
    for (const auto& [name, p] : lm.points) {
        const float pt[3] = {p[0], p[1], p[2]};
        ev::Persistent v(makeFloat32Array(pt, 3));
        pts.set(name, v.get());
    }
    obj.set("points", pts.build());
    return obj.build();
}

// Number-valued property setters: absent / non-number keys leave the default.
template <typename T>
void readNum(Value obj, const char* key, T& out) {
    Value v = ev::getProperty(obj, key);
    if (ev::isNumber(v)) out = static_cast<T>(ev::toDouble(v));
}

// Rig.autoRig's options bag: method + the shared smoothing / pruning keys,
// then the per-method blocks { voxel, boneHeat, bbw }. Each block is read
// independently so a caller can set all three and switch `method`.
bromesh::WeightingOptions weightingOptionsFromObject(Value in) {
    bromesh::WeightingOptions wo;
    if (!ev::isObject(in)) return wo;
    Rooted opts(in);  // re-read at every field: each read allocates
    Value mVal = ev::getProperty(opts, "method");
    if (ev::isString(mVal)) wo.method = bromesh::parseWeightingMethod(ev::toUtf8(mVal).c_str());
    readNum(opts, "smoothIterations", wo.smoothIterations);
    readNum(opts, "smoothAlpha", wo.smoothAlpha);
    readNum(opts, "minWeight", wo.minWeight);

    Rooted v(ev::getProperty(opts, "voxel"));
    if (ev::isObject(v)) {
        readNum(v, "maxResolution", wo.voxel.maxResolution);
        readNum(v, "maxInfluences", wo.voxel.maxInfluences);
        readNum(v, "falloffPower", wo.voxel.falloffPower);
        readNum(v, "minWeight", wo.voxel.minWeight);
        readNum(v, "smoothIterations", wo.voxel.smoothIterations);
        readNum(v, "smoothAlpha", wo.voxel.smoothAlpha);
    }
    v.p.set(ev::getProperty(opts, "boneHeat"));
    if (ev::isObject(v)) {
        readNum(v, "maxInfluences", wo.boneHeat.maxInfluences);
        readNum(v, "minWeight", wo.boneHeat.minWeight);
        readNum(v, "heatStrength", wo.boneHeat.heatStrength);
        readNum(v, "solverTol", wo.boneHeat.solverTol);
        readNum(v, "solverMaxIter", wo.boneHeat.solverMaxIter);
    }
    v.p.set(ev::getProperty(opts, "bbw"));
    if (ev::isObject(v)) {
        readNum(v, "maxInfluences", wo.bbw.maxInfluences);
        readNum(v, "minWeight", wo.bbw.minWeight);
        readNum(v, "anchorsPerBone", wo.bbw.anchorsPerBone);
        readNum(v, "eps", wo.bbw.eps);
        readNum(v, "maxIter", wo.bbw.maxIter);
    }
    return wo;
}

// Rig.generateLocomotionCycle's params: a gait name string, or
// { strideLength, cycleDuration, footLiftHeight, keyframesPerCycle,
//   bodyBobAmplitude, armSwingAmplitude, forwardAxis, upAxis,
//   gait: string | { name, phases, dutyFactor } }.
bromesh::LocomotionParams locomotionParamsFromValue(Value in) {
    bromesh::LocomotionParams p;
    if (ev::isString(in)) {
        p.gait.name = ev::toUtf8(in);
        return p;
    }
    if (!ev::isObject(in)) return p;
    Rooted v(in);  // re-read at every field: each read allocates
    readNum(v, "strideLength", p.strideLength);
    readNum(v, "cycleDuration", p.cycleDuration);
    readNum(v, "footLiftHeight", p.footLiftHeight);
    readNum(v, "keyframesPerCycle", p.keyframesPerCycle);
    readNum(v, "bodyBobAmplitude", p.bodyBobAmplitude);
    readNum(v, "armSwingAmplitude", p.armSwingAmplitude);
    for (auto [key, dst] : {std::pair<const char*, float*>{"forwardAxis", p.forwardAxis},
                            std::pair<const char*, float*>{"upAxis", p.upAxis}}) {
        std::vector<float> axis = toFloatVector(ev::getProperty(v, key));
        if (axis.size() >= 3) { dst[0] = axis[0]; dst[1] = axis[1]; dst[2] = axis[2]; }
    }
    Rooted g(ev::getProperty(v, "gait"));
    if (ev::isString(g)) {
        p.gait.name = ev::toUtf8(g);
    } else if (ev::isObject(g)) {
        if (Value n = ev::getProperty(g, "name"); ev::isString(n)) p.gait.name = ev::toUtf8(n);
        Value ph = ev::getProperty(g, "phases");
        if (ev::isObject(ph)) p.gait.phases = toFloatVector(ph);
        readNum(g, "dutyFactor", p.gait.dutyFactor);
    }
    return p;
}

// A bone given as a plain object: name, parent, localT|translation,
// localR|rotation, localS|scale, inverseBind|inverseBindMatrix. Every field is
// read right before it is converted, with the object rooted across them.
void readBoneObject(Value in, bromesh::Bone& b) {
    Rooted o(in);
    auto either = [&](const char* k1, const char* k2) {
        Value v = ev::getProperty(o, k1);
        return ev::isUndefined(v) ? ev::getProperty(o, k2) : v;
    };
    auto floats = [&](const char* k1, const char* k2, float* dst, size_t n) {
        Value v = either(k1, k2);
        if (ev::isUndefined(v)) return;
        std::vector<float> f = toFloatVector(v);
        if (f.size() >= n) std::copy(f.begin(), f.begin() + n, dst);
    };
    if (Value v = ev::getProperty(o, "name"); ev::isString(v)) b.name = ev::toUtf8(v);
    if (Value v = ev::getProperty(o, "parent"); ev::isNumber(v)) b.parent = static_cast<int>(ev::toDouble(v));
    floats("localT", "translation", b.localT, 3);
    floats("localR", "rotation", b.localR, 4);
    floats("localS", "scale", b.localS, 3);
    floats("inverseBind", "inverseBindMatrix", b.inverseBind, 16);
}

} // namespace

void initRiggingCore(HostClass& skinCls, HostClass& skelCls, HostClass& jointCls,
                     HostClass& rigCls, HostClass& voxelCls) {
    // =========================================================================
    // SkinData Class
    // =========================================================================
    skinCls.install("SkinData", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostSkinData>();
        if (!a.empty() && ev::isObject(a[0])) {
            // Each field is read right before it is converted: a read (and a
            // conversion) allocates, which would stale an earlier read.
            Rooted opts(a[0]);
            auto either = [&](const char* k1, const char* k2) {
                Value v = ev::getProperty(opts, k1);
                return ev::isUndefined(v) ? ev::getProperty(opts, k2) : v;
            };
            if (Value v = either("boneWeights", "weights"); !ev::isUndefined(v)) {
                h->skin.boneWeights = toFloatVector(v);
            }
            if (Value v = either("boneIndices", "indices"); !ev::isUndefined(v)) {
                h->skin.boneIndices = toUint32Vector(v);
            }
            if (Value v = ev::getProperty(opts, "inverseBindMatrices"); !ev::isUndefined(v)) {
                h->skin.inverseBindMatrices = toFloatVector(v);
            }
            Value bcVal = ev::getProperty(opts, "boneCount");
            if (ev::isNumber(bcVal)) {
                h->skin.boneCount = static_cast<size_t>(ev::toDouble(bcVal));
            } else if (!h->skin.inverseBindMatrices.empty()) {
                h->skin.boneCount = h->skin.inverseBindMatrices.size() / 16;
            }
        }
        return g_skinDataClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("boneWeights", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::undefined();
            return makeFloat32Array(s->skin.boneWeights.data(), s->skin.boneWeights.size());
        });
        // Alias getters. They used to ev::call the VALUE of boneWeights (a
        // Float32Array, not a function), which threw and handed the error
        // object back as the property's value.
        proto.accessor("weights", [](Value self, std::span<const Value>) -> Value {
            return ev::getProperty(self, "boneWeights");
        });

        proto.accessor("boneIndices", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::undefined();
            std::vector<uint32_t> idx32(s->skin.boneIndices.begin(), s->skin.boneIndices.end());
            return makeUint32Array(idx32.data(), idx32.size());
        });
        proto.accessor("indices", [](Value self, std::span<const Value>) -> Value {
            return ev::getProperty(self, "boneIndices");
        });

        proto.accessor("inverseBindMatrices", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::undefined();
            return makeFloat32Array(s->skin.inverseBindMatrices.data(), s->skin.inverseBindMatrices.size());
        });

        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            return ev::fromDouble(s ? static_cast<double>(s->skin.boneCount) : 0.0);
        });

        proto.accessor("vertexCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            return ev::fromDouble(s ? static_cast<double>(s->skin.boneWeights.size() / 4) : 0.0);
        });

        proto.accessor("maxWeights", [](Value, std::span<const Value>) -> Value {
            return ev::fromDouble(4.0);
        });

        proto.def("normalize", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::throwTypeError("SkinData.normalize: not an instance");
            for (size_t v = 0; v < s->skin.boneWeights.size(); v += 4) {
                float sum = s->skin.boneWeights[v] + s->skin.boneWeights[v+1] +
                            s->skin.boneWeights[v+2] + s->skin.boneWeights[v+3];
                if (sum > 1e-6f) {
                    float inv = 1.0f / sum;
                    s->skin.boneWeights[v] *= inv;
                    s->skin.boneWeights[v+1] *= inv;
                    s->skin.boneWeights[v+2] *= inv;
                    s->skin.boneWeights[v+3] *= inv;
                }
            }
            return self;
        });

        proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::throwTypeError("SkinData.clone: not an instance");
            return wrapSkinData(s->skin);
        });

        proto.def("validate", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::throwTypeError("SkinData.validate: not an instance");
            if (!a.empty()) {
                auto* m = unwrapMesh(a[0]);
                if (m && !m->mesh.empty()) {
                    auto v = bromesh::validateSkin(m->mesh, s->skin);
                    return makeSkinValidationObject(v);
                }
            }
            auto v = bromesh::validateSkin(s->skin);
            return makeSkinValidationObject(v);
        });
    });

    skinCls.setStatic("validate", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        bromesh::MeshData mesh;
        bromesh::SkinData skin;
        if (a.size() > 0) {
            auto* m = unwrapMesh(a[0]);
            if (m) mesh = m->mesh;
        }
        if (a.size() > 1) {
            auto* s = unwrapSkinData(a[1]);
            if (s) skin = s->skin;
        }
        auto v = bromesh::validateSkin(mesh, skin);
        return makeSkinValidationObject(v);
    }, 2, "validate"));

    skinCls.setStatic("transfer", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("SkinData.transfer: targetMesh, sourceMesh, sourceSkin required");
        auto* tgt = unwrapMesh(a[0]);
        auto* src = unwrapMesh(a[1]);
        auto* skin = unwrapSkinData(a[2]);
        if (!tgt || !src || !skin) return ev::throwTypeError("SkinData.transfer: invalid arguments");
        float maxDist = a.size() > 3 ? static_cast<float>(numAt(a, 3)) : 0.0f;
        auto res = bromesh::transferSkinWeights(tgt->mesh, src->mesh, skin->skin, maxDist);
        return wrapSkinData(res);
    }, 4, "transfer"));

    // =========================================================================
    // Joint Class
    // =========================================================================
    jointCls.install("Joint", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostJoint>();
        if (!a.empty() && ev::isObject(a[0])) {
            readBoneObject(a[0], h->bone);
            Value idxVal = ev::getProperty(a[0], "index");  // a[0]: rooted slot
            if (ev::isNumber(idxVal)) h->index = static_cast<int>(ev::toDouble(idxVal));
        }
        return g_jointClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("name",
            [](Value self, std::span<const Value>) -> Value {
                auto* j = unwrapJoint(self);
                return ev::fromUtf8(j ? j->bone.name : "");
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* j = unwrapJoint(self);
                if (j && !a.empty()) j->bone.name = ev::toUtf8(a[0]);
                return ev::undefined();
            });

        proto.accessor("parent",
            [](Value self, std::span<const Value>) -> Value {
                auto* j = unwrapJoint(self);
                return ev::fromDouble(j ? static_cast<double>(j->bone.parent) : -1.0);
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* j = unwrapJoint(self);
                if (j && !a.empty()) j->bone.parent = i32At(a, 0);
                return ev::undefined();
            });

        proto.accessor("index",
            [](Value self, std::span<const Value>) -> Value {
                auto* j = unwrapJoint(self);
                return ev::fromDouble(j ? static_cast<double>(j->index) : -1.0);
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* j = unwrapJoint(self);
                if (j && !a.empty()) j->index = i32At(a, 0);
                return ev::undefined();
            });

        proto.accessor("localT", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.localT, 3);
        });

        proto.accessor("localR", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.localR, 4);
        });

        proto.accessor("localS", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.localS, 3);
        });

        proto.accessor("inverseBind", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.inverseBind, 16);
        });
    });

    // =========================================================================
    // Skeleton Class
    // =========================================================================
    skelCls.install("Skeleton", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostSkeleton>();
        if (!a.empty() && ev::isObject(a[0])) {
            // a[0] is the rooted slot; the lists and their elements are
            // rooted here, since every read allocates.
            Rooted bones(ev::getProperty(a[0], "bones"));
            if (ev::isObject(bones)) {
                Value lenVal = ev::getProperty(bones, "length");
                size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
                h->skeleton.bones.reserve(n);
                for (size_t i = 0; i < n; ++i) {
                    Value bVal = ev::getElement(bones, static_cast<uint32_t>(i));
                    bromesh::Bone b;
                    if (auto* j = unwrapJoint(bVal)) {
                        b = j->bone;
                    } else if (ev::isObject(bVal)) {
                        readBoneObject(bVal, b);
                    }
                    h->skeleton.bones.push_back(b);
                }
            }

            Rooted sockets(ev::getProperty(a[0], "sockets"));
            if (ev::isObject(sockets)) {
                Value lenVal = ev::getProperty(sockets, "length");
                size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
                for (size_t i = 0; i < n; ++i) {
                    Rooted s(ev::getElement(sockets, static_cast<uint32_t>(i)));
                    if (!ev::isObject(s)) continue;
                    bromesh::Socket sock;
                    if (Value v = ev::getProperty(s, "name"); ev::isString(v)) sock.name = ev::toUtf8(v);
                    Value boneVal = ev::getProperty(s, "boneIndex");
                    if (ev::isUndefined(boneVal)) boneVal = ev::getProperty(s, "bone");
                    if (ev::isNumber(boneVal)) sock.bone = static_cast<int>(ev::toDouble(boneVal));
                    if (Value v = ev::getProperty(s, "offset"); !ev::isUndefined(v)) {
                        auto off = toFloatVector(v);
                        if (off.size() >= 16) {
                            for (int k = 0; k < 16; ++k) sock.offset[k] = off[k];
                        }
                    }
                    h->skeleton.sockets.push_back(sock);
                }
            }
        }
        return g_skeletonClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            return ev::fromDouble(s ? static_cast<double>(s->skeleton.bones.size()) : 0.0);
        });

        proto.accessor("socketCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            return ev::fromDouble(s ? static_cast<double>(s->skeleton.sockets.size()) : 0.0);
        });

        proto.accessor("bones", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(s->skeleton.bones.size(), [&](size_t i) {
                return wrapJoint(s->skeleton.bones[i], static_cast<int>(i));
            });
        });

        proto.accessor("sockets", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(s->skeleton.sockets.size(), [&](size_t i) {
                const auto& sock = s->skeleton.sockets[i];
                ObjectBuilder obj;
                obj.set("name", sock.name);
                obj.set("bone", static_cast<double>(sock.bone));
                obj.set("boneIndex", static_cast<double>(sock.bone));
                ev::Persistent off(makeFloat32Array(sock.offset, 16));
                obj.set("offset", off.get());
                return obj.build();
            });
        });

        proto.def("findBone", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            std::string name = ev::toUtf8(a[0]);
            return ev::fromDouble(static_cast<double>(s->skeleton.findBone(name)));
        });

        proto.def("findBoneBySuffix", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            std::string suffix = ev::toUtf8(a[0]);
            for (size_t i = 0; i < s->skeleton.bones.size(); ++i) {
                const std::string& bname = s->skeleton.bones[i].name;
                if (bname.size() >= suffix.size() &&
                    bname.compare(bname.size() - suffix.size(), suffix.size(), suffix) == 0) {
                    return ev::fromDouble(static_cast<double>(i));
                }
            }
            return ev::fromDouble(-1.0);
        });

        proto.def("findSocket", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            std::string name = ev::toUtf8(a[0]);
            for (size_t i = 0; i < s->skeleton.sockets.size(); ++i) {
                if (s->skeleton.sockets[i].name == name) return ev::fromDouble(static_cast<double>(i));
            }
            return ev::fromDouble(-1.0);
        });

        proto.def("addSocket", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            bromesh::Socket sock;
            if (ev::isObject(a[0]) && a.size() == 1) {
                // Each read is consumed before the next allocating one.
                sock.name = ev::toUtf8(ev::getProperty(a[0], "name"));
                Value bVal = ev::getProperty(a[0], "bone");
                if (ev::isUndefined(bVal)) bVal = ev::getProperty(a[0], "boneIndex");
                sock.bone = ev::isNumber(bVal) ? static_cast<int>(ev::toDouble(bVal)) : 0;
                auto off = toFloatVector(ev::getProperty(a[0], "offset"));
                if (off.size() >= 16) for (int i = 0; i < 16; ++i) sock.offset[i] = off[i];
            } else {
                sock.name = ev::toUtf8(a[0]);
                sock.bone = a.size() > 1 ? i32At(a, 1) : 0;
                if (a.size() > 2) {
                    auto off = toFloatVector(a[2]);
                    if (off.size() >= 16) for (int i = 0; i < 16; ++i) sock.offset[i] = off[i];
                }
            }
            s->skeleton.sockets.push_back(sock);
            return ev::fromDouble(static_cast<double>(s->skeleton.sockets.size() - 1));
        });

        // addRigifySockets() -> number — append the standard attachment
        // sockets (hands, feet, head, spine) for a Rigify / Mixamo-named
        // skeleton, matching "ORG-"/"DEF-"/bare/"mixamorig:" bone spellings.
        // Returns how many were added. Dropped by the bronze port
        // (bro docs/transition-drift.md H7).
        proto.def("addRigifySockets", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return ev::throwTypeError("Skeleton.addRigifySockets: not an instance");
            return ev::fromDouble(static_cast<double>(bromesh::addRigifySockets(s->skeleton)));
        });

        proto.def("bindPose", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return ev::throwTypeError("Skeleton.bindPose: not an instance");
            return wrapPose(bromesh::bindPose(s->skeleton));
        });

        proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return ev::throwTypeError("Skeleton.clone: not an instance");
            return wrapSkeleton(s->skeleton);
        });
    });

    skelCls.setStatic("fromBones", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapSkeleton(bromesh::Skeleton{});
        ObjectBuilder opts;
        opts.set("bones", a[0]);
        const Value args[1] = {opts.build()};
        ev::CallResult r = ev::call(g_skeletonClass.constructor(), ev::undefined(), std::span<const Value>(args, 1));
        return r.thrown ? ev::throwValue(r.value) : r.value;
    }, 1, "fromBones"));

    // =========================================================================
    // SkeletonRig (Rig / RigSpec) Class
    // =========================================================================
    rigCls.install("SkeletonRig", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostSkeletonRig>();
        if (!a.empty()) {
            if (ev::isString(a[0])) {
                std::string type = ev::toUtf8(a[0]);
                if (type == "humanoid") h->spec = bromesh::builtinHumanoidSpec();
                else if (type == "quadruped") h->spec = bromesh::builtinQuadrupedSpec();
                else h->spec.name = type;
            } else if (ev::isObject(a[0])) {
                Value tVal = ev::getProperty(a[0], "type");
                if (ev::isUndefined(tVal)) tVal = ev::getProperty(a[0], "name");
                if (ev::isString(tVal)) {
                    std::string t = ev::toUtf8(tVal);
                    if (t == "humanoid") h->spec = bromesh::builtinHumanoidSpec();
                    else if (t == "quadruped") h->spec = bromesh::builtinQuadrupedSpec();
                    else h->spec.name = t;
                }
            }
        } else {
            h->spec = bromesh::builtinHumanoidSpec();
        }
        return g_skeletonRigClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("type", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromUtf8(r ? r->spec.name : "");
        });
        proto.accessor("name", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromUtf8(r ? r->spec.name : "");
        });
        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromDouble(r ? static_cast<double>(r->spec.bones.size()) : 0.0);
        });
        proto.accessor("landmarkCount", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromDouble(r ? static_cast<double>(r->spec.landmarks.size()) : 0.0);
        });
        proto.accessor("socketCount", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromDouble(r ? static_cast<double>(r->spec.sockets.size()) : 0.0);
        });
        proto.accessor("symmetric", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromBool(r && r->spec.symmetric);
        });
        proto.def("toJSON", 0, [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromUtf8(r ? bromesh::serializeRigSpecJSON(r->spec) : "");
        });
        proto.def("landmarkNames", 0, [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            if (!r) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(r->spec.landmarks.size(), [&](size_t i) {
                return ev::fromUtf8(r->spec.landmarks[i].name);
            });
        });
        proto.def("boneNames", 0, [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            if (!r) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(r->spec.bones.size(), [&](size_t i) {
                return ev::fromUtf8(r->spec.bones[i].name);
            });
        });
    });

    // Aliases on globalThis
    rigCls.alias("RigSpec");
    rigCls.alias("Rig");

    // Static Rig methods
    auto bindRigStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        rigCls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindRigStatic("spec", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapSkeletonRig(bromesh::RigSpec{});
        std::string name = ev::toUtf8(a[0]);
        return wrapSkeletonRig(bromesh::builtinRigSpec(name));
    });

    bindRigStatic("specFromJSON", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapSkeletonRig(bromesh::RigSpec{});
        std::string json = ev::toUtf8(a[0]);
        return wrapSkeletonRig(bromesh::parseRigSpecJSON(json));
    });

    bindRigStatic("specFromFile", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.specFromFile: path required");
        std::string path = ev::toUtf8(a[0]);
        return wrapSkeletonRig(bromesh::loadRigSpecFile(path));
    });

    bindRigStatic("detectHumanoid", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.detectHumanoid: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.detectHumanoid: mesh required");
        auto lm = bromesh::detectHumanoidLandmarks(m->mesh);
        return landmarksToObject(lm);
    });

    bindRigStatic("detectLandmarks", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.detectLandmarks: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.detectLandmarks: mesh required");
        auto lm = bromesh::detectHumanoidLandmarks(m->mesh);
        return landmarksToObject(lm);
    });

    bindRigStatic("detectQuadruped", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.detectQuadruped: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.detectQuadruped: mesh required");
        auto lm = bromesh::detectQuadrupedLandmarks(m->mesh);
        return landmarksToObject(lm);
    });

    bindRigStatic("missingLandmarks", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return hostArrayOf(std::span<const Value>{});
        auto* r = unwrapSkeletonRig(a[0]);
        bromesh::RigSpec spec = r ? r->spec : bromesh::builtinHumanoidSpec();
        bromesh::Landmarks lm = landmarksFromObject(a[1]);
        auto missing = bromesh::missingLandmarks(spec, lm);
        return hostArrayOf(missing.size(), [&](size_t i) {
            return ev::fromUtf8(missing[i]);
        });
    });

    bindRigStatic("fitSkeleton", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("Rig.fitSkeleton: spec, landmarks, and mesh required");
        auto* r = unwrapSkeletonRig(a[0]);
        bromesh::RigSpec spec = r ? r->spec : bromesh::builtinHumanoidSpec();
        bromesh::Landmarks lm = landmarksFromObject(a[1]);
        auto* m = unwrapMesh(a[2]);
        if (!m) return ev::throwTypeError("Rig.fitSkeleton: mesh required");
        return wrapSkeleton(bromesh::fitSkeleton(spec, lm, m->mesh));
    });

    bindRigStatic("autoRig", 4, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.autoRig: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.autoRig: mesh must be a Mesh instance");

        bromesh::RigSpec spec = bromesh::builtinHumanoidSpec();
        bromesh::Landmarks lm;
        bool hasLandmarks = false;
        Rooted opts(ev::undefined());  // read field by field below

        if (a.size() >= 3 && unwrapSkeletonRig(a[1])) {
            spec = unwrapSkeletonRig(a[1])->spec;
            if (ev::isObject(a[2])) {
                lm = landmarksFromObject(a[2]);
                hasLandmarks = true;
            }
            if (a.size() > 3 && ev::isObject(a[3])) opts.p.set(a[3]);
        } else if (a.size() > 1 && ev::isObject(a[1])) {
            opts.p.set(a[1]);
            Value specVal = ev::getProperty(opts, "spec");
            if (auto* r = unwrapSkeletonRig(specVal)) spec = r->spec;
            else {
                Value typeVal = ev::getProperty(opts, "rigType");
                if (ev::isString(typeVal)) {
                    spec = bromesh::builtinRigSpec(ev::toUtf8(typeVal));
                    if (spec.name.empty()) spec.name = ev::toUtf8(typeVal);
                }
            }
            Value lmVal = ev::getProperty(opts, "landmarks");
            if (ev::isObject(lmVal)) {
                lm = landmarksFromObject(lmVal);
                hasLandmarks = true;
            }
        }

        if (!hasLandmarks) {
            if (spec.name == "quadruped") lm = bromesh::detectQuadrupedLandmarks(m->mesh);
            else lm = bromesh::detectHumanoidLandmarks(m->mesh);
        }

        bromesh::WeightingOptions wopts = weightingOptionsFromObject(opts);

        auto res = bromesh::autoRig(m->mesh, spec, lm, wopts);
        ObjectBuilder obj;
        {
            ev::Persistent skel(wrapSkeleton(std::move(res.skeleton)));
            obj.set("skeleton", skel.get());
            ev::Persistent skin(wrapSkinData(std::move(res.skin)));
            obj.set("skin", skin.get());
        }
        obj.set("methodUsed", res.methodUsed == bromesh::WeightingMethod::BBW ? "bbw" :
                              res.methodUsed == bromesh::WeightingMethod::BoneHeat ? "boneHeat" : "voxelBind");
        obj.set("missingLandmarks", hostArrayOf(res.missingLandmarks.size(), [&](size_t i) {
            return ev::fromUtf8(res.missingLandmarks[i]);
        }));
        obj.set("warnings", hostArrayOf(res.warnings.size(), [&](size_t i) {
            return ev::fromUtf8(res.warnings[i]);
        }));
        return obj.build();
    });

    bindRigStatic("generateLocomotionCycle", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Rig.generateLocomotionCycle: skeleton and spec required");
        auto* s = unwrapSkeleton(a[0]);
        auto* r = unwrapSkeletonRig(a[1]);
        if (!s || !r) return ev::throwTypeError("Rig.generateLocomotionCycle: invalid arguments");
        bromesh::LocomotionParams params = locomotionParamsFromValue(a.size() > 2 ? a[2] : ev::undefined());
        auto anim = bromesh::generateLocomotionCycle(s->skeleton, r->spec, params);
        return wrapAnimation(std::move(anim));
    });

    bindRigStatic("transferWeights", 3, [](Value, std::span<const Value> a) -> Value {
        return callMethod(g_skinDataClass.constructor(), "transfer", a);
    });

    // =========================================================================
    // VoxelChunk Class
    // =========================================================================
    voxelCls.install("VoxelChunk", 4, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        int dx = r.getInt(0, 16);
        int dy = r.getInt(1, 16);
        int dz = r.getInt(2, 16);
        float cs = static_cast<float>(r.getDouble(3, 1.0));
        auto h = std::make_unique<HostVoxelChunk>();
        h->chunk = std::make_unique<bromesh::VoxelChunk>(dx > 0 ? dx : 16, dy > 0 ? dy : 16, dz > 0 ? dz : 16, cs > 0.0f ? cs : 1.0f);
        return g_voxelChunkClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("sizeX", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->sizeX()) : 0.0);
        });
        proto.accessor("sizeY", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->sizeY()) : 0.0);
        });
        proto.accessor("sizeZ", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->sizeZ()) : 0.0);
        });
        proto.accessor("cellSize", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->cellSize()) : 1.0);
        });
        proto.accessor("isDirty",
            [](Value self, std::span<const Value>) -> Value {
                auto* v = unwrapVoxelChunk(self);
                return ev::fromBool(v && v->chunk ? v->chunk->isDirty() : false);
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* v = unwrapVoxelChunk(self);
                if (v && v->chunk) {
                    if (!a.empty() && ev::toBool(a[0])) v->chunk->markDirty();
                    else v->chunk->clearDirty();
                }
                return ev::undefined();
            });

        proto.def("clearDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->clearDirty();
            return self;
        });
        proto.def("markDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->markDirty();
            return self;
        });

        proto.def("set", 4, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.set: not an instance");
            ArgReader r(a);
            v->chunk->setVoxel(r.getInt(0, 0), r.getInt(1, 0), r.getInt(2, 0), static_cast<uint8_t>(r.getInt(3, 0)));
            return self;
        });
        proto.def("setVoxel", 4, [](Value self, std::span<const Value> a) -> Value {
            return callMethod(self, "set", a);
        });

        proto.def("get", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.get: not an instance");
            ArgReader r(a);
            return ev::fromDouble(static_cast<double>(v->chunk->getVoxel(r.getInt(0, 0), r.getInt(1, 0), r.getInt(2, 0))));
        });
        proto.def("getVoxel", 3, [](Value self, std::span<const Value> a) -> Value {
            return callMethod(self, "get", a);
        });

        proto.def("fill", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.fill: not an instance");
            v->chunk->fill(static_cast<uint8_t>(i32At(a, 0)));
            return self;
        });

        proto.def("markDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->markDirty();
            return self;
        });
        proto.def("clearDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->clearDirty();
            return self;
        });

        proto.def("data", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.data: not an instance");
            size_t sz = static_cast<size_t>(v->chunk->sizeX()) * v->chunk->sizeY() * v->chunk->sizeZ();
            return makeUint8Array(v->chunk->data(), sz);
        });

        proto.def("setData", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk || a.empty()) return ev::throwTypeError("VoxelChunk.setData: not an instance or empty");
            Rooted selfP(self);  // a plain-array argument's reads allocate
            std::vector<uint8_t> d = toUint8Vector(a[0]);
            size_t sz = static_cast<size_t>(v->chunk->sizeX()) * v->chunk->sizeY() * v->chunk->sizeZ();
            std::memcpy(v->chunk->data(), d.data(), std::min(sz, d.size()));
            v->chunk->markDirty();
            return selfP.get();
        });

        proto.def("buildMesh", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.buildMesh: not an instance");
            const float* pal = nullptr;
            int count = 0;
            std::vector<float> palVec;
            if (!a.empty()) {
                palVec = toFloatVector(a[0]);
                if (!palVec.empty()) {
                    pal = palVec.data();
                    count = a.size() > 1 ? i32At(a, 1) : static_cast<int>(palVec.size() / 4);
                }
            }
            return wrapMesh(v->chunk->buildMesh(pal, count));
        });
        proto.def("toMesh", 0, [](Value self, std::span<const Value> a) -> Value {
            return callMethod(self, "buildMesh", a);
        });
    });
}

} // namespace bromesh::api
