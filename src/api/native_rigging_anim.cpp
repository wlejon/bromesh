#include "host_mesh_internal.h"

namespace bromesh::api {

HostClass g_poseClass;
HostClass g_animationClass;

namespace {

Value animChannelsToObject(const bromesh::Animation& anim) {
    return hostArrayOf(anim.channels.size(), [&](size_t i) {
        const auto& ch = anim.channels[i];
        ObjectBuilder obj;
        obj.set("boneIndex", static_cast<double>(ch.boneIndex));
        obj.set("path", ch.path == bromesh::AnimChannel::Path::Translation ? "translation" :
                        ch.path == bromesh::AnimChannel::Path::Rotation ? "rotation" : "scale");
        obj.set("interp", ch.interp == bromesh::AnimChannel::Interp::Step ? "step" :
                          ch.interp == bromesh::AnimChannel::Interp::CubicSpline ? "cubicspline" : "linear");
        ev::Persistent t(makeFloat32Array(ch.times.data(), ch.times.size()));
        obj.set("times", t.get());
        ev::Persistent v(makeFloat32Array(ch.values.data(), ch.values.size()));
        obj.set("values", v.get());
        return obj.build();
    });
}

// Most bones a Pose(boneCount) allocates: 10 floats each.
constexpr double kMaxBones = 1048576.0;

// The skeleton/pose pair the world-matrix and IK walks index unchecked:
// every parent must name a bone of the skeleton (or be negative, a root), and
// the pose must carry at least one TRS record per skeleton bone.
bool rigFits(const bromesh::Skeleton& s, const bromesh::Pose& p, const char* fn) {
    const size_t bones = s.bones.size();
    if (p.data.size() < bones * 10) {
        ev::throwRangeError(std::string(fn) + ": the pose has " + std::to_string(p.data.size() / 10) +
                            " bones, the skeleton " + std::to_string(bones));
        return false;
    }
    for (size_t i = 0; i < bones; ++i) {
        const int parent = s.bones[i].parent;
        if (parent >= 0 && static_cast<size_t>(parent) >= bones) {
            ev::throwRangeError(std::string(fn) + ": bone " + std::to_string(i) + " has parent " +
                                std::to_string(parent) + ", out of range for " + std::to_string(bones) +
                                " bones");
            return false;
        }
    }
    return true;
}

// A bone mask is read one byte per bone of the pose.
bool maskFits(const std::vector<uint8_t>& mask, size_t bones, const char* fn) {
    if (mask.empty() || mask.size() >= bones) return true;
    ev::throwRangeError(std::string(fn) + ": the mask has " + std::to_string(mask.size()) +
                        " entries, the pose " + std::to_string(bones) + " bones");
    return false;
}

} // namespace

void initRiggingAnim(HostClass& poseCls, HostClass& animCls, HostClass& meshCls,
                     ObjectBuilder& ikBuilder) {
    // =========================================================================
    // Pose Class
    // =========================================================================
    poseCls.install("Pose", 2, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostPose>();
        if (!a.empty()) {
            if (auto* skel = unwrapSkeleton(a[0])) {
                h->pose = bromesh::bindPose(skel->skeleton);
            } else if (ev::isNumber(a[0])) {
                size_t bc = 0;
                if (!countArg(a, 0, "Pose: boneCount", 0, kMaxBones, bc)) return ev::undefined();
                h->pose.data.assign(bc * 10, 0.0f);
                for (size_t i = 0; i < bc; ++i) {
                    h->pose.data[i * 10 + 6] = 1.0f; // rw = 1
                    h->pose.data[i * 10 + 7] = 1.0f; // sx = 1
                    h->pose.data[i * 10 + 8] = 1.0f; // sy = 1
                    h->pose.data[i * 10 + 9] = 1.0f; // sz = 1
                }
            } else {
                std::vector<float> d = toFloatVector(a[0]);
                size_t bc = 0;
                if (!countArg(a, 1, "Pose: boneCount", 0, kMaxBones, bc)) return ev::undefined();
                if (d.size() < bc * 10) d.resize(bc * 10, 0.0f);
                h->pose.data = std::move(d);
            }
        }
        return g_poseClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("data",
            [](Value self, std::span<const Value>) -> Value {
                auto* p = unwrapPose(self);
                if (!p) return ev::undefined();
                return makeFloat32Array(p->pose.data.data(), p->pose.data.size());
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* p = unwrapPose(self);
                if (!p || a.empty()) return ev::undefined();
                p->pose.data = toFloatVector(a[0]);
                return ev::undefined();
            });

        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPose(self);
            return ev::fromDouble(p ? static_cast<double>(p->pose.boneCount()) : 0.0);
        });

        proto.def("computeWorldMatrices", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPose(self);
            if (!p) return ev::throwTypeError("Pose.computeWorldMatrices: not a Pose instance");
            if (a.empty()) return ev::throwTypeError("Pose.computeWorldMatrices: Skeleton required");
            auto* s = unwrapSkeleton(a[0]);
            if (!s) return ev::throwTypeError("Pose.computeWorldMatrices: argument must be a Skeleton");
            if (!rigFits(s->skeleton, p->pose, "Pose.computeWorldMatrices")) return ev::undefined();
            std::vector<float> outWorld;
            bromesh::computeWorldMatrices(s->skeleton, p->pose, outWorld);
            return makeFloat32Array(outWorld.data(), outWorld.size());
        });

        proto.def("computeSkinningMatrices", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPose(self);
            if (!p) return ev::throwTypeError("Pose.computeSkinningMatrices: not a Pose instance");
            if (a.empty()) return ev::throwTypeError("Pose.computeSkinningMatrices: Skeleton required");
            auto* s = unwrapSkeleton(a[0]);
            if (!s) return ev::throwTypeError("Pose.computeSkinningMatrices: argument must be a Skeleton");
            if (!rigFits(s->skeleton, p->pose, "Pose.computeSkinningMatrices")) return ev::undefined();
            std::vector<float> outSkinning;
            bromesh::computeSkinningMatrices(s->skeleton, p->pose, outSkinning);
            return makeFloat32Array(outSkinning.data(), outSkinning.size());
        });

        proto.def("socketWorld", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPose(self);
            if (!p || a.size() < 2) return ev::null();
            auto* s = unwrapSkeleton(a[0]);
            if (!s) return ev::null();
            std::string name = ev::toUtf8(a[1]);
            if (!rigFits(s->skeleton, p->pose, "Pose.socketWorld")) return ev::undefined();
            auto res = bromesh::socketWorldMatrix(s->skeleton, p->pose, name);
            if (!res.has_value()) return ev::null();
            return makeFloat32Array(res->data(), 16);
        });

        proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPose(self);
            if (!p) return ev::throwTypeError("Pose.clone: not a Pose instance");
            return wrapPose(p->pose);
        });
    });

    // Static blend on Pose
    poseCls.setStatic("blend", hostFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("Pose.blend: a, b, weight required");
        auto* pa = unwrapPose(a[0]);
        auto* pb = unwrapPose(a[1]);
        if (!pa || !pb) return ev::throwTypeError("Pose.blend: arguments must be Pose instances");
        float w = static_cast<float>(numAt(a, 2));
        const uint8_t* mask = nullptr;
        std::vector<uint8_t> maskVec;
        if (a.size() > 3 && !ev::isUndefined(a[3])) {
            maskVec = toUint8Vector(a[3]);
            if (!maskFits(maskVec, pa->pose.boneCount(), "Pose.blend")) return ev::undefined();
            if (!maskVec.empty()) mask = maskVec.data();
        }
        bromesh::blendPoses(pa->pose, pb->pose, w, mask);
        return a[0];
    }, 4, "blend"));

    poseCls.setStatic("blendN", hostFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Pose.blendN: poses and weights required");
        std::vector<const bromesh::Pose*> posePtrs;
        if (ev::isObject(a[0])) {
            size_t n = 0;
            if (!listLength(ev::getProperty(a[0], "length"), "Pose.blendN: poses", n)) return ev::undefined();
            for (size_t i = 0; i < n; ++i) {
                Value elem = ev::getElement(a[0], static_cast<uint32_t>(i));
                auto* p = unwrapPose(elem);
                if (p) posePtrs.push_back(&p->pose);
            }
        }
        if (posePtrs.empty()) return ev::throwTypeError("Pose.blendN: no valid poses provided");

        std::vector<float> weights = toFloatVector(a[1]);
        if (weights.size() != posePtrs.size()) {
            return ev::throwTypeError("Pose.blendN: weights count must match poses count");
        }

        const uint8_t* mask = nullptr;
        std::vector<uint8_t> maskVec;
        if (a.size() > 2 && !ev::isUndefined(a[2])) {
            maskVec = toUint8Vector(a[2]);
            if (!maskFits(maskVec, posePtrs[0]->boneCount(), "Pose.blendN")) return ev::undefined();
            if (!maskVec.empty()) mask = maskVec.data();
        }

        bromesh::Pose out;
        bromesh::blendPosesN(posePtrs.data(), weights.data(), posePtrs.size(), out, mask);
        return wrapPose(std::move(out));
    }, 3, "blendN"));

    // =========================================================================
    // AnimationClip Class
    // =========================================================================
    animCls.install("AnimationClip", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostAnimation>();
        if (!a.empty() && ev::isObject(a[0])) {
            // a[0] is the rooted slot; the channel list and each channel are
            // rooted here, since every read allocates.
            Value nameVal = ev::getProperty(a[0], "name");
            if (ev::isString(nameVal)) h->animation.name = ev::toUtf8(nameVal);
            Value durVal = ev::getProperty(a[0], "duration");
            if (ev::isNumber(durVal)) h->animation.duration = static_cast<float>(ev::toDouble(durVal));

            Rooted chsVal(ev::getProperty(a[0], "channels"));
            if (ev::isObject(chsVal)) {
                size_t n = 0;
                if (!listLength(ev::getProperty(chsVal, "length"), "AnimationClip: channels", n)) {
                    return ev::undefined();
                }
                for (size_t i = 0; i < n; ++i) {
                    Rooted chVal(ev::getElement(chsVal, static_cast<uint32_t>(i)));
                    if (ev::isObject(chVal)) {
                        bromesh::AnimChannel ch;
                        Value biVal = ev::getProperty(chVal, "boneIndex");
                        if (ev::isNumber(biVal)) ch.boneIndex = satInt(ev::toDouble(biVal));
                        Value pVal = ev::getProperty(chVal, "path");
                        if (ev::isString(pVal)) {
                            std::string ps = ev::toUtf8(pVal);
                            if (ps == "rotation") ch.path = bromesh::AnimChannel::Path::Rotation;
                            else if (ps == "scale") ch.path = bromesh::AnimChannel::Path::Scale;
                            else ch.path = bromesh::AnimChannel::Path::Translation;
                        }
                        Value inVal = ev::getProperty(chVal, "interp");
                        if (ev::isString(inVal)) {
                            std::string is = ev::toUtf8(inVal);
                            if (is == "step" || is == "STEP") ch.interp = bromesh::AnimChannel::Interp::Step;
                            else if (is == "cubic" || is == "cubicspline" || is == "CUBICSPLINE") ch.interp = bromesh::AnimChannel::Interp::CubicSpline;
                            else ch.interp = bromesh::AnimChannel::Interp::Linear;
                        }
                        if (Value v = ev::getProperty(chVal, "times"); !ev::isUndefined(v)) ch.times = toFloatVector(v);
                        if (Value v = ev::getProperty(chVal, "values"); !ev::isUndefined(v)) ch.values = toFloatVector(v);
                        // Sampling reads `stride` values per key (three
                        // stride-runs per key for cubicspline) unchecked.
                        const size_t stride = ch.path == bromesh::AnimChannel::Path::Rotation ? 4 : 3;
                        const size_t perKey =
                            stride * (ch.interp == bromesh::AnimChannel::Interp::CubicSpline ? 3 : 1);
                        if (ch.values.size() < ch.times.size() * perKey) {
                            return ev::throwRangeError(
                                "AnimationClip: channel " + std::to_string(i) + " has " +
                                std::to_string(ch.times.size()) + " keys but " + std::to_string(ch.values.size()) +
                                " values; it needs " + std::to_string(perKey) + " per key");
                        }
                        h->animation.channels.push_back(std::move(ch));
                    }
                }
            }
        }
        return g_animationClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("name",
            [](Value self, std::span<const Value>) -> Value {
                auto* a = unwrapAnimation(self);
                return ev::fromUtf8(a ? a->animation.name : "");
            },
            [](Value self, std::span<const Value> v) -> Value {
                auto* a = unwrapAnimation(self);
                if (a && !v.empty()) a->animation.name = ev::toUtf8(v[0]);
                return ev::undefined();
            });

        proto.accessor("duration",
            [](Value self, std::span<const Value>) -> Value {
                auto* a = unwrapAnimation(self);
                return ev::fromDouble(a ? static_cast<double>(a->animation.duration) : 0.0);
            },
            [](Value self, std::span<const Value> v) -> Value {
                auto* a = unwrapAnimation(self);
                if (a && !v.empty()) a->animation.duration = static_cast<float>(numAt(v, 0));
                return ev::undefined();
            });

        proto.accessor("channelCount", [](Value self, std::span<const Value>) -> Value {
            auto* a = unwrapAnimation(self);
            return ev::fromDouble(a ? static_cast<double>(a->animation.channels.size()) : 0.0);
        });

        proto.accessor("channels", [](Value self, std::span<const Value>) -> Value {
            auto* a = unwrapAnimation(self);
            if (!a) return hostArrayOf(std::span<const Value>{});
            return animChannelsToObject(a->animation);
        });

        proto.def("evaluate", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* anim = unwrapAnimation(self);
            if (!anim) return ev::throwTypeError("AnimationClip.evaluate: not an instance");
            if (a.empty()) return ev::throwTypeError("AnimationClip.evaluate: skeleton required");
            auto* skel = unwrapSkeleton(a[0]);
            if (!skel) return ev::throwTypeError("AnimationClip.evaluate: argument must be a Skeleton");
            float t = a.size() > 1 ? static_cast<float>(numAt(a, 1)) : 0.0f;
            bool loop = true;
            if (a.size() > 2) {
                if (ev::isObject(a[2])) {
                    Value lv = ev::getProperty(a[2], "loop");
                    if (ev::isBool(lv)) loop = ev::toBool(lv);
                } else {
                    loop = boolAt(a, 2);
                }
            }
            bromesh::Pose p = bromesh::evaluateAnimation(skel->skeleton, anim->animation, t, loop);
            return wrapPose(std::move(p));
        });

        proto.def("evaluateInto", 4, [](Value self, std::span<const Value> a) -> Value {
            auto* anim = unwrapAnimation(self);
            if (!anim) return ev::throwTypeError("AnimationClip.evaluateInto: not an instance");
            if (a.size() < 3) return ev::throwTypeError("AnimationClip.evaluateInto: skeleton, time, pose required");
            auto* skel = unwrapSkeleton(a[0]);
            if (!skel) return ev::throwTypeError("AnimationClip.evaluateInto: skeleton required");
            float t = static_cast<float>(numAt(a, 1));
            bool loop = true;
            Value poseVal = ev::undefined();
            if (unwrapPose(a[2])) {
                poseVal = a[2];
                if (a.size() > 3) loop = boolAt(a, 3);
            } else {
                if (ev::isObject(a[2])) {
                    Value lv = ev::getProperty(a[2], "loop");
                    if (ev::isBool(lv)) loop = ev::toBool(lv);
                } else {
                    loop = boolAt(a, 2);
                }
                if (a.size() > 3) poseVal = a[3];
            }
            auto* p = unwrapPose(poseVal);
            if (!p) return ev::throwTypeError("AnimationClip.evaluateInto: pose required");
            bromesh::evaluateAnimationInto(skel->skeleton, anim->animation, t, loop, p->pose);
            return poseVal;
        });
    });

    // =========================================================================
    // IK Solvers (IK namespace)
    // =========================================================================
    ikBuilder.def("twoBone", 7, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 6) return ev::throwTypeError("IK.twoBone: skel, pose, root, mid, end, target required");
        auto* skel = unwrapSkeleton(a[0]);
        auto* pose = unwrapPose(a[1]);
        if (!skel || !pose) return ev::throwTypeError("IK.twoBone: skel and pose required");
        if (!rigFits(skel->skeleton, pose->pose, "IK.twoBone")) return ev::undefined();
        int root = satInt(ev::toDouble(a[2]));
        int mid = satInt(ev::toDouble(a[3]));
        int end = satInt(ev::toDouble(a[4]));
        auto tVec = toFloatVector(a[5]);
        if (tVec.size() < 3) return ev::throwTypeError("IK.twoBone: target must have [x, y, z]");
        float target[3] = {tVec[0], tVec[1], tVec[2]};
        float pole[3] = {0, 0, 0};
        const float* polePtr = nullptr;
        if (a.size() > 6 && !ev::isUndefined(a[6]) && !ev::isNull(a[6])) {
            auto pVec = toFloatVector(a[6]);
            if (pVec.size() >= 3) {
                pole[0] = pVec[0]; pole[1] = pVec[1]; pole[2] = pVec[2];
                polePtr = pole;
            }
        }
        bool ok = bromesh::solveTwoBoneIK(skel->skeleton, pose->pose, root, mid, end, target, polePtr);
        return ev::fromBool(ok);
    });

    ikBuilder.def("FABRIK", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 4) return ev::throwTypeError("IK.FABRIK: skel, pose, chain, target required");
        auto* skel = unwrapSkeleton(a[0]);
        auto* pose = unwrapPose(a[1]);
        if (!skel || !pose) return ev::throwTypeError("IK.FABRIK: skel and pose required");
        if (!rigFits(skel->skeleton, pose->pose, "IK.FABRIK")) return ev::undefined();
        std::vector<int> chain;
        if (ev::isObject(a[2])) {
            // A chain is a handful of bone indices; a longer `length` is a
            // bad argument, not a loop of billions of reads.
            size_t n = lengthValue(ev::getProperty(a[2], "length"));
            if (n > skel->skeleton.bones.size()) {
                return ev::throwRangeError("IK.FABRIK: the chain has " + std::to_string(n) +
                                           " entries, more than the skeleton's " +
                                           std::to_string(skel->skeleton.bones.size()) + " bones");
            }
            for (size_t i = 0; i < n; ++i) {
                chain.push_back(satInt(ev::toDouble(ev::getElement(a[2], static_cast<uint32_t>(i)))));
            }
        }
        auto tVec = toFloatVector(a[3]);
        if (tVec.size() < 3) return ev::throwTypeError("IK.FABRIK: target must have [x, y, z]");
        float target[3] = {tVec[0], tVec[1], tVec[2]};
        int iters = 10;
        float tol = 1e-3f;
        if (a.size() > 4 && ev::isObject(a[4])) {
            if (!countField(a[4], "iterations", "IK.FABRIK: iterations", 0, kMaxIterations, iters)) {
                return ev::undefined();
            }
            Value toV = ev::getProperty(a[4], "tolerance");
            if (ev::isNumber(toV)) tol = static_cast<float>(ev::toDouble(toV));
        }
        bool ok = bromesh::solveFABRIK(skel->skeleton, pose->pose, chain, target, iters, tol);
        return ev::fromBool(ok);
    });

    ikBuilder.def("lookAt", 5, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 4) return ev::throwTypeError("IK.lookAt: skel, pose, bone, target required");
        auto* skel = unwrapSkeleton(a[0]);
        auto* pose = unwrapPose(a[1]);
        if (!skel || !pose) return ev::throwTypeError("IK.lookAt: skel and pose required");
        if (!rigFits(skel->skeleton, pose->pose, "IK.lookAt")) return ev::undefined();
        int bone = satInt(ev::toDouble(a[2]));
        auto tVec = toFloatVector(a[3]);
        if (tVec.size() < 3) return ev::throwTypeError("IK.lookAt: target must have [x, y, z]");
        float target[3] = {tVec[0], tVec[1], tVec[2]};
        float fwd[3] = {0, 0, 1};
        float up[3] = {0, 1, 0};
        const float* fwdPtr = nullptr;
        const float* upPtr = nullptr;
        if (a.size() > 4 && ev::isObject(a[4])) {
            // Each read is consumed before the next one allocates.
            Value fV = ev::getProperty(a[4], "forward");
            if (ev::isObject(fV)) {
                auto fVec = toFloatVector(fV);
                if (fVec.size() >= 3) { fwd[0] = fVec[0]; fwd[1] = fVec[1]; fwd[2] = fVec[2]; fwdPtr = fwd; }
            }
            Value uV = ev::getProperty(a[4], "up");
            if (ev::isObject(uV)) {
                auto uVec = toFloatVector(uV);
                if (uVec.size() >= 3) { up[0] = uVec[0]; up[1] = uVec[1]; up[2] = uVec[2]; upPtr = up; }
            }
        }
        bool ok = bromesh::solveLookAt(skel->skeleton, pose->pose, bone, target, fwdPtr, upPtr);
        return ev::fromBool(ok);
    });

    // Object-form wrappers over the positional solvers. Every field read
    // allocates, so the receiver and each field are rooted as they are read,
    // and only turned back into plain Values for the call itself.
    auto forwardFields = [](Value self, Value opts, const char* method,
                            std::initializer_list<const char*> keys,
                            const std::function<void(ObjectBuilder&, Value)>& extra) -> Value {
        Rooted recv(self);
        Rooted o(opts);
        std::vector<Rooted> fields;
        fields.reserve(keys.size() + 1);
        for (const char* k : keys) fields.emplace_back(ev::getProperty(o, k));
        if (extra) {
            ObjectBuilder b;
            extra(b, o);
            fields.emplace_back(b.build());
        }
        std::vector<Value> args;
        args.reserve(fields.size());
        for (const auto& f : fields) args.push_back(f.get());
        return callMethod(recv, method, args);
    };

    ikBuilder.def("solveTwoBone", 1, [forwardFields](Value self, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isObject(a[0])) return ev::fromBool(false);
        return forwardFields(self, a[0], "twoBone",
                             {"skel", "pose", "root", "mid", "end", "targetPos", "poleVector"}, nullptr);
    });

    ikBuilder.def("solveFabrik", 1, [forwardFields](Value self, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isObject(a[0])) return ev::fromBool(false);
        return forwardFields(self, a[0], "FABRIK", {"skel", "pose", "chain", "targetPos"},
                             [](ObjectBuilder& b, Value o) {
                                 Rooted src(o);
                                 b.set("iterations", ev::getProperty(src, "maxIterations"));
                                 b.set("tolerance", ev::getProperty(src, "tolerance"));
                             });
    });

    ikBuilder.def("solveLookAt", 1, [forwardFields](Value self, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isObject(a[0])) return ev::fromBool(false);
        return forwardFields(self, a[0], "lookAt", {"skel", "pose", "bone", "targetPos"},
                             [](ObjectBuilder& b, Value o) {
                                 Rooted src(o);
                                 b.set("forward", ev::getProperty(src, "forward"));
                                 b.set("up", ev::getProperty(src, "up"));
                             });
    });

    // =========================================================================
    // Retargeting
    // =========================================================================
    animCls.setStatic("retarget", hostFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("AnimationClip.retarget: anim, srcSkel, dstSkel required");
        auto* anim = unwrapAnimation(a[0]);
        auto* src = unwrapSkeleton(a[1]);
        auto* dst = unwrapSkeleton(a[2]);
        if (!anim || !src || !dst) return ev::throwTypeError("AnimationClip.retarget: invalid arguments");
        auto retargeted = bromesh::retargetAnimation(anim->animation, src->skeleton, dst->skeleton);
        return wrapAnimation(std::move(retargeted));
    }, 3, "retarget"));

    animCls.alias("Animation");
    animCls.alias("SkeletalAnimation");

    // =========================================================================
    // Mesh Rigging Extensions
    // =========================================================================
    meshCls.prototype(); // ensures proto exists
    ObjectBuilder meshProto(g_meshClass.prototype());

    meshProto.def("applySkinning", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.applySkinning: not a Mesh instance");
        if (a.size() < 2) return ev::throwTypeError("Mesh.applySkinning: skinData and poseMatrices required");
        auto* s = unwrapSkinData(a[0]);
        if (!s) return ev::throwTypeError("Mesh.applySkinning: first argument must be SkinData");
        Rooted selfP(self);  // a plain-array argument's reads allocate
        std::vector<float> mats = toFloatVector(a[1]);
        if (mats.size() < s->skin.boneCount * 16) {
            return ev::throwTypeError("Mesh.applySkinning: not enough matrix floats for bones");
        }
        // The matrices are joint matrices (world x inverseBind, as
        // pose.computeSkinningMatrices returns); the inverse binds are in them.
        bromesh::applySkinning(m->mesh, s->skin, mats.data());
        return selfP.get();
    });

    meshProto.def("applyMorphTarget", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.applyMorphTarget: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.applyMorphTarget: target or name required");
        Rooted selfP(self);  // every read below may allocate
        if (ev::isObject(a[0])) {
            Rooted opt(a[0]);
            bromesh::MorphTarget mt;
            Value nVal = ev::getProperty(opt, "name");
            if (ev::isString(nVal)) mt.name = ev::toUtf8(nVal);
            Value dpVal = ev::getProperty(opt, "deltaPositions");
            if (!ev::isUndefined(dpVal)) mt.deltaPositions = toFloatVector(dpVal);
            Value dnVal = ev::getProperty(opt, "deltaNormals");
            if (!ev::isUndefined(dnVal)) mt.deltaNormals = toFloatVector(dnVal);
            float w = 1.0f;
            if (a.size() > 1 && ev::isNumber(a[1])) w = static_cast<float>(ev::toDouble(a[1]));
            else {
                Value wVal = ev::getProperty(opt, "weight");
                if (ev::isNumber(wVal)) w = static_cast<float>(ev::toDouble(wVal));
            }
            bromesh::applyMorphTarget(m->mesh, mt, w);
            return selfP.get();
        }
        if (a.size() < 3) return ev::throwTypeError("Mesh.applyMorphTarget: name, deltaPositions, weight required");
        bromesh::MorphTarget mt;
        mt.name = ev::toUtf8(a[0]);
        mt.deltaPositions = toFloatVector(a[1]);
        if (a.size() > 3) {
            mt.deltaNormals = toFloatVector(a[2]);
            float w = static_cast<float>(numAt(a, 3));
            bromesh::applyMorphTarget(m->mesh, mt, w);
        } else {
            float w = static_cast<float>(numAt(a, 2));
            bromesh::applyMorphTarget(m->mesh, mt, w);
        }
        return selfP.get();
    });

#if BROMESH_HAS_GLTF
    meshProto.def("saveGLTF", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.saveGLTF: not a Mesh instance");
        if (a.empty()) return ev::throwTypeError("Mesh.saveGLTF: path required");
        std::string path = resolveMeshWritePath(ev::toUtf8(a[0]));
        const bromesh::SkinData* skinPtr = nullptr;
        const bromesh::Skeleton* skelPtr = nullptr;
        std::vector<bromesh::Animation> anims;
        if (a.size() > 1 && ev::isObject(a[1])) {
            Rooted opts(a[1]);
            if (auto* s = unwrapSkinData(ev::getProperty(opts, "skin"))) skinPtr = &s->skin;
            if (auto* k = unwrapSkeleton(ev::getProperty(opts, "skeleton"))) skelPtr = &k->skeleton;
            Rooted anVal(ev::getProperty(opts, "animations"));
            if (ev::isObject(anVal)) {
                size_t n = 0;
                if (!listLength(ev::getProperty(anVal, "length"), "Mesh.saveGLTF: opts.animations", n)) {
                    return ev::undefined();
                }
                for (size_t i = 0; i < n; ++i) {
                    Value item = ev::getElement(anVal, static_cast<uint32_t>(i));
                    if (auto* an = unwrapAnimation(item)) anims.push_back(an->animation);
                }
            }
        }
        bool ok = bromesh::saveGLTF(m->mesh, skinPtr, skelPtr, anims, path);
        return ev::fromBool(ok);
    });

    meshCls.setStatic("loadGLTF", hostFunction([](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.loadGLTF: path required");
        std::string path = resolveMeshPath(ev::toUtf8(a[0]));
        auto scene = bromesh::loadGLTF(path);
        ObjectBuilder res;
        res.set("meshes", hostArrayOf(scene.meshes.size(), [&](size_t i) {
            return wrapMesh(std::move(scene.meshes[i]));
        }));
        res.set("skins", hostArrayOf(scene.skins.size(), [&](size_t i) {
            return wrapSkinData(std::move(scene.skins[i]));
        }));
        res.set("skeletons", hostArrayOf(scene.skeletons.size(), [&](size_t i) {
            return wrapSkeleton(std::move(scene.skeletons[i]));
        }));
        res.set("animations", hostArrayOf(scene.animations.size(), [&](size_t i) {
            return wrapAnimation(std::move(scene.animations[i]));
        }));
        res.set("meshSkeleton", hostArrayOf(scene.meshSkeleton.size(), [&](size_t i) {
            return ev::fromDouble(static_cast<double>(scene.meshSkeleton[i]));
        }));
        res.set("animationSkeleton", hostArrayOf(scene.animationSkeleton.size(), [&](size_t i) {
            return ev::fromDouble(static_cast<double>(scene.animationSkeleton[i]));
        }));
        return res.build();
    }, 1, "loadGLTF"));
#endif

}

void ensureRiggingClassesInstalled() {
    // Per thread, like ensureMeshClassesInstalled.
    static thread_local bool installed = false;
    if (installed) return;
    installed = true;

    ensureMeshClassesInstalled();

    ObjectBuilder ikObj;
    initRiggingCore(g_skinDataClass, g_skeletonClass, g_jointClass, g_skeletonRigClass, g_voxelChunkClass);
    initRiggingAnim(g_poseClass, g_animationClass, g_meshClass, ikObj);

    // Register IK on globalThis
    ev::registerGlobal("IK", ikObj.get());
    ev::GlobalValue gt = ev::globalValue("globalThis");
    if (gt.found && !gt.value.isUndefined() && ev::isObject(gt.value)) {
        ev::setProperty(gt.value, "IK", ikObj.get());
    }
}

} // namespace bromesh::api
