#pragma once

#include "embed/embed.h"
#include "host_class.h"
#include "object_builder.h"
#include "arg_reader.h"

#include <bromesh/mesh_data.h>
#include <bromesh/primitives/primitives.h>
#include <bromesh/primitives/par_primitives.h>
#include <bromesh/manipulation/transform.h>
#include <bromesh/manipulation/normals.h>
#include <bromesh/manipulation/weld.h>
#include <bromesh/manipulation/repair.h>
#include <bromesh/manipulation/simplify.h>
#include <bromesh/manipulation/subdivide.h>
#include <bromesh/manipulation/smooth.h>
#include <bromesh/manipulation/remesh.h>
#include <bromesh/manipulation/shrinkwrap.h>
#include <cstring>
#include <bromesh/manipulation/split_components.h>
#include <bromesh/manipulation/merge.h>
#include <bromesh/manipulation/skin.h>
#include <bromesh/manipulation/skin_transfer.h>
#include <bromesh/manipulation/sweep.h>
#include <bromesh/csg/boolean.h>
#include <bromesh/analysis/bbox.h>
#include <bromesh/analysis/sample.h>
#include <bromesh/analysis/raycast.h>
#include <bromesh/analysis/bvh.h>
#include <bromesh/analysis/intersect.h>
#include <bromesh/analysis/bake.h>
#include <bromesh/analysis/bake_texture.h>
#include <bromesh/analysis/bake_transfer.h>
#include <bromesh/analysis/convex_decomposition.h>
#include <bromesh/uv/projection.h>
#include <bromesh/uv/unwrap.h>
#include <bromesh/uv/uv_metrics.h>
#include <bromesh/optimization/meshlets.h>
#include <bromesh/optimization/analyze.h>
#include <bromesh/optimization/optimize.h>
#include <bromesh/optimization/spatial.h>
#include <bromesh/optimization/strips.h>
#include <bromesh/optimization/encode.h>
#include <bromesh/optimization/progressive.h>
#include <bromesh/manipulation/polygon.h>
#include <bromesh/reconstruction/reconstruct.h>
#include <bromesh/isosurface/marching_cubes.h>
#include <bromesh/isosurface/surface_nets.h>
#include <bromesh/isosurface/dual_contouring.h>
#include <bromesh/isosurface/transvoxel.h>
#include <bromesh/voxel/greedy_mesh.h>
#include <bromesh/voxel/voxel_chunk.h>
#include <bromesh/animation/pose.h>
#include <bromesh/animation/ik.h>
#include <bromesh/animation/retarget.h>
#include <bromesh/animation/locomotion.h>
#include <bromesh/rigging/rig_spec.h>
#include <bromesh/rigging/landmarks.h>
#include <bromesh/rigging/landmark_detect.h>
#include <bromesh/rigging/skeleton_fit.h>
#include <bromesh/rigging/auto_rig.h>
#include <bromesh/rigging/skin_validate.h>
#if BROMESH_HAS_GLTF
#include <bromesh/io/gltf.h>
#endif
#include <bromesh/isosurface/jit/sdf_node.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace bromesh::api {

inline constexpr uint32_t kHostMeshTag        = 0x4D455348u; // 'MESH'
inline constexpr uint32_t kHostMeshBVHTag     = 0x4D425648u; // 'MBVH'
inline constexpr uint32_t kHostProgressiveTag = 0x504D5348u; // 'PMSH'
inline constexpr uint32_t kHostSkinDataTag    = 0x534B494Eu; // 'SKIN'
inline constexpr uint32_t kHostSkeletonTag    = 0x534B454Cu; // 'SKEL'
inline constexpr uint32_t kHostJointTag       = 0x4A4F4E54u; // 'JONT'
inline constexpr uint32_t kHostSkeletonRigTag = 0x53524947u; // 'SRIG'
inline constexpr uint32_t kHostPoseTag        = 0x504F5345u; // 'POSE'
inline constexpr uint32_t kHostAnimationTag   = 0x414E494Du; // 'ANIM'
inline constexpr uint32_t kHostVoxelChunkTag  = 0x564F5843u; // 'VOXC'
inline constexpr uint32_t kHostSdfGraphTag    = 0x53444647u; // 'SDFG'

struct HostSdfGraph {
    bromesh::SdfGraph graph;
    uint32_t tag = kHostSdfGraphTag;
};

struct HostMesh {
    bromesh::MeshData mesh;
    uint32_t tag = kHostMeshTag;
};

struct HostMeshBVH {
    std::unique_ptr<bromesh::MeshBVH> bvh;
    bromesh::MeshData meshCopy;
    uint32_t tag = kHostMeshBVHTag;
};

struct HostProgressiveMesh {
    std::unique_ptr<bromesh::ProgressiveMesh> pm;
    uint32_t tag = kHostProgressiveTag;
};

struct HostSkinData {
    bromesh::SkinData skin;
    uint32_t tag = kHostSkinDataTag;
};

struct HostSkeleton {
    bromesh::Skeleton skeleton;
    uint32_t tag = kHostSkeletonTag;
};

struct HostJoint {
    bromesh::Bone bone;
    int index = -1;
    uint32_t tag = kHostJointTag;
};

struct HostSkeletonRig {
    bromesh::RigSpec spec;
    uint32_t tag = kHostSkeletonRigTag;
};

struct HostPose {
    bromesh::Pose pose;
    uint32_t tag = kHostPoseTag;
};

struct HostAnimation {
    bromesh::Animation animation;
    uint32_t tag = kHostAnimationTag;
};

struct HostVoxelChunk {
    std::unique_ptr<bromesh::VoxelChunk> chunk;
    uint32_t tag = kHostVoxelChunkTag;
};

// ---------------------------------------------------------------------------
// Host Classes
// ---------------------------------------------------------------------------
extern HostClass g_meshClass;
extern HostClass g_meshBvhClass;
extern HostClass g_progressiveMeshClass;
extern HostClass g_capsuleFieldClass;
extern HostClass g_lsystemClass;
extern HostClass g_polyMeshClass;
extern HostClass g_sdfGraphClass;

extern HostClass g_skinDataClass;
extern HostClass g_skeletonClass;
extern HostClass g_jointClass;
extern HostClass g_skeletonRigClass;
extern HostClass g_poseClass;
extern HostClass g_animationClass;
extern HostClass g_voxelChunkClass;

// ---------------------------------------------------------------------------
// Unwrap Helpers
// ---------------------------------------------------------------------------
inline HostMesh* unwrapMesh(Value v) {
    void* ptr = g_meshClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostMesh*>(ptr);
    return (h && h->tag == kHostMeshTag) ? h : nullptr;
}

inline HostMeshBVH* unwrapBVH(Value v) {
    void* ptr = g_meshBvhClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostMeshBVH*>(ptr);
    return (h && h->tag == kHostMeshBVHTag) ? h : nullptr;
}

inline HostProgressiveMesh* unwrapPM(Value v) {
    void* ptr = g_progressiveMeshClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostProgressiveMesh*>(ptr);
    return (h && h->tag == kHostProgressiveTag) ? h : nullptr;
}

inline HostSkinData* unwrapSkinData(Value v) {
    void* ptr = g_skinDataClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostSkinData*>(ptr);
    return (h && h->tag == kHostSkinDataTag) ? h : nullptr;
}

inline HostSkeleton* unwrapSkeleton(Value v) {
    void* ptr = g_skeletonClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostSkeleton*>(ptr);
    return (h && h->tag == kHostSkeletonTag) ? h : nullptr;
}

inline HostJoint* unwrapJoint(Value v) {
    void* ptr = g_jointClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostJoint*>(ptr);
    return (h && h->tag == kHostJointTag) ? h : nullptr;
}

inline HostSkeletonRig* unwrapSkeletonRig(Value v) {
    void* ptr = g_skeletonRigClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostSkeletonRig*>(ptr);
    return (h && h->tag == kHostSkeletonRigTag) ? h : nullptr;
}

inline HostPose* unwrapPose(Value v) {
    void* ptr = g_poseClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostPose*>(ptr);
    return (h && h->tag == kHostPoseTag) ? h : nullptr;
}

inline HostAnimation* unwrapAnimation(Value v) {
    void* ptr = g_animationClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostAnimation*>(ptr);
    return (h && h->tag == kHostAnimationTag) ? h : nullptr;
}

inline HostVoxelChunk* unwrapVoxelChunk(Value v) {
    void* ptr = g_voxelChunkClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostVoxelChunk*>(ptr);
    return (h && h->tag == kHostVoxelChunkTag) ? h : nullptr;
}

inline HostSdfGraph* unwrapSdfGraph(Value v) {
    void* ptr = g_sdfGraphClass.unwrap(v);
    if (!ptr) return nullptr;
    auto* h = static_cast<HostSdfGraph*>(ptr);
    return (h && h->tag == kHostSdfGraphTag) ? h : nullptr;
}

// ---------------------------------------------------------------------------
// Wrap Helpers
// ---------------------------------------------------------------------------
inline Value wrapMesh(bromesh::MeshData mesh) {
    auto h = std::make_unique<HostMesh>();
    h->mesh = std::move(mesh);
    return g_meshClass.createInstance(std::move(h));
}

inline Value wrapSkinData(bromesh::SkinData skin) {
    auto h = std::make_unique<HostSkinData>();
    h->skin = std::move(skin);
    return g_skinDataClass.createInstance(std::move(h));
}

inline Value wrapSkeleton(bromesh::Skeleton skel) {
    auto h = std::make_unique<HostSkeleton>();
    h->skeleton = std::move(skel);
    return g_skeletonClass.createInstance(std::move(h));
}

inline Value wrapPose(bromesh::Pose pose) {
    auto h = std::make_unique<HostPose>();
    h->pose = std::move(pose);
    return g_poseClass.createInstance(std::move(h));
}

inline Value wrapAnimation(bromesh::Animation anim) {
    auto h = std::make_unique<HostAnimation>();
    h->animation = std::move(anim);
    return g_animationClass.createInstance(std::move(h));
}

inline Value wrapJoint(bromesh::Bone bone, int index = -1) {
    auto h = std::make_unique<HostJoint>();
    h->bone = std::move(bone);
    h->index = index;
    return g_jointClass.createInstance(std::move(h));
}

inline Value wrapSkeletonRig(bromesh::RigSpec spec) {
    auto h = std::make_unique<HostSkeletonRig>();
    h->spec = std::move(spec);
    return g_skeletonRigClass.createInstance(std::move(h));
}

inline Value wrapSdfGraph(bromesh::SdfGraph graph) {
    auto h = std::make_unique<HostSdfGraph>();
    h->graph = std::move(graph);
    return g_sdfGraphClass.createInstance(std::move(h));
}

// ---------------------------------------------------------------------------
// TypedArray Creation and Conversion
// ---------------------------------------------------------------------------
inline Value makeFloat32Array(const float* data, size_t count) {
    Value arr = ev::createTypedArray(ev::elements::Float32, static_cast<uint32_t>(count));
    if (data && count > 0) {
        ev::fillTypedArray(arr, std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(data), count * sizeof(float)));
    }
    return arr;
}

inline Value makeUint32Array(const uint32_t* data, size_t count) {
    Value arr = ev::createTypedArray(ev::elements::Uint32, static_cast<uint32_t>(count));
    if (data && count > 0) {
        ev::fillTypedArray(arr, std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(data), count * sizeof(uint32_t)));
    }
    return arr;
}

inline Value makeUint8Array(const uint8_t* data, size_t count) {
    Value arr = ev::createTypedArray(ev::elements::Uint8, static_cast<uint32_t>(count));
    if (data && count > 0) {
        ev::fillTypedArray(arr, std::span<const uint8_t>(data, count));
    }
    return arr;
}

inline bool readFloat32Array(Value val, const float*& outData, size_t& outCount) {
    outData = nullptr;
    outCount = 0;
    ev::TypedArrayInfo info = ev::typedArrayInfo(val);
    if (!info || info.elementKind != ev::elements::Float32) return false;
    outData = reinterpret_cast<const float*>(info.data);
    outCount = info.elementCount;
    return true;
}

inline bool readUint32Array(Value val, const uint32_t*& outData, size_t& outCount) {
    outData = nullptr;
    outCount = 0;
    ev::TypedArrayInfo info = ev::typedArrayInfo(val);
    if (!info || info.elementKind != ev::elements::Uint32) return false;
    outData = reinterpret_cast<const uint32_t*>(info.data);
    outCount = info.elementCount;
    return true;
}

inline bool readUint8Array(Value val, const uint8_t*& outData, size_t& outCount) {
    outData = nullptr;
    outCount = 0;
    ev::TypedArrayInfo info = ev::typedArrayInfo(val);
    if (!info || (info.elementKind != ev::elements::Uint8 && info.elementKind != ev::elements::Uint8Clamped)) {
        return false;
    }
    outData = reinterpret_cast<const uint8_t*>(info.data);
    outCount = info.elementCount;
    return true;
}

inline std::vector<float> toFloatVector(Value val) {
    const float* fptr = nullptr;
    size_t fcount = 0;
    if (readFloat32Array(val, fptr, fcount)) {
        return std::vector<float>(fptr, fptr + fcount);
    }
    std::vector<float> result;
    if (ev::isObject(val)) {
        Value lenVal = ev::getProperty(val, "length");
        if (ev::isNumber(lenVal)) {
            size_t n = static_cast<size_t>(ev::toDouble(lenVal));
            result.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                Value elem = ev::getElement(val, static_cast<uint32_t>(i));
                if (ev::isNumber(elem)) {
                    result.push_back(static_cast<float>(ev::toDouble(elem)));
                } else if (ev::isObject(elem)) {
                    // Possible [x, y, z] or {x, y, z}
                    Value x = ev::getProperty(elem, "x");
                    Value y = ev::getProperty(elem, "y");
                    Value z = ev::getProperty(elem, "z");
                    if (!ev::isUndefined(x) && !ev::isUndefined(y) && !ev::isUndefined(z)) {
                        result.push_back(static_cast<float>(ev::toDouble(x)));
                        result.push_back(static_cast<float>(ev::toDouble(y)));
                        result.push_back(static_cast<float>(ev::toDouble(z)));
                    } else {
                        Value e0 = ev::getElement(elem, 0);
                        Value e1 = ev::getElement(elem, 1);
                        Value e2 = ev::getElement(elem, 2);
                        if (!ev::isUndefined(e0) && !ev::isUndefined(e1) && !ev::isUndefined(e2)) {
                            result.push_back(static_cast<float>(ev::toDouble(e0)));
                            result.push_back(static_cast<float>(ev::toDouble(e1)));
                            result.push_back(static_cast<float>(ev::toDouble(e2)));
                        }
                    }
                }
            }
        }
    }
    return result;
}

inline std::vector<uint32_t> toUint32Vector(Value val) {
    const uint32_t* uptr = nullptr;
    size_t ucount = 0;
    if (readUint32Array(val, uptr, ucount)) {
        return std::vector<uint32_t>(uptr, uptr + ucount);
    }
    std::vector<uint32_t> result;
    if (ev::isObject(val)) {
        Value lenVal = ev::getProperty(val, "length");
        if (ev::isNumber(lenVal)) {
            size_t n = static_cast<size_t>(ev::toDouble(lenVal));
            result.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                Value elem = ev::getElement(val, static_cast<uint32_t>(i));
                result.push_back(static_cast<uint32_t>(ev::toDouble(elem)));
            }
        }
    }
    return result;
}

inline std::vector<uint8_t> toUint8Vector(Value val) {
    const uint8_t* uptr = nullptr;
    size_t ucount = 0;
    if (readUint8Array(val, uptr, ucount)) {
        return std::vector<uint8_t>(uptr, uptr + ucount);
    }
    std::vector<uint8_t> result;
    if (ev::isObject(val)) {
        Value lenVal = ev::getProperty(val, "length");
        if (ev::isNumber(lenVal)) {
            size_t n = static_cast<size_t>(ev::toDouble(lenVal));
            result.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                Value elem = ev::getElement(val, static_cast<uint32_t>(i));
                result.push_back(static_cast<uint8_t>(ev::toDouble(elem)));
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Module Installation Entry Points
// ---------------------------------------------------------------------------
void ensureMeshClassesInstalled();
void ensureRiggingClassesInstalled();

void initMeshCore(ObjectBuilder& proto, HostClass& cls);
void initMeshOps(ObjectBuilder& proto, HostClass& cls);
void initMeshAnalysis(ObjectBuilder& proto, HostClass& cls);

void initMeshPlants(HostClass& cls);
void initMeshIo(ObjectBuilder& proto, HostClass& cls);

void initMeshBvh(HostClass& cls);
void initProgressiveMesh(HostClass& cls);
void initCapsuleField(HostClass& cls);
void initLSystem(HostClass& cls);
void initPolyMesh(HostClass& cls);
void initMeshSdf(ObjectBuilder& proto, HostClass& cls);
void initSdfGraph(HostClass& cls);

// A file path as the host's resolver sees it (api.h setPathResolver), or as
// given when no resolver is set (native_mesh_io.cpp).
std::string resolveMeshPath(const std::string& path);
// The same for a file about to be written: its parent directory is what
// gets resolved, since the file itself does not exist yet.
std::string resolveMeshWritePath(const std::string& path);

void initRiggingCore(HostClass& skinCls, HostClass& skelCls, HostClass& jointCls,
                     HostClass& rigCls, HostClass& voxelCls);
void initRiggingAnim(HostClass& poseCls, HostClass& animCls, HostClass& meshCls,
                     ObjectBuilder& ikBuilder);

Value makeMeshNamespace();
Value makeRiggingNamespace();

} // namespace bromesh::api
