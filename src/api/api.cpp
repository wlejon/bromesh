#include "api.h"
#include "host_mesh_internal.h"

namespace bromesh::api {

namespace {

// globalThis (undefined when the realm has none), rooted: every step of an
// install allocates, and a raw Value is only good until the next allocation.
ev::Persistent rootedGlobalThis() {
    ev::GlobalValue gt = ev::globalValue("globalThis");
    return ev::Persistent(gt.found && ev::isObject(gt.value) ? gt.value : ev::undefined());
}

// The realm's `bro` object, created (and published on globalThis) when no
// earlier installer made it.
ev::Persistent rootedBro(const ev::Persistent& globalThis) {
    {
        ev::GlobalValue g = ev::globalValue("bro");
        if (g.found && ev::isObject(g.value)) return ev::Persistent(g.value);
    }
    if (ev::isObject(globalThis.get())) {
        Value candidate = ev::getProperty(globalThis.get(), "bro");
        if (ev::isObject(candidate)) return ev::Persistent(candidate);
    }
    ev::Persistent bro(ev::createObject());
    ev::registerGlobal("bro", bro.get());
    if (ev::isObject(globalThis.get())) ev::setProperty(globalThis.get(), "bro", bro.get());
    return bro;
}

// Publish a constructor as a global binding and on globalThis.
void publishGlobal(const ev::Persistent& globalThis, const char* name, const HostClass& cls) {
    if (ev::isObject(globalThis.get())) ev::setProperty(globalThis.get(), name, cls.constructor());
    ev::registerGlobal(name, cls.constructor());
}

}  // namespace

Value makeMeshNamespace() {
    ObjectBuilder ns;
    ns.set("Mesh", g_meshClass.constructor());
    ns.set("MeshBVH", g_meshBvhClass.constructor());
    ns.set("ProgressiveMesh", g_progressiveMeshClass.constructor());
    ns.set("CapsuleField", g_capsuleFieldClass.constructor());
    ns.set("LSystem", g_lsystemClass.constructor());
    ns.set("PolyMesh", g_polyMeshClass.constructor());
    ns.set("SDFGraph", g_sdfGraphClass.constructor());

    // Forward static methods of Mesh onto bro.mesh
    const char* factories[] = {
        "box", "sphere", "cylinder", "capsule", "cone", "plane", "torus",
        "icosahedron", "dodecahedron", "octahedron", "tetrahedron", "disk",
#if BROMESH_HAS_PAR_SHAPES
        "geodesicSphere", "rock",
#endif
        "blob", "tube", "sweep", "bezierSweep", "heightmapGrid",
        "leafCard", "flower", "bladeStrip", "bladePath",
        "spaceColonize", "thickenBranches", "meshBranches",
        "placeLeavesOnBranches", "scatterLeaves", "tree",
        "capsuleField", "capsuleFieldFromSegments", "packAnchors",
        "parseLSystem", "lsystemToBranches",
        "merge", "booleanUnion", "booleanDifference", "booleanIntersection",
        "convexHull", "marchingCubes", "surfaceNets", "dualContouring",
        "transvoxel", "greedyMesh", "polygon2D", "polygon3D", "reconstruct",
        "createSDF", "marchingCubesSDF", "surfaceNetsSDF",
        "loadOBJ", "loadPLY", "loadSTL", "loadVOX", "loadFBX", "loadSplatPLY", "saveSplatPLY"
#if BROMESH_HAS_GLTF
        , "loadGLTF"
#endif
#if BROMESH_HAS_DRACO
        , "decodeDraco", "encodeDraco"
#endif
    };
    for (const char* f : factories) {
        // Re-fetch the constructor each round: ns.set allocates, and a raw
        // Value is only good until the next allocation.
        Value fn = ev::getProperty(g_meshClass.constructor(), f);
        if (ev::isFunction(fn)) {
            ns.set(f, fn);
        }
    }
    return ns.build();
}

Value makeRiggingNamespace() {
    ObjectBuilder ns;
    ns.set("SkinData", g_skinDataClass.constructor());
    ns.set("Skeleton", g_skeletonClass.constructor());
    ns.set("Joint", g_jointClass.constructor());
    ns.set("SkeletonRig", g_skeletonRigClass.constructor());
    ns.set("RigSpec", g_skeletonRigClass.constructor());
    ns.set("Rig", g_skeletonRigClass.constructor());
    ns.set("Pose", g_poseClass.constructor());
    ns.set("AnimationClip", g_animationClass.constructor());
    ns.set("Animation", g_animationClass.constructor());
    ns.set("SkeletalAnimation", g_animationClass.constructor());
    ns.set("VoxelChunk", g_voxelChunkClass.constructor());

    auto ikGlobal = ev::globalValue("IK");
    if (ikGlobal.found && ev::isObject(ikGlobal.value)) {
        ns.set("IK", ikGlobal.value);
    }

    const char* rigMethods[] = {
        "specFromFile", "detectHumanoid", "detectLandmarks", "detectQuadruped",
        "missingLandmarks", "fitSkeleton", "autoRig", "transferWeights"
    };
    for (const char* m : rigMethods) {
        // Re-read the constructor from its root each round (ns.set allocates).
        Value fn = ev::getProperty(g_skeletonRigClass.constructor(), m);
        if (ev::isFunction(fn)) {
            ns.set(m, fn);
        }
    }
    return ns.build();
}

void installMesh() {
    ensureMeshClassesInstalled();

    ev::Persistent globalThis = rootedGlobalThis();
    ev::Persistent bro = rootedBro(globalThis);

    // Mount bro.mesh
    Value meshVal = makeMeshNamespace();
    bro.set(ev::setProperty(bro.get(), "mesh", meshVal));

    publishGlobal(globalThis, "Mesh", g_meshClass);
    publishGlobal(globalThis, "MeshBVH", g_meshBvhClass);
    publishGlobal(globalThis, "ProgressiveMesh", g_progressiveMeshClass);
    publishGlobal(globalThis, "CapsuleField", g_capsuleFieldClass);
    publishGlobal(globalThis, "LSystem", g_lsystemClass);
    publishGlobal(globalThis, "PolyMesh", g_polyMeshClass);
    publishGlobal(globalThis, "SDFGraph", g_sdfGraphClass);
}

bool isMeshValue(Value v) {
    return unwrapMesh(v) != nullptr;
}

bool takeMeshData(Value v, bromesh::MeshData& out) {
    HostMesh* h = unwrapMesh(v);
    if (!h) return false;
    out = std::move(h->mesh);
    h->mesh = bromesh::MeshData();
    return true;
}

Value makeMeshValue(bromesh::MeshData mesh) {
    ensureMeshClassesInstalled();
    return wrapMesh(std::move(mesh));
}

const bromesh::MeshData* meshDataOf(Value v) {
    HostMesh* h = unwrapMesh(v);
    return h ? &h->mesh : nullptr;
}

bromesh::MeshData* meshDataOfMut(Value v) {
    HostMesh* h = unwrapMesh(v);
    return h ? &h->mesh : nullptr;
}

const bromesh::SkinData* skinDataOf(Value v) {
    HostSkinData* h = unwrapSkinData(v);
    return h ? &h->skin : nullptr;
}

const bromesh::Skeleton* skeletonOf(Value v) {
    HostSkeleton* h = unwrapSkeleton(v);
    return h ? &h->skeleton : nullptr;
}

const bromesh::Animation* animationOf(Value v) {
    HostAnimation* h = unwrapAnimation(v);
    return h ? &h->animation : nullptr;
}

void installRigging() {
    ensureRiggingClassesInstalled();

    ev::Persistent globalThis = rootedGlobalThis();
    ev::Persistent broP = rootedBro(globalThis);

    // Mount bro.rigging
    Value rigVal = makeRiggingNamespace();
    broP.set(ev::setProperty(broP.get(), "rigging", rigVal));

    // Mesh.loadGLTF is a rigging-side static (it returns skins and
    // skeletons), so it exists only now — after installMesh built bro.mesh.
    // Forward it late, the way the mesh-side statics were forwarded early.
#if BROMESH_HAS_GLTF
    {
        Value meshNs = ev::getProperty(broP.get(), "mesh");
        if (ev::isObject(meshNs)) {
            ev::Persistent meshNsP(meshNs);
            Value fn = ev::getProperty(g_meshClass.constructor(), "loadGLTF");
            if (ev::isFunction(fn)) ev::setProperty(meshNsP.get(), "loadGLTF", fn);
        }
    }
#endif

    publishGlobal(globalThis, "SkinData", g_skinDataClass);
    publishGlobal(globalThis, "Skeleton", g_skeletonClass);
    publishGlobal(globalThis, "Joint", g_jointClass);
    publishGlobal(globalThis, "SkeletonRig", g_skeletonRigClass);
    publishGlobal(globalThis, "RigSpec", g_skeletonRigClass);
    publishGlobal(globalThis, "Rig", g_skeletonRigClass);
    publishGlobal(globalThis, "Pose", g_poseClass);
    publishGlobal(globalThis, "AnimationClip", g_animationClass);
    publishGlobal(globalThis, "Animation", g_animationClass);
    publishGlobal(globalThis, "SkeletalAnimation", g_animationClass);
    publishGlobal(globalThis, "VoxelChunk", g_voxelChunkClass);
}

} // namespace bromesh::api
