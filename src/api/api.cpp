#include "api.h"
#include "host_mesh_internal.h"

namespace bromesh::api {

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

    Value rigCtor = g_skeletonRigClass.constructor();
    const char* rigMethods[] = {
        "specFromFile", "detectHumanoid", "detectLandmarks", "detectQuadruped",
        "missingLandmarks", "fitSkeleton", "autoRig", "transferWeights"
    };
    for (const char* m : rigMethods) {
        Value fn = ev::getProperty(rigCtor, m);
        if (ev::isFunction(fn)) {
            ns.set(m, fn);
        }
    }
    return ns.build();
}

void installMesh() {
    ensureMeshClassesInstalled();

    Value globalThisVal = ev::undefined();
    auto gt = ev::globalValue("globalThis");
    if (gt.found && ev::isObject(gt.value)) {
        globalThisVal = gt.value;
    }

    Value broVal = ev::globalValue("bro").found ? ev::globalValue("bro").value : ev::undefined();
    if (!ev::isObject(broVal)) {
        if (!ev::isUndefined(globalThisVal)) {
            Value candidate = ev::getProperty(globalThisVal, "bro");
            if (ev::isObject(candidate)) {
                broVal = candidate;
            }
        }
    }
    if (!ev::isObject(broVal)) {
        broVal = ev::createObject();
        ev::registerGlobal("bro", broVal);
        if (!ev::isUndefined(globalThisVal)) {
            ev::setProperty(globalThisVal, "bro", broVal);
        }
    }

    ev::Persistent broP(broVal);

    // Mount bro.mesh
    Value meshVal = makeMeshNamespace();
    broP.set(ev::setProperty(broP.get(), "mesh", meshVal));

    if (!ev::isUndefined(globalThisVal)) {
        ev::setProperty(globalThisVal, "Mesh", g_meshClass.constructor());
        ev::setProperty(globalThisVal, "MeshBVH", g_meshBvhClass.constructor());
        ev::setProperty(globalThisVal, "ProgressiveMesh", g_progressiveMeshClass.constructor());
        ev::setProperty(globalThisVal, "CapsuleField", g_capsuleFieldClass.constructor());
        ev::setProperty(globalThisVal, "LSystem", g_lsystemClass.constructor());
        ev::setProperty(globalThisVal, "PolyMesh", g_polyMeshClass.constructor());
        ev::setProperty(globalThisVal, "SDFGraph", g_sdfGraphClass.constructor());
    }
    ev::registerGlobal("Mesh", g_meshClass.constructor());
    ev::registerGlobal("MeshBVH", g_meshBvhClass.constructor());
    ev::registerGlobal("ProgressiveMesh", g_progressiveMeshClass.constructor());
    ev::registerGlobal("CapsuleField", g_capsuleFieldClass.constructor());
    ev::registerGlobal("LSystem", g_lsystemClass.constructor());
    ev::registerGlobal("PolyMesh", g_polyMeshClass.constructor());
    ev::registerGlobal("SDFGraph", g_sdfGraphClass.constructor());
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

void installRigging() {
    ensureRiggingClassesInstalled();

    Value globalThisVal = ev::undefined();
    auto gt = ev::globalValue("globalThis");
    if (gt.found && ev::isObject(gt.value)) {
        globalThisVal = gt.value;
    }

    Value broVal = ev::globalValue("bro").found ? ev::globalValue("bro").value : ev::undefined();
    if (!ev::isObject(broVal)) {
        if (!ev::isUndefined(globalThisVal)) {
            Value candidate = ev::getProperty(globalThisVal, "bro");
            if (ev::isObject(candidate)) {
                broVal = candidate;
            }
        }
    }
    if (!ev::isObject(broVal)) {
        broVal = ev::createObject();
        ev::registerGlobal("bro", broVal);
        if (!ev::isUndefined(globalThisVal)) {
            ev::setProperty(globalThisVal, "bro", broVal);
        }
    }

    ev::Persistent broP(broVal);

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

    if (!ev::isUndefined(globalThisVal)) {
        ev::setProperty(globalThisVal, "SkinData", g_skinDataClass.constructor());
        ev::setProperty(globalThisVal, "Skeleton", g_skeletonClass.constructor());
        ev::setProperty(globalThisVal, "Joint", g_jointClass.constructor());
        ev::setProperty(globalThisVal, "SkeletonRig", g_skeletonRigClass.constructor());
        ev::setProperty(globalThisVal, "RigSpec", g_skeletonRigClass.constructor());
        ev::setProperty(globalThisVal, "Rig", g_skeletonRigClass.constructor());
        ev::setProperty(globalThisVal, "Pose", g_poseClass.constructor());
        ev::setProperty(globalThisVal, "AnimationClip", g_animationClass.constructor());
        ev::setProperty(globalThisVal, "Animation", g_animationClass.constructor());
        ev::setProperty(globalThisVal, "SkeletalAnimation", g_animationClass.constructor());
        ev::setProperty(globalThisVal, "VoxelChunk", g_voxelChunkClass.constructor());
    }
    ev::registerGlobal("SkinData", g_skinDataClass.constructor());
    ev::registerGlobal("Skeleton", g_skeletonClass.constructor());
    ev::registerGlobal("Joint", g_jointClass.constructor());
    ev::registerGlobal("SkeletonRig", g_skeletonRigClass.constructor());
    ev::registerGlobal("RigSpec", g_skeletonRigClass.constructor());
    ev::registerGlobal("Rig", g_skeletonRigClass.constructor());
    ev::registerGlobal("Pose", g_poseClass.constructor());
    ev::registerGlobal("AnimationClip", g_animationClass.constructor());
    ev::registerGlobal("Animation", g_animationClass.constructor());
    ev::registerGlobal("SkeletalAnimation", g_animationClass.constructor());
    ev::registerGlobal("VoxelChunk", g_voxelChunkClass.constructor());
}

} // namespace bromesh::api
