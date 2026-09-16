#include "bromesh/api.h"
#include "eval/eval.h"
#include "embed/embed.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

namespace ev = bronze::embed;
using bronze::Value;

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Running bromesh Standalone Bronze API Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    // 1. Install Mesh and Rigging APIs into Bronze realm
    std::cout << "[1/4] Installing Mesh and Rigging APIs..." << std::endl;
    bromesh::api::installMesh();
    bromesh::api::installRigging();

    // 2. Verify global mountings
    std::cout << "[2/4] Verifying global mountings..." << std::endl;
    auto gBro = ev::globalValue("bro");
    assert(gBro.found);
    assert(ev::isObject(gBro.value));

    auto meshNs = ev::getProperty(gBro.value, "mesh");
    assert(ev::isObject(meshNs));

    auto rigNs = ev::getProperty(gBro.value, "rigging");
    assert(ev::isObject(rigNs));

    // Verify global constructors
    auto gMesh = ev::globalValue("Mesh");
    assert(gMesh.found && ev::isFunction(gMesh.value));

    auto gBVH = ev::globalValue("MeshBVH");
    assert(gBVH.found && ev::isFunction(gBVH.value));

    auto gPM = ev::globalValue("ProgressiveMesh");
    assert(gPM.found && ev::isFunction(gPM.value));

    auto gSkel = ev::globalValue("Skeleton");
    assert(gSkel.found && ev::isFunction(gSkel.value));

    auto gJoint = ev::globalValue("Joint");
    assert(gJoint.found && ev::isFunction(gJoint.value));

    auto gRig = ev::globalValue("SkeletonRig");
    assert(gRig.found && ev::isFunction(gRig.value));

    auto gPose = ev::globalValue("Pose");
    assert(gPose.found && ev::isFunction(gPose.value));

    auto gAnim = ev::globalValue("AnimationClip");
    assert(gAnim.found && ev::isFunction(gAnim.value));

    auto gSkin = ev::globalValue("SkinData");
    assert(gSkin.found && ev::isFunction(gSkin.value));

    auto gIK = ev::globalValue("IK");
    assert(gIK.found && ev::isObject(gIK.value));

    std::cout << "  Mounting verification passed." << std::endl;

    // 3. Direct embed API calls
    std::cout << "[3/4] Testing direct embed API calls..." << std::endl;
    // Call Mesh.box(1, 1, 1)
    auto boxFn = ev::getProperty(gMesh.value, "box");
    assert(ev::isFunction(boxFn));

    const Value boxArgs[3] = {ev::fromDouble(1.0), ev::fromDouble(1.0), ev::fromDouble(1.0)};
    auto boxRes = ev::call(boxFn, gMesh.value, std::span<const Value>(boxArgs, 3));
    assert(!boxRes.thrown);
    assert(ev::isObject(boxRes.value));

    Value boxMesh = boxRes.value;
    Value vCount = ev::getProperty(boxMesh, "vertexCount");
    assert(ev::isNumber(vCount));
    assert(ev::toDouble(vCount) == 24.0);

    Value tCount = ev::getProperty(boxMesh, "triangleCount");
    assert(ev::isNumber(tCount));
    assert(ev::toDouble(tCount) == 12.0);

    Value isManifold = ev::getProperty(boxMesh, "isManifold");
    assert(ev::isFunction(isManifold));
    auto maniRes = ev::call(isManifold, boxMesh, {});
    assert(!maniRes.thrown && ev::toBool(maniRes.value));

    Value volFn = ev::getProperty(boxMesh, "volume");
    assert(ev::isFunction(volFn));
    auto volRes = ev::call(volFn, boxMesh, {});
    assert(!volRes.thrown && ev::isNumber(volRes.value));
    assert(std::fabs(ev::toDouble(volRes.value) - 8.0) < 1e-2);

    std::cout << "  Direct embed calls passed (box volume: " << ev::toDouble(volRes.value) << ")." << std::endl;

    // 4. Test Bronze evalScript
    std::cout << "[4/4] Testing JS execution via Bronze eval..." << std::endl;
    const char* testScript = R"JS(
        // Mesh creation and transforms
        const m = Mesh.box(0.5, 0.5, 0.5);
        if (m.vertexCount !== 24) throw new Error("Expected 24 vertices");
        if (m.triangleCount !== 12) throw new Error("Expected 12 triangles");
        if (!m.isManifold()) throw new Error("Box should be manifold");

        m.translate(1, 2, 3);
        const b = m.bounds();
        if (b.length !== 6) throw new Error("Bounds should have 6 components");

        // MeshBVH raycasting
        const bvh = new MeshBVH(m);
        if (bvh.empty) throw new Error("BVH should not be empty");
        const hit = bvh.raycast([1, 10, 3], [0, -1, 0], 100);
        if (!hit || !hit.hit) throw new Error("Ray should hit box");

        // Clone and merge operation
        const m2 = m.clone();
        m2.translate(0.5, 0, 0);
        const merged = Mesh.merge([m, m2]);
        if (merged.triangleCount !== m.triangleCount * 2) throw new Error("Mesh merge failed");

        // Boolean union operation (returns valid mesh; non-empty if manifold is compiled in)
        const unionMesh = Mesh.booleanUnion(m, m2);
        if (!unionMesh || typeof unionMesh.triangleCount !== "number") throw new Error("Boolean union returned invalid mesh");

        // Skeleton and Pose
        const skel = new Skeleton({
            bones: [
                { name: "root", parent: -1, localT: [0, 0, 0] },
                { name: "spine", parent: 0, localT: [0, 1, 0] },
                { name: "head", parent: 1, localT: [0, 1, 0] }
            ]
        });
        if (skel.boneCount !== 3) throw new Error("Skeleton should have 3 bones");
        if (skel.findBone("spine") !== 1) throw new Error("Bone index mismatch");

        const pose = skel.bindPose();
        if (pose.boneCount !== 3) throw new Error("Pose should have 3 bones");

        const mats = pose.computeWorldMatrices(skel);
        if (mats.length !== 3 * 16) throw new Error("World matrices length mismatch");

        // SkinData
        const skin = new SkinData({
            boneWeights: new Float32Array([1, 0, 0, 0, 0.5, 0.5, 0, 0]),
            boneIndices: new Uint32Array([0, 0, 0, 0, 0, 1, 0, 0]),
            boneCount: 2
        });
        skin.normalize();
        if (skin.boneCount !== 2) throw new Error("Skin boneCount mismatch");

        // IK
        const ikRes = IK.twoBone({
            rootPos: [0, 0, 0],
            midPos: [0, 1, 0],
            endPos: [0, 2, 0],
            targetPos: [1, 1, 0]
        });
        if (!ikRes || !ikRes.midPos) throw new Error("IK twoBone failed");

        // VoxelChunk
        const chunk = new VoxelChunk(8, 8, 8, 1.0);
        chunk.set(2, 2, 2, 1);
        if (chunk.get(2, 2, 2) !== 1) throw new Error("VoxelChunk get/set mismatch");

        "SUCCESS";
    )JS";

    auto evalRes = bronze::eval::evalScript(testScript);
    if (evalRes.thrown) {
        std::cerr << "Eval threw error: " << ev::toUtf8(evalRes.value) << std::endl;
        return 1;
    }
    std::cout << "  Bronze eval returned: " << ev::toUtf8(evalRes.value) << std::endl;
    assert(ev::toUtf8(evalRes.value) == "SUCCESS");

    std::cout << "========================================" << std::endl;
    std::cout << "All bromesh Bronze API tests PASSED!" << std::endl;
    std::cout << "========================================" << std::endl;
    return 0;
}
