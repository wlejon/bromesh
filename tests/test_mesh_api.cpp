#include "bromesh/api.h"
#include "eval/eval.h"
#include "embed/embed.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

namespace ev = bronze::embed;
using bronze::Value;

// assert() is compiled out of the Release build this test runs in.
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::cerr << "FAIL: " #cond " (" __FILE__ ":" << __LINE__ << ")\n"; \
            return 1;                                                            \
        }                                                                        \
    } while (0)

// tests/test_mesh_api_restored.cpp — the Mesh / Skeleton members the bronze
// port dropped (bro docs/transition-drift.md row H7).
void bromeshTestRestoredSurface();
// Same file: the bindings fixed in the GC-rooting audit.
void bromeshTestBindingFixes();
// tests/test_mesh_api_validation.cpp — counts, sizes and indices throw.
void bromeshTestValidation();

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
    // Every Value that outlives an allocating call is held in a Persistent
    // (embed.h's GC contract), so the test is itself clean under
    // BRONZE_GC_STRESS=1.
    ev::Persistent bro(ev::globalValue("bro").value);
    CHECK(ev::isObject(bro.get()));
    {
        Value v = ev::getProperty(bro.get(), "mesh");
        CHECK(ev::isObject(v));
        v = ev::getProperty(bro.get(), "rigging");
        CHECK(ev::isObject(v));
    }

    // Global constructors (IK is a namespace object, not a constructor).
    for (const char* name : {"Mesh", "MeshBVH", "ProgressiveMesh", "Skeleton", "Joint",
                             "SkeletonRig", "Pose", "AnimationClip", "SkinData"}) {
        auto g = ev::globalValue(name);
        if (!g.found || !ev::isFunction(g.value)) {
            std::cerr << "FAIL: global " << name << " is not a constructor" << std::endl;
            return 1;
        }
    }
    {
        auto g = ev::globalValue("IK");
        CHECK(g.found && ev::isObject(g.value));
    }

    std::cout << "  Mounting verification passed." << std::endl;

    // 3. Direct embed API calls
    std::cout << "[3/4] Testing direct embed API calls..." << std::endl;
    ev::Persistent meshCtor(ev::globalValue("Mesh").value);
    ev::Persistent boxFn(ev::getProperty(meshCtor.get(), "box"));
    CHECK(ev::isFunction(boxFn.get()));

    const Value boxArgs[3] = {ev::fromDouble(1.0), ev::fromDouble(1.0), ev::fromDouble(1.0)};
    auto boxRes = ev::call(boxFn.get(), meshCtor.get(), std::span<const Value>(boxArgs, 3));
    CHECK(!boxRes.thrown);
    CHECK(ev::isObject(boxRes.value));
    ev::Persistent boxMesh(boxRes.value);

    Value vCount = ev::getProperty(boxMesh.get(), "vertexCount");
    CHECK(ev::isNumber(vCount) && ev::toDouble(vCount) == 24.0);

    Value tCount = ev::getProperty(boxMesh.get(), "triangleCount");
    CHECK(ev::isNumber(tCount) && ev::toDouble(tCount) == 12.0);

    ev::Persistent isManifold(ev::getProperty(boxMesh.get(), "isManifold"));
    CHECK(ev::isFunction(isManifold.get()));
    auto maniRes = ev::call(isManifold.get(), boxMesh.get(), {});
    CHECK(!maniRes.thrown && ev::toBool(maniRes.value));

    ev::Persistent volFn(ev::getProperty(boxMesh.get(), "volume"));
    CHECK(ev::isFunction(volFn.get()));
    auto volRes = ev::call(volFn.get(), boxMesh.get(), {});
    CHECK(!volRes.thrown && ev::isNumber(volRes.value));
    const double vol = ev::toDouble(volRes.value);
    CHECK(std::fabs(vol - 8.0) < 1e-2);

    std::cout << "  Direct embed calls passed (box volume: " << vol << ")." << std::endl;

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
        if (!b || !b.min || b.min.length !== 3 || !b.max || b.max.length !== 3) throw new Error("Bounds should have min and max with 3 components");

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
        const ikRes = IK.twoBone(skel, pose, 0, 1, 2, [1, 1, 0]);
        if (typeof ikRes !== "boolean") throw new Error("IK twoBone failed");

        // VoxelChunk
        const chunk = new VoxelChunk(8, 8, 8, 1.0);
        chunk.set(2, 2, 2, 1);
        if (chunk.get(2, 2, 2) !== 1) throw new Error("VoxelChunk get/set mismatch");

        // SDFGraph and JIT SDF meshing
        const sdf = Mesh.createSDF();
        const s0 = sdf.sphere(1.0);
        const b0 = sdf.box([0.7, 0.7, 0.7]);
        const u0 = sdf.opSmoothUnion(s0, b0, 0.2);
        sdf.setRoot(u0);
        if (sdf.size !== 3) throw new Error("SDFGraph size should be 3");

        const sdfMesh = sdf.marchingCubes({ dims: [32, 32, 32], bounds: { min: [-1.5, -1.5, -1.5], max: [1.5, 1.5, 1.5] } });
        if (!sdfMesh || sdfMesh.triangleCount === 0) throw new Error("SDF marching cubes returned empty mesh");
        if (sdfMesh.volume() <= 0) throw new Error("SDF marching cubes volume must be positive");

        const snMesh = Mesh.surfaceNetsSDF(sdf, { dims: [32, 32, 32], bounds: { min: [-1.5, -1.5, -1.5], max: [1.5, 1.5, 1.5] } });
        if (!snMesh || snMesh.triangleCount === 0) throw new Error("SDF surface nets returned empty mesh");
        if (snMesh.volume() <= 0) throw new Error("SDF surface nets volume must be positive");

        "SUCCESS";
    )JS";

    auto evalRes = bronze::eval::evalScript(testScript);
    if (evalRes.thrown) {
        std::cerr << "Eval threw error: " << ev::toUtf8(evalRes.value) << std::endl;
        return 1;
    }
    const std::string evalOut = ev::toUtf8(evalRes.value);
    std::cout << "  Bronze eval returned: " << evalOut << std::endl;
    CHECK(evalOut == "SUCCESS");

    bromeshTestRestoredSurface();
    bromeshTestBindingFixes();
    bromeshTestValidation();

    std::cout << "========================================" << std::endl;
    std::cout << "All bromesh Bronze API tests PASSED!" << std::endl;
    std::cout << "========================================" << std::endl;
    return 0;
}
