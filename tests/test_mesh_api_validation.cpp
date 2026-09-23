// Argument validation in the bromesh JS bindings. Counts, sizes, segments,
// resolutions and indices used to go through a bare static_cast of the JS
// double, so -1 wrapped to ~1.8e19 (then sized an allocation) and NaN or a
// fraction was undefined behaviour. They now throw: TypeError for a
// non-number, RangeError for NaN, a fraction or a value out of range.
#include "eval/eval.h"
#include "embed/embed.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {
namespace ev = bronze::embed;
}

void bromeshTestValidation() {
    std::cout << "[validation] Exercising count / index validation..." << std::endl;

    const char* script = R"JS(
        function expectThrow(kind, label, fn) {
            let err = null;
            try { fn(); } catch (e) { err = e; }
            if (err === null) throw new Error(label + ": expected a " + kind + ", nothing was thrown");
            if (err.name !== kind) throw new Error(label + ": expected a " + kind + ", got " + err.name + ": " + err.message);
        }

        // ── primitives: negative, NaN, fractional, non-number ──────────────
        expectThrow("RangeError", "sphere segments -1", () => Mesh.sphere(1, -1));
        expectThrow("RangeError", "sphere segments NaN", () => Mesh.sphere(1, NaN));
        expectThrow("RangeError", "sphere segments 2.5", () => Mesh.sphere(1, 2.5));
        expectThrow("RangeError", "sphere segments 2 (< 3)", () => Mesh.sphere(1, 2));
        expectThrow("RangeError", "sphere segments 1e12", () => Mesh.sphere(1, 1e12));
        expectThrow("TypeError", "sphere segments '8'", () => Mesh.sphere(1, "8"));
        expectThrow("RangeError", "torus tubeSegments -3", () => Mesh.torus(1, 0.25, 16, -3));
        expectThrow("RangeError", "geodesicSphere subdivisions 40", () => Mesh.geodesicSphere(1, 40));
        if (Mesh.sphere(1, 8, 6).triangleCount === 0) throw new Error("a valid sphere is empty");

        // ── ops: iteration counts, 0 is a no-op ─────────────────────────────
        const cube = Mesh.box(1, 1, 1);
        const tris = cube.triangleCount;
        cube.subdivideMidpoint(0);
        if (cube.triangleCount !== tris) throw new Error("subdivideMidpoint(0) should be a no-op");
        expectThrow("RangeError", "subdivideLoop -1", () => Mesh.box(1, 1, 1).subdivideLoop(-1));
        expectThrow("RangeError", "subdivideLoop 9", () => Mesh.box(1, 1, 1).subdivideLoop(9));
        expectThrow("RangeError", "smooth iterations -2", () => Mesh.box(1, 1, 1).smooth(0.5, -2));
        expectThrow("RangeError", "fillHoles -1", () => Mesh.box(1, 1, 1).fillHoles(-1));
        expectThrow("RangeError", "bakeAOToTexture width -512", () => Mesh.box(1, 1, 1).bakeAOToTexture(-512, 64, 8));
        expectThrow("RangeError", "bakeAmbientOcclusion rays 0", () => Mesh.box(1, 1, 1).bakeAmbientOcclusion(0));
        expectThrow("RangeError", "shrinkwrap mode 7", () => Mesh.box(1, 1, 1).shrinkwrap(Mesh.box(2, 2, 2), 7));
        expectThrow("RangeError", "convexDecomposition maxHulls -1",
                    () => Mesh.box(1, 1, 1).convexDecomposition({ maxHulls: -1 }));

        // ── analysis: meshlets, sampling, transvoxel lod ────────────────────
        expectThrow("RangeError", "sampleSurface -5", () => Mesh.box(1, 1, 1).sampleSurface(-5));
        expectThrow("RangeError", "buildMeshlets maxVertices 1024",
                    () => Mesh.box(1, 1, 1).buildMeshlets({ maxVertices: 1024 }));
        const G = 9;
        const field = new Float32Array(G * G * G);
        for (let z = 0; z < G; z++) for (let y = 0; y < G; y++) for (let x = 0; x < G; x++)
            field[(z * G + y) * G + x] = y - 4.3;
        expectThrow("RangeError", "transvoxel lod -1", () => Mesh.transvoxel(field, G, -1));
        expectThrow("RangeError", "transvoxel lod 4 (stride 16 > 8 cells)", () => Mesh.transvoxel(field, G, 4));
        expectThrow("RangeError", "transvoxel lod 1.5", () => Mesh.transvoxel(field, G, 1.5));
        expectThrow("RangeError", "transvoxel neighbour lod 31",
                    () => Mesh.transvoxel(field, G, 0, [31, -1, -1, -1, -1, -1]));
        if (Mesh.transvoxel(field, G, 3).triangleCount === 0) throw new Error("transvoxel lod 3 should mesh");

        // ── compression: stripify / decode reject what meshopt only asserts ─
        expectThrow("RangeError", "stripify index >= vertexCount", () => Mesh.stripify([0, 1, 5], 3));
        expectThrow("RangeError", "stripify vertexCount -1", () => Mesh.stripify([0, 1, 2], -1));
        const enc = Mesh.box(1, 1, 1).encode();
        const back = Mesh.decode(enc);
        if (back.triangleCount !== 12) throw new Error("decode round trip lost triangles: " + back.triangleCount);
        expectThrow("RangeError", "decode vertexSize 6", () => Mesh.decode(Object.assign({}, enc, { vertexSize: 6 })));
        expectThrow("RangeError", "decode vertexCount -1", () => Mesh.decode(Object.assign({}, enc, { vertexCount: -1 })));
        expectThrow("RangeError", "decode 4e9 vertices",
                    () => Mesh.decode(Object.assign({}, enc, { vertexCount: 4e9 })));

        // ── plants ──────────────────────────────────────────────────────────
        expectThrow("RangeError", "tube sides 2", () => Mesh.tube([[0, 0, 0], [0, 1, 0]], 0.1, 2));
        expectThrow("RangeError", "leafCard widthSegments 0", () => Mesh.leafCard("oval", { widthSegments: 0 }));
        expectThrow("RangeError", "flower petalCount -6", () => Mesh.flower({ petalCount: -6 }));
        expectThrow("RangeError", "blob nsub 20", () => Mesh.blob({ nsub: 20 }));
        const ls = new LSystem("F");
        ls.addRule("F", "FF");
        expectThrow("RangeError", "LSystem.derive -1", () => ls.derive(-1));
        expectThrow("RangeError", "LSystem.derive 1000", () => ls.derive(1000));
        if (ls.derive(3) !== "FFFFFFFF") throw new Error("LSystem.derive(3) gave " + ls.derive(3));

        // ── SDF graph: node ids must exist (a self-reference recursed forever)
        const sdf = Mesh.createSDF();
        const s0 = sdf.sphere(1);
        expectThrow("RangeError", "opUnion with its own id", () => sdf.opUnion(s0, 1));
        expectThrow("RangeError", "translate node -1", () => sdf.translate(-1, [0, 0, 0]));
        expectThrow("RangeError", "marchingCubes dims -4", () => sdf.marchingCubes({ dims: -4 }));
        expectThrow("RangeError", "marchingCubes dims 513", () => sdf.marchingCubes({ dims: [513, 8, 8] }));
        expectThrow("RangeError", "marchingCubes resolution 1e6", () => sdf.marchingCubes({ resolution: 1e6 }));

        // ── PolyMesh: vertex references index the vertex table unchecked ─────
        const pm = PolyMesh.fromPolygon([0, 0, 0, 1, 0, 0, 1, 1, 0], [0, 0, 1]);
        expectThrow("RangeError", "addFace vertex 99", () => pm.addFace([0, 1, 99]));
        expectThrow("RangeError", "fromMeshData index 7", () => PolyMesh.fromMeshData([0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 7]));
        expectThrow("RangeError", "fromPolygons offset past polyVerts",
                    () => PolyMesh.fromPolygons([0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 2], [0, 9]));

        // ── rigging: bone counts, parents, pose/skeleton agreement ──────────
        expectThrow("RangeError", "Pose(-1)", () => new Pose(-1));
        expectThrow("RangeError", "Skeleton parent 5 of 2", () => new Skeleton({ bones: [
            { name: "a", parent: -1 }, { name: "b", parent: 5 } ] }));
        const skel = new Skeleton({ bones: [
            { name: "root", parent: -1, localT: [0, 0, 0] },
            { name: "mid",  parent: 0,  localT: [0, 1, 0] },
            { name: "end",  parent: 1,  localT: [0, 1, 0] } ] });
        expectThrow("RangeError", "computeWorldMatrices with a 1-bone pose", () => new Pose(1).computeWorldMatrices(skel));
        expectThrow("RangeError", "IK.twoBone with a 1-bone pose", () => IK.twoBone(skel, new Pose(1), 0, 1, 2, [1, 1, 0]));
        const pose = skel.bindPose();
        expectThrow("RangeError", "Pose.blend with a short mask", () => Pose.blend(pose, skel.bindPose(), 0.5, [1]));
        expectThrow("RangeError", "AnimationClip values short of times", () => new AnimationClip({ channels: [
            { boneIndex: 0, path: "translation", times: [0, 1], values: [0, 0, 0] } ] }));
        expectThrow("RangeError", "VoxelChunk size -8", () => new VoxelChunk(-8, 8, 8));
        expectThrow("RangeError", "VoxelChunk.set material 300", () => new VoxelChunk(4, 4, 4).set(0, 0, 0, 300));
        const vc = new VoxelChunk(4, 4, 4);
        vc.fill(1);
        expectThrow("RangeError", "buildMesh paletteCount past the palette", () => vc.buildMesh([1, 1, 1, 1], 5));
        if (vc.buildMesh([1, 0, 0, 1, 0, 1, 0, 1], 2).triangleCount === 0) throw new Error("voxel mesh is empty");

        "SUCCESS";
    )JS";

    auto res = bronze::eval::evalScript(script);
    if (res.thrown) {
        std::cerr << "  validation script threw: " << ev::toUtf8(res.value) << std::endl;
        std::exit(1);
    }
    if (ev::toUtf8(res.value) != "SUCCESS") {
        std::cerr << "  validation script returned: " << ev::toUtf8(res.value) << std::endl;
        std::exit(1);
    }
    std::cout << "  validation OK." << std::endl;
}
