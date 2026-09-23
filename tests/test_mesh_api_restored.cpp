// Behavioural coverage for the Mesh / Skeleton members that the QuickJS →
// bronze port dropped or narrowed (bro docs/transition-drift.md rows H7 and
// the mesh rows of build/binding-audit/shape_all.txt).
//
// Linked into test_mesh_api; called from its main(). The checks run real
// geometry rather than probing for names, so a stub would still fail them.

#include "eval/eval.h"
#include "embed/embed.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {
namespace ev = bronze::embed;
} // namespace

// Failures exit the process: assert() is a no-op in the Release configuration
// this test actually runs in, so a thrown script must not be allowed to pass.
void bromeshTestRestoredSurface() {
    std::cout << "[restored] Exercising the members the bronze port dropped..." << std::endl;

    const char* script = R"JS(
        // ── computeTangents(): 4 floats per vertex, finite, non-degenerate ──
        const t = Mesh.box(1, 1, 1);
        t.projectUVs("box", 1);
        const tan = t.computeTangents();
        if (!(tan instanceof Float32Array)) throw new Error("computeTangents: not a Float32Array");
        if (tan.length !== t.vertexCount * 4) throw new Error("computeTangents: expected 4 floats per vertex, got " + tan.length);
        let anyNonZero = false;
        for (let i = 0; i < tan.length; i++) {
            if (!Number.isFinite(tan[i])) throw new Error("computeTangents: non-finite at " + i);
            if (tan[i] !== 0) anyNonZero = true;
        }
        if (!anyNonZero) throw new Error("computeTangents: all-zero tangents");

        // ── simplifyWithAttributes(): returns this, drops triangles ─────────
        const s = Mesh.sphere(1, 32, 32);
        const before = s.triangleCount;
        const ret = s.simplifyWithAttributes(0.3, 0.01, 1.0, 0.5);
        if (ret !== s) throw new Error("simplifyWithAttributes should return this");
        if (s.triangleCount >= before) throw new Error("simplifyWithAttributes did not simplify: " + before + " -> " + s.triangleCount);
        if (s.triangleCount === 0) throw new Error("simplifyWithAttributes collapsed the mesh");

        // ── subdivideMidpoint(): 1 -> 4 triangles per iteration, no smoothing
        const md = Mesh.box(1, 1, 1);
        const mdTris = md.triangleCount;
        const mdBefore = md.bounds();
        const mdRet = md.subdivideMidpoint(1);
        if (mdRet !== md) throw new Error("subdivideMidpoint should return this");
        if (md.triangleCount !== mdTris * 4) throw new Error("subdivideMidpoint: expected " + (mdTris * 4) + " triangles, got " + md.triangleCount);
        // Midpoint subdivision must NOT shrink the surface the way Loop does:
        // the bounding box is unchanged, only the density rises.
        const mdAfter = md.bounds();
        for (let i = 0; i < 3; i++) {
            if (Math.abs(mdAfter.min[i] - mdBefore.min[i]) > 1e-4 ||
                Math.abs(mdAfter.max[i] - mdBefore.max[i]) > 1e-4) {
                throw new Error("subdivideMidpoint moved the surface on axis " + i);
            }
        }

        // ── hasSelfIntersections / findSelfIntersections ────────────────────
        const clean = Mesh.box(1, 1, 1);
        if (clean.hasSelfIntersections() !== false) throw new Error("a box should not self-intersect");
        const pairs = clean.findSelfIntersections();
        if (!Array.isArray(pairs)) throw new Error("findSelfIntersections should return an array");
        if (pairs.length !== 0) throw new Error("a box reported " + pairs.length + " self-intersections");

        // A box merged with a copy of itself shifted a quarter width overlaps
        // itself, so the pair list must be non-empty and carry triA/triB.
        const a0 = Mesh.box(1, 1, 1);
        const b0 = Mesh.box(1, 1, 1);
        b0.translate(0.25, 0, 0);
        const overlap = Mesh.merge([a0, b0]);
        if (!overlap.hasSelfIntersections()) throw new Error("overlapping boxes should self-intersect");
        const op = overlap.findSelfIntersections();
        if (op.length === 0) throw new Error("findSelfIntersections found no pairs for overlapping boxes");
        if (typeof op[0].triA !== "number" || typeof op[0].triB !== "number") throw new Error("findSelfIntersections entries need triA/triB");

        // ── intersectsMesh(other) ───────────────────────────────────────────
        const c0 = Mesh.box(1, 1, 1);
        const c1 = Mesh.box(1, 1, 1);
        c1.translate(0.25, 0, 0);
        if (c0.intersectsMesh(c1) !== true) throw new Error("overlapping boxes should intersect");
        const c2 = Mesh.box(1, 1, 1);
        c2.translate(10, 0, 0);
        if (c0.intersectsMesh(c2) !== false) throw new Error("distant boxes should not intersect");
        let threw = false;
        try { c0.intersectsMesh(42); } catch (e) { threw = true; }
        if (!threw) throw new Error("intersectsMesh should reject a non-Mesh argument");

        // ── convexDecomposition(optionsObject) ──────────────────────────────
        // The object form is what the pre-transition binding read; the
        // positional form has to keep working alongside it.
        const cd = Mesh.sphere(1, 16, 16);
        const hullsObj = cd.convexDecomposition({ maxHulls: 4, maxVerticesPerHull: 16, resolution: 10000, minVolumePerHull: 0.001 });
        if (!Array.isArray(hullsObj)) throw new Error("convexDecomposition should return an array");
        if (hullsObj.length === 0) throw new Error("convexDecomposition returned no hulls");
        if (hullsObj.length > 4) throw new Error("convexDecomposition ignored maxHulls: " + hullsObj.length);
        const hullsPos = cd.convexDecomposition(4, 16, 10000, 0.001);
        if (!Array.isArray(hullsPos) || hullsPos.length === 0) throw new Error("positional convexDecomposition regressed");

        // ── buildMeshlets(): bounds carry the cull cone, not just a radius ──
        const ml = Mesh.sphere(1, 24, 24);
        const mls = ml.buildMeshlets({ maxVertices: 64, maxTriangles: 124, coneWeight: 0.5 });
        if (!Array.isArray(mls) || mls.length === 0) throw new Error("buildMeshlets returned nothing");
        const bnd = mls[0].bounds;
        for (const key of ["center", "coneApex", "coneAxis"]) {
            if (!Array.isArray(bnd[key])) throw new Error("buildMeshlets bounds." + key + " missing");
            if (bnd[key].length !== 3) throw new Error("buildMeshlets bounds." + key + " should have 3 components");
            for (const v of bnd[key]) if (!Number.isFinite(v)) throw new Error("bounds." + key + " non-finite");
        }
        if (typeof bnd.radius !== "number" || typeof bnd.coneCutoff !== "number") throw new Error("buildMeshlets bounds lost radius/coneCutoff");
        // A unit sphere's meshlet centres live on/near the surface, so the
        // cone axis must be a real direction rather than the zero vector.
        const axLen = Math.hypot(bnd.coneAxis[0], bnd.coneAxis[1], bnd.coneAxis[2]);
        if (axLen < 1e-3) throw new Error("buildMeshlets coneAxis is degenerate");

        // ── Skeleton.addRigifySockets() ─────────────────────────────────────
        const rig = new Skeleton({
            bones: [
                { name: "DEF-spine",        parent: -1, localT: [0, 0, 0] },
                { name: "DEF-hand.L",       parent: 0,  localT: [1, 0, 0] },
                { name: "DEF-hand.R",       parent: 0,  localT: [-1, 0, 0] },
                { name: "DEF-head",         parent: 0,  localT: [0, 1, 0] }
            ]
        });
        const added = rig.addRigifySockets();
        if (typeof added !== "number") throw new Error("addRigifySockets should return a count");
        if (added <= 0) throw new Error("addRigifySockets added nothing to a Rigify-named skeleton");

        // A skeleton with no recognisable bone names gets nothing, and the
        // call still answers with 0 rather than throwing.
        const plain = new Skeleton({ bones: [{ name: "zzz", parent: -1, localT: [0, 0, 0] }] });
        if (plain.addRigifySockets() !== 0) throw new Error("addRigifySockets matched an unrelated skeleton");

        "SUCCESS";
    )JS";

    auto res = bronze::eval::evalScript(script);
    if (res.thrown) {
        std::cerr << "  restored-surface script threw: " << ev::toUtf8(res.value) << std::endl;
        std::exit(1);
    }
    if (ev::toUtf8(res.value) != "SUCCESS") {
        std::cerr << "  restored-surface script returned: " << ev::toUtf8(res.value) << std::endl;
        std::exit(1);
    }
    std::cout << "  restored Mesh/Skeleton members OK." << std::endl;
}

// Bindings fixed during the GC-rooting audit. Every path here reads several
// properties or calls back into JS, so under BRONZE_GC_STRESS=1 a raw Value
// held across one of those allocations would crash or misread.
void bromeshTestBindingFixes() {
    std::cout << "[fixes] Exercising the bindings fixed in the GC audit..." << std::endl;

    const char* script = R"JS(
        // ── SkinData weights/indices aliases return the arrays ─────────────
        const sd = new SkinData({ weights: [1, 0, 0, 0, 0.5, 0.5, 0, 0], indices: [0, 0, 0, 0, 0, 1, 0, 0], boneCount: 2 });
        if (!(sd.weights instanceof Float32Array)) throw new Error("SkinData.weights is not a Float32Array: " + sd.weights);
        if (sd.weights.length !== 8 || sd.weights[4] !== 0.5) throw new Error("SkinData.weights has the wrong contents");
        if (!(sd.indices instanceof Uint32Array)) throw new Error("SkinData.indices is not a Uint32Array: " + sd.indices);
        if (sd.indices[5] !== 1) throw new Error("SkinData.indices has the wrong contents");

        // ── alias methods forward to their targets and propagate throws ────
        const box = Mesh.box(1, 1, 1);  // half-extents: a 2x2x2 cube
        const bb = box.computeBBox();
        if (!bb || Math.abs(bb.max[0] - bb.min[0] - 2) > 1e-4) throw new Error("computeBBox did not forward to bounds()");
        if (Math.abs(box.computeVolume() - 8) > 1e-3) throw new Error("computeVolume did not forward to volume()");
        let threw = false;
        try { box.union(42); } catch (e) { threw = true; }
        if (!threw) throw new Error("union(non-mesh) should throw, not return the error");

        // ── transvoxel honours neighborLods ────────────────────────────────
        // The plane y = 4.3 crosses the +X face. A coarser +X neighbour
        // (LOD 1, 2-cell grid) snaps that face's vertices onto its grid, so
        // they drop to y = 4 while the no-neighbour chunk stays flat at 4.3.
        const G = 9;
        const field = new Float32Array(G * G * G);
        for (let z = 0; z < G; z++) for (let y = 0; y < G; y++) for (let x = 0; x < G; x++)
            field[(z * G + y) * G + x] = y - 4.3;
        const plain = Mesh.transvoxel(field, G, 0, [-1, -1, -1, -1, -1, -1], 0, 1);
        const seam = Mesh.transvoxel(field, G, 0, [1, -1, -1, -1, -1, -1], 0, 1);
        if (plain.triangleCount === 0) throw new Error("transvoxel produced no surface");
        if (Math.abs(plain.bounds().min[1] - 4.3) > 1e-3) throw new Error("transvoxel plane is not at y=4.3");
        if (Math.abs(seam.bounds().min[1] - 4.0) > 1e-3)
            throw new Error("transvoxel ignored neighborLods (minY " + seam.bounds().min[1] + ")");
        threw = false;
        try { Mesh.transvoxel(new Float32Array(8), G); } catch (e) { threw = true; }
        if (!threw) throw new Error("transvoxel should reject a field smaller than gridSize^3");

        // ── IK option-object wrappers reach the positional solvers ─────────
        const skel = new Skeleton({ bones: [
            { name: "root", parent: -1, localT: [0, 0, 0] },
            { name: "mid",  parent: 0,  localT: [0, 1, 0] },
            { name: "end",  parent: 1,  localT: [0, 1, 0] }
        ]});
        const pose = skel.bindPose();
        const r1 = IK.solveTwoBone({ skel, pose, root: 0, mid: 1, end: 2, targetPos: [1, 1, 0], poleVector: [0, 0, 1] });
        if (typeof r1 !== "boolean") throw new Error("solveTwoBone should return a boolean");
        const r2 = IK.solveFabrik({ skel, pose, chain: [0, 1, 2], targetPos: [0.5, 1.5, 0], maxIterations: 20, tolerance: 1e-4 });
        if (typeof r2 !== "boolean") throw new Error("solveFabrik should return a boolean");
        const r3 = IK.solveLookAt({ skel, pose, bone: 0, targetPos: [1, 0, 0], forward: [0, 1, 0], up: [0, 0, 1] });
        if (typeof r3 !== "boolean") throw new Error("solveLookAt should return a boolean");
        threw = false;
        try { IK.solveTwoBone({ skel, pose, root: 0, mid: 1, end: 2 }); } catch (e) { threw = true; }
        if (!threw) throw new Error("solveTwoBone without targetPos should throw");

        // ── applyMorphTarget object form returns this and moves vertices ───
        const mm = Mesh.box(1, 1, 1);
        const n = mm.vertexCount;
        const dp = new Float32Array(n * 3);
        for (let i = 0; i < n; i++) dp[i * 3 + 1] = 1;
        const y0 = mm.bounds().max[1];
        if (mm.applyMorphTarget({ name: "up", deltaPositions: dp, weight: 0.5 }) !== mm)
            throw new Error("applyMorphTarget should return this");
        if (Math.abs(mm.bounds().max[1] - (y0 + 0.5)) > 1e-4) throw new Error("applyMorphTarget weight not applied");

        "SUCCESS";
    )JS";

    auto res = bronze::eval::evalScript(script);
    if (res.thrown) {
        std::cerr << "  binding-fixes script threw: " << ev::toUtf8(res.value) << std::endl;
        std::exit(1);
    }
    if (ev::toUtf8(res.value) != "SUCCESS") {
        std::cerr << "  binding-fixes script returned: " << ev::toUtf8(res.value) << std::endl;
        std::exit(1);
    }
    std::cout << "  binding fixes OK." << std::endl;
}
