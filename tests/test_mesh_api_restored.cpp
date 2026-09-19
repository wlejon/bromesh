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
