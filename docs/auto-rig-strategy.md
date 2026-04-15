# bromesh auto-rigging: strategy

*Status: strategy document. Not an implementation plan. A fresh session should read this, pick the section it judges most appropriate to tackle first, then produce its own implementation plan for that section.*

---

## 1. Context

bromesh already has the entire **runtime** side of skeletal animation: glTF skin/animation load, pose evaluation with all three interpolation modes, pose blending with bone masks, skinning matrix composition, two-bone IK, FABRIK, look-at, name-based animation retargeting, Rigify-style socket discovery, weight transfer between meshes, shrinkwrap, cross-mesh normal/AO bake.

The remaining gap is **authoring**: you can't yet go from `MeshData` → `(Skeleton, SkinData)` inside bromesh. Today, that step happens in Blender (Rigify) or via external services (Mixamo, AccuRIG). Everything after it is in-engine.

Closing this gap removes the last external dependency in the text-to-animated-character pipeline. That matters for three reasons:

1. **Pipeline integration.** A bromesh auto-rig runs in-process. No containerized Blender, no Mixamo web round-trip, no licensing, no account, no network. Part of the game build, not adjacent to it.
2. **Scale.** Producing hundreds of character variants procedurally with a per-character Blender/Mixamo round-trip is either slow (network) or fragile (headless Blender). In-process rigging is milliseconds per character, trivially parallelizable.
3. **Generality.** Rigify is humanoid-first with quadruped support bolted on. Mixamo is humanoid-only. If we want spiders, centipedes, tentacle monsters, or custom skeleton topologies, the existing tools fight us. A generic data-driven auto-rigger doesn't.

This is not about competing with paid software on features. It's about owning the whole pipeline.

---

## 2. Goals

**Hard requirements for the finished system:**

- Accept a mesh + a declarative **rig spec** (skeleton template + landmark schema) + optional user-marked landmarks, and produce `(Skeleton, SkinData)` ready to drive with the existing `animation/` module.
- Support arbitrary skeleton topologies via the rig spec, including at least: biped (2 legs + arms), quadruped (4 legs), hexapod (6 legs), octopod (8 legs). Leg count is just a parameter of the spec, not a special case.
- Produce weights good enough that skinning deformations look right on typical character meshes without manual weight painting. "Right" = visually indistinguishable from Blender's default auto-weights at normal viewing distance; no catastrophic collapses at joint bends.
- Handle messy topology (non-manifold, holes, floating components, MeshyAI-style blob output) gracefully. Fall back from high-quality algorithms to more robust ones when the mesh can't support them.
- Emit sockets (named attachment points with offsets) as part of the output, derived from the rig spec. The existing `addRigifySockets` helper becomes one instance of a more general facility.
- Pair with a **procedural locomotion** module that generates walk/run cycles for any rigged character, adapting to leg count, stride length, and gait pattern. This is what makes "arbitrary creatures" viable without authoring animations per species.

**Quality bar — what SOTA means here:**

- Weighting: geodesic voxel binding (Dionne & de Lasa 2013, Maya's default) as the primary algorithm. Bounded biharmonic weights (Jacobson et al. 2011) as an optional higher-quality path if worth the QP solver. Surface bone heat (Baran & Popović 2007, Blender's default) as a cheap fallback for clean topology.
- Skeleton fitting: landmark-driven warping of a rig template, where landmarks come from either user clicks or geometric heuristics (extrema detection, symmetry axes, principal component analysis). A plug point for an ML landmark detector when/if one is trained.
- Locomotion: footfall-pattern-driven cycle generation with IK foot placement on uneven ground. Hand-authored cycles still supported via glTF import; procedural fills gaps.

**Non-negotiable:**

The output is indistinguishable from a hand-rigged character to the rest of the engine. It's a standard `Skeleton` + `SkinData` pair. No special-case runtime code. This means iteration and replacement of algorithms is safe — they just produce better-quality versions of the same data.

---

## 3. Non-goals

- **Animation authoring UI.** We don't build Blender. Animations are authored externally (Blender, Mixamo, motion capture) or generated procedurally (the locomotion module). We don't support keyframing by hand in an editor.
- **Control rigs.** Rigify's IK/FK switches, pole targets, custom bone shapes, and Python drivers are *authoring-time* ergonomics for animators dragging handles. They don't exist at runtime — Blender bakes them out on export. We don't need them. Runtime IK is already handled by our solvers.
- **Physically-based simulation.** Cloth, hair, jiggle physics, ragdoll — these are adjacent domains handled by Jolt or a dedicated solver, not the rigger.
- **ML-based full automation on day one.** A trained skeleton-placement network (RigNet-style) is ideal but out of scope for the initial system. The architecture must leave a clean plug-in point for it, but ship with landmark-driven placement first.
- **Face rigging.** Facial blendshapes and jaw/eye bones are a distinct problem from body rigging and deserve their own treatment. The body rigger should produce a head bone; face work lives in a separate module.

---

## 4. Architecture overview

Four conceptual sub-systems. Each is independently replaceable.

### 4.1 Rig spec (data model)

A declarative description of a skeleton template. Not a humanoid hard-coded structure; a generic schema that can describe any creature.

Contents of a rig spec:
- **Bone list.** Names, parent indices, relative rest positions (as fractions of a reference bounding box, or as offsets from named landmarks).
- **Landmark schema.** Named points on the mesh the template aligns to (e.g. "pelvis", "wrist_L", "ankle_R", "crown"). Each bone's rest position is expressed in terms of landmarks: "midpoint(pelvis, chest)", "ankle_R + 0.3 * (ankle_R - knee_R)", etc.
- **Chain annotations.** Which bones form IK chains (arms, legs, spines, tails). Used later by the locomotion module and by the solver to pick appropriate IK algorithms.
- **Socket declarations.** Named attachment points with bone + offset.
- **Symmetry hints.** Left/right pairs, so symmetry can be enforced during fitting.

Stored on disk as a simple text format (JSON or TOML). Rig specs for humanoid, quadruped, hexapod, and octopod ship with bromesh. Users can author custom ones.

This is the single most important design decision. Get the schema right and everything downstream is straightforward.

### 4.2 Landmark detection

Three modes, used interchangeably:

1. **User-marked.** The user clicks landmarks on the mesh in a host application. Fastest to ship, most reliable, ~30s of user interaction per character. Ship this first.
2. **Heuristic.** Geometric analysis finds landmarks from shape alone: principal component analysis for main axis, extrema finding for head/hands/feet, symmetry-plane detection, medial axis skeleton for rough bone layout. Works for T-pose humanoids; degrades on non-T-pose or complex creatures.
3. **Learned.** A trained network emits landmark positions from mesh geometry. Best accuracy, requires a model. Ship the plug point; train or adopt a model later.

The interface is the same for all three: given a mesh and a rig spec's landmark schema, return `{landmarkName: position}`. Downstream code doesn't care which mode produced them.

### 4.3 Skeleton fitter

Takes a rig spec + resolved landmark positions + mesh, produces a `Skeleton` whose bones sit inside the mesh at anatomically plausible positions.

Steps:
- Resolve each bone's rest position from its landmark expression.
- Project bones inward from the surface (so they sit inside the mesh volume, not on the skin).
- Apply symmetry constraints if the spec declares them.
- Compute inverse bind matrices from final rest transforms.
- Validate: no zero-length bones, no disconnected chains, all expected bones present.

Not algorithmically hard once the rig spec + landmarks are clean. Bulk of the work is in the spec.

### 4.4 Weighting

Three algorithm choices, ranked by quality:

1. **Bounded biharmonic weights (BBW).** Solves a constrained QP minimizing a biharmonic energy. Smoothly varying, respects mesh interior, provably non-negative, sums to 1. Highest quality, especially for overlapping bone influence regions. Requires a QP solver (not currently in bromesh). Maybe 800 LOC including solver, or a third-party dependency.
2. **Geodesic voxel binding.** Voxelize the mesh, compute per-voxel distance fields from each bone's line segment using geodesic distance in the voxel grid (not Euclidean — respects mesh interior cavities). Sample back onto vertices, apply falloff, normalize, keep top 4 bones. Robust to non-manifold meshes because it operates on voxels, not surface topology. Quality within a few percent of BBW for typical characters. **This is the primary choice.** Reuses bromesh's existing voxelization.
3. **Surface bone heat.** Solve heat equation on mesh surface with each bone as a heat source. Blender's default. Requires manifold mesh. Fastest of the three. Fallback for clean topology when voxel binding is overkill.

The system picks automatically: check manifoldness → bone heat if clean, voxel binding otherwise. BBW is opt-in for users who want maximum quality and have a clean mesh.

### 4.5 Procedural locomotion (companion system)

Not strictly part of rigging, but the rigger's output is useless for non-humanoid creatures without a way to animate them. Motion capture doesn't exist for spiders.

Generates walk / trot / run cycles from:
- A gait pattern (which feet are lifted at which phase — a per-leg phase offset in [0, 1]).
- Stride length and cycle duration.
- The rig spec's leg chain annotations.

Algorithm:
- Each cycle, each foot follows a trajectory: planted phase (stationary on ground), lift phase (arc in the direction of travel).
- Foot world positions drive two-bone IK (or FABRIK for longer chains) via bromesh's existing solvers.
- Spine/neck get a small counter-sway synthesized from hip motion.
- Head optionally drives look-at toward a target.

Ships with default gait patterns for biped, quadruped, hexapod, octopod. Users can define custom gaits.

---

## 5. Design principles

**Data-driven, not hard-coded.** No function anywhere says "if humanoid then...". Everything flows from the rig spec. Adding a new creature type = authoring a rig spec, not writing code.

**Graceful degradation.** Each algorithm has a cheaper fallback. Users shouldn't hit a cliff where bromesh refuses to rig a mesh because a high-end algorithm's preconditions aren't met.

**Same output shape as hand-rigged.** `Skeleton` + `SkinData` + optional `Animation`. No parallel code paths downstream.

**Parallelizable at the character level.** Rigging 1000 characters should scale linearly across cores. No global state, no shared caches that need invalidation.

**Testable deterministic output.** Given the same inputs, produce bit-identical output. No randomness without explicit seeds. This is the prerequisite for regression testing weight quality as we swap algorithms.

---

## 6. Phased roadmap

Phases are **dependency-ordered**, not chronological. Earlier phases unblock later ones. Each phase is shippable on its own — bromesh is useful at every step.

### Phase 0: rig spec schema and format

Define the data model for rig specs. Write the humanoid spec. No code yet uses it, but other phases target this schema.

### Phase 1: user-marked landmarks + skeleton fitter + voxel weighting

The first end-to-end slice. User marks landmarks, fitter places a humanoid skeleton, voxel-binding computes weights. Output goes through existing glTF save + existing runtime and plays the existing animation library. This is the "replace Blender for humanoids" deliverable.

### Phase 2: additional rig specs

Author quadruped, hexapod, octopod rig specs. Verify Phase 1's code needs zero changes to handle them — if it does, the rig spec schema needs another iteration.

### Phase 3: procedural locomotion

Gait-driven cycle generator. Wired through existing IK solvers. Ships default gaits for 2/4/6/8 legs.

### Phase 4: quality upgrades

Bounded biharmonic weights as an opt-in. Manifold-safe surface bone heat fallback. Weight post-processing (smoothing, outlier rejection).

### Phase 5: heuristic landmark detection

Geometric auto-landmark for humanoids + quadrupeds. Removes the user-click step for standard cases while keeping it available for custom creatures.

### Phase 6: ML landmark detection plug-in

Define the interface for an external landmark detector. Ship a reference ONNX-based implementation if a suitable model is found or trained. Must be optional — bromesh still builds without it.

### Phase 7: face rig (separate scope, named here for completeness)

Blendshape + jaw/eye bone rigging. Logically part of character authoring but substantially different techniques. Plan separately.

---

## 7. Key algorithm trade-offs

### 7.1 Voxel-grid resolution for binding

The voxel resolution trades quality for speed and memory. 64³ is fast, OK quality for simple shapes. 128³ is the sweet spot for humanoids. 256³ captures fine detail (fingers, face) at 8x the cost.

Strategy: adaptive — resolution scales with bounding box, tuned to yield ~5-10 voxels per bone cross-section. Fingers need fine voxels; torsos don't.

### 7.2 Surface vs. voxel heat

Surface heat needs a manifold mesh; fails catastrophically on non-manifold input. Voxel heat shrugs off topology problems. Given MeshyAI blobs are frequently non-manifold, voxel is the right default. Surface is the optimization.

### 7.3 Top-N weights per vertex

glTF standard is 4 weights. bromesh's runtime uses 4. Dual-quaternion skinning and some high-quality pipelines use 8. Stay at 4 for compatibility; the quality ceiling of 4-weight skinning isn't where we're bottlenecked.

### 7.4 Authored animations vs. procedural

Authored (mocap, keyframed) are richer and higher quality for humanoids. Procedural is the only option for creatures with no existing animation library. Both must coexist; the rig spec + skeleton shape is the same either way.

### 7.5 Symmetry enforcement

Symmetric characters (L/R) benefit from enforced symmetry during skeleton fitting and weighting — it eliminates small asymmetries that cause visual tells. But some characters are asymmetric by design (one-armed, prosthetics, alien). Symmetry is a per-rig-spec opt-in, not global.

---

## 8. What reuses what already exists in bromesh

Significant infrastructure is already in place. This is not a from-scratch build.

- **Voxelization** for geodesic voxel weighting: `VoxelChunk`, marching cubes / surface nets.
- **BVH + closest-point + raycast** for landmark projection, symmetry detection, weight post-processing.
- **Skeleton/SkinData/Animation** data types and runtime: entirely done.
- **IK solvers** (two-bone, FABRIK, look-at) for procedural locomotion.
- **Pose evaluation + blending** for mixing authored and procedural motion.
- **glTF save with skin + animations** for exporting finished rigged characters.
- **Socket metadata** on `Skeleton` for equipment attachment.

The new code is concentrated in: rig spec parser, skeleton fitter, landmark detection, voxel weighting, locomotion. Everything else hangs off what's already built.

---

## 9. Success criteria

For the system as a whole:

1. A MeshyAI humanoid blob → fully rigged character with playable walk cycle, no Blender, no human intervention beyond marking landmarks once (Phase 1). Measured: visually acceptable deformation at shoulders, elbows, knees, hips during a standard walk-cycle test.
2. Same for a quadruped (wolf, horse) (Phase 2+3). No code changes from the humanoid path — only a different rig spec.
3. Same for a hexapod and octopod (Phase 2+3). Spiders walk.
4. Produced rigs round-trip through glTF export and reload into bromesh or any external glTF viewer without loss.
5. Rigging a character takes under 2 seconds wall-clock for a 20k-triangle mesh on a modern desktop.
6. Rigging is deterministic: same input produces byte-identical output across runs.

For individual phases, the success bar is "does it pass tests covering the quality and correctness claims above, using meshes representative of the target use cases."

---

## 10. Open questions a fresh session should resolve

These are decisions I don't want to prescribe:

1. **Rig spec format: JSON, TOML, or custom.** Trade-off is tool ergonomics vs. dependencies. bromesh currently has no parser for any of these; adding one is a minor dependency choice.
2. **Landmark marking UI: part of bromesh, or host-app responsibility.** bromesh is a headless library. A landmark-marking tool needs 3D viewport interaction, which is the host's domain (bro's UI, Blender add-on, command-line with numeric input, etc.). The library's job is probably to accept a dict of `{name: position}` and nothing else. But worth deciding whether bromesh ships a reference landmark-marking tool or leaves it entirely to callers.
3. **QP solver for BBW: dependency or from-scratch.** BBW needs a constrained QP solve. Options: OSQP (small, permissively licensed, a few thousand LOC), Eigen's own QP routines, or implement primal active-set ourselves. Not a Phase-1 decision but worth thinking about before Phase 4.
4. **Voxel weighting: single-thread first or parallelized from day one.** At 128³ resolution with 30+ bones it's not small work. Worth deciding early since it affects code structure.
5. **How much of the rig spec is Turing-complete expression evaluation.** "midpoint(pelvis, chest)" is a tiny expression language. The scope of that language determines how much authoring flexibility rig specs have. Too small = every creature needs special-case landmarks. Too big = it's a programming language inside a config file. Find the minimum expressive set.
6. **Where to ship reference rig specs.** Bromesh repo, separate repo, or bundled as headers/strings. Affects distribution but not architecture.

A fresh session should pick Phase 0 or Phase 1 (most likely Phase 0 — the rig spec schema is the architectural keystone) and produce a concrete implementation plan for it. Phase 0 is mostly design work and data format selection; Phase 1 is the first substantial code. Whichever phase is picked first, earlier phases in the roadmap have no code dependencies to worry about, so the choice is mostly about what's more useful to have settled first.

---

## 11. Relationship to the broader character pipeline

For context, the broader pipeline this lives inside:

```
text/image → MeshyAI (external, generation)
           → bromesh (cleanup, UV, detail bake)                    [done]
           → MeshyAI (external, texturing)
           → bromesh auto-rig                                       [this document]
           → bromesh animation (IK, blending, retargeting, sockets) [done]
           → bromesh runtime (weight transfer for equipment, etc.)  [done]
           → glTF export for external use                           [done]
```

MeshyAI (generation + texturing) is the only external dependency we don't plan to replace — it's a vision/generative problem rather than a math one. Everything else is bromesh-native once the auto-rigger lands.
