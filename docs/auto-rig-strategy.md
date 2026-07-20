# bromesh auto-rigging: strategy

*Status: strategy document — the design rationale behind the shipped `rigging/`
module, kept as the north star for changes to it. Section 6 tracks what is built
and what is not; section 10 records which design questions are settled. Sections
1-5 and 7 describe the reasoning and remain the reference for why the subsystem
is shaped the way it is.*

---

## 1. Context

bromesh already has the entire **runtime** side of skeletal animation: glTF skin/animation load, pose evaluation with all three interpolation modes, pose blending with bone masks, skinning matrix composition, two-bone IK, FABRIK, look-at, name-based animation retargeting, Rigify-style socket discovery, weight transfer between meshes, shrinkwrap, cross-mesh normal/AO bake.

The gap this document set out to close was **authoring**: going from `MeshData` → `(Skeleton, SkinData)` inside bromesh instead of in Blender (Rigify) or via external services (Mixamo, AccuRIG). That path now exists as the `rigging/` module.

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

## 6. Roadmap and current state

The roadmap was dependency-ordered, not chronological, with each step shippable
on its own. Everything through heuristic landmark detection is now implemented.

**Built:**

- **Rig spec schema and format** — `rigging/rig_spec.h`. JSON in and out
  (`parseRigSpecJSON` / `serializeRigSpecJSON` / `loadRigSpecFile`), bone
  positions as landmark expressions.
- **Landmarks + skeleton fitter + voxel weighting** — the end-to-end slice.
  `rigging/landmarks.h`, `rigging/skeleton_fit.h`, `rigging/voxel_bind.h`,
  driven by `autoRig()` in `rigging/auto_rig.h`; output exports through the
  existing glTF save and plays through the existing animation runtime.
- **Additional rig specs** — quadruped, hexapod, octopod, alongside the
  humanoid. Bundled as builtins and as data under `data/rig_specs/`; they run
  the humanoid code path unchanged.
- **Procedural locomotion** — `animation/locomotion.h`, gait-driven cycle
  generation on top of the IK solvers, with default gaits for 2/4/6/8 legs.
- **Quality upgrades** — bounded biharmonic weights as an opt-in
  (`rigging/bbw.h`, OSQP-backed), surface bone heat as the manifold default
  (`rigging/bone_heat.h`), Laplacian weight post-processing
  (`rigging/weight_smooth.h`) and `rigging/skin_validate.h`.
- **Heuristic landmark detection** — `rigging/landmark_detect.h`, geometric
  auto-landmarking that removes the user-click step for standard cases while
  leaving hand-authored landmarks available.

**Not built:**

- **ML landmark detection plug-in.** An interface for an external detector plus
  a reference ONNX-based implementation. Must stay optional — bromesh builds
  without it.
- **Face rig.** Blendshape + jaw/eye bone rigging. Logically part of character
  authoring but substantially different techniques; a separate scope.

---

## 7. Key algorithm trade-offs

### 7.1 Voxel-grid resolution for binding

The voxel resolution trades quality for speed and memory. 64³ is fast, OK quality for simple shapes. 128³ is the sweet spot for humanoids. 256³ captures fine detail (fingers, face) at 8x the cost.

Strategy: adaptive — resolution scales with bounding box, tuned to yield ~5-10 voxels per bone cross-section. Fingers need fine voxels; torsos don't.

### 7.2 Surface vs. voxel heat

Surface heat needs a manifold mesh; fails catastrophically on non-manifold input. Voxel heat shrugs off topology problems. Neither is the default: `WeightingMethod::Auto` inspects manifoldness per mesh (`autoSelectWeightingMethod`) and picks surface bone heat when the mesh supports it, voxel-bind otherwise. That keeps the quality of surface heat for clean input without a cliff on MeshyAI-style blobs.

### 7.3 Top-N weights per vertex

glTF standard is 4 weights. bromesh's runtime uses 4. Dual-quaternion skinning and some high-quality pipelines use 8. Stay at 4 for compatibility; the quality ceiling of 4-weight skinning isn't where we're bottlenecked.

### 7.4 Authored animations vs. procedural

Authored (mocap, keyframed) are richer and higher quality for humanoids. Procedural is the only option for creatures with no existing animation library. Both must coexist; the rig spec + skeleton shape is the same either way.

### 7.5 Symmetry enforcement

Symmetric characters (L/R) benefit from enforced symmetry during skeleton fitting and weighting — it eliminates small asymmetries that cause visual tells. But some characters are asymmetric by design (one-armed, prosthetics, alien). Symmetry is a per-rig-spec opt-in, not global.

---

## 8. What the rigging module reuses

The rigging module is a thin layer over infrastructure that already existed:

- **Voxelization** for geodesic voxel weighting: `VoxelChunk`, marching cubes / surface nets.
- **BVH + closest-point + raycast** for landmark projection, symmetry detection, weight post-processing.
- **Skeleton/SkinData/Animation** data types and runtime: entirely done.
- **IK solvers** (two-bone, FABRIK, look-at) for procedural locomotion.
- **Pose evaluation + blending** for mixing authored and procedural motion.
- **glTF save with skin + animations** for exporting finished rigged characters.
- **Socket metadata** on `Skeleton` for equipment attachment.

The new code is concentrated in: rig spec parser, skeleton fitter, landmark detection, weighting (voxel-bind, bone heat, BBW), weight smoothing, locomotion. Everything else hangs off what was already built.

---

## 9. Success criteria

For the system as a whole:

1. A MeshyAI humanoid blob → fully rigged character with playable walk cycle, no Blender, no human intervention beyond marking landmarks once. Measured: visually acceptable deformation at shoulders, elbows, knees, hips during a standard walk-cycle test.
2. Same for a quadruped (wolf, horse). No code changes from the humanoid path — only a different rig spec.
3. Same for a hexapod and octopod. Spiders walk.
4. Produced rigs round-trip through glTF export and reload into bromesh or any external glTF viewer without loss.
5. Rigging a character takes under 2 seconds wall-clock for a 20k-triangle mesh on a modern desktop.
6. Rigging is deterministic: same input produces byte-identical output across runs.

For any individual change, the success bar is "does it pass tests covering the quality and correctness claims above, using meshes representative of the target use cases."

---

## 10. Resolved decisions and remaining questions

Settled:

1. **Rig spec format: JSON**, with a hand-rolled parser in `rig_spec.cpp` — no
   parser dependency added.
2. **QP solver for BBW: OSQP**, as an optional submodule. BBW compiles to a
   no-op returning failure when the submodule is absent.
3. **Rig spec expression language: a fixed, closed set** — `landmark:A`,
   `mid:A,B`, `lerp:A,B,t`, `offset:A,dx,dy,dz`. No control flow, no
   user-defined functions.
4. **Reference rig specs ship in this repo**, both as compiled builtins
   (`builtinHumanoidSpec()` and friends) and as JSON under `data/rig_specs/`,
   so callers can start from data without a rebuild.
5. **Landmark marking UI is the host's responsibility.** The library accepts a
   `{name: position}` dict and nothing else; `landmark_detect.h` covers the
   standard cases without any UI.

Still open:

1. **Parallelizing weight computation.** All weighting methods are currently
   single-threaded. At 128³ voxel resolution with 30+ bones this is the bulk of
   auto-rig wall-clock, and it is embarrassingly parallel per bone.
2. **ML landmark detection.** The plug-in interface is undesigned; a suitable
   model has not been selected or trained.

---

## 11. Relationship to the broader character pipeline

For context, the broader pipeline this lives inside:

```
text/image → MeshyAI (external, generation)
           → bromesh (cleanup, UV, detail bake)                    [done]
           → MeshyAI (external, texturing)
           → bromesh auto-rig                                       [done]
           → bromesh animation (IK, blending, retargeting, sockets) [done]
           → bromesh runtime (weight transfer for equipment, etc.)  [done]
           → glTF export for external use                           [done]
```

MeshyAI (generation + texturing) is the only external dependency we don't plan to replace — it's a vision/generative problem rather than a math one. Everything else is bromesh-native once the auto-rigger lands.
