#pragma once

// Leaf clumps: foliage as distinct masses of leaves rather than a uniform
// spray along the twigs. A clump is a lumpy ball (a few overlapping lobes)
// hung at a twig end; its leaves sit in a shell just under the lobes' surface
// and face out of it, so each clump lights as one rounded form and the sky
// shows between clumps.
//
// Placement is a pure function of (clump, its index, the leaf's global index,
// the options): bro.scene's GPU scatter (scene_procedural.glsl, clump mode)
// runs the same hash and the same steps, so a scatter node draws leaf for
// leaf what `placeLeavesInClumps` returns.

#include <cstdint>
#include <vector>

#include "bromesh/procedural/leaf_scatter.h"

#include <bromath/vec.h>

namespace bromesh {

/// One clump of leaf mass.
struct LeafClump {
    bromath::Vec3 center{};
    float radius = 0.5f;
    /// The direction of the twig the clump grows from (unit).
    bromath::Vec3 dir{0.0f, 1.0f, 0.0f};
    /// Multiplier on `LeafClumpOptions::perArea` (shaded clumps carry fewer leaves).
    float density = 1.0f;
    /// The caller's index of the branch segment the clump hangs from; carried
    /// through untouched (broflora maps it to the segment's wind sway).
    int32_t segment = -1;
};

struct LeafClumpOptions {
    /// Leaves per unit^2 of clump surface (4 pi r^2) at density 1.
    float perArea = 40.0f;
    /// 0: a leaf faces straight out of its lobe; 1: it faces up.
    float upBias = 0.3f;
    /// How far a leaf's normal strays from that, as a fraction of it.
    float tiltJitter = 0.6f;
    /// Depth of the shell the leaves sit in, as a fraction of the lobe radius
    /// (0: all on the surface).
    float shell = 0.4f;
    float baseScale = 1.0f;
    /// +/- fraction of `baseScale` per leaf.
    float scaleJitter = 0.25f;
    uint32_t seed = 0;
};

/// Lobes per clump: lobe 0 at the centre, the others offset from it.
inline constexpr int kLeafClumpLobes = 3;

/// The hash both sides draw from (Chris Wellons' triple32-lite).
inline uint32_t leafClumpHash(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352du;
    x ^= x >> 15; x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

/// Lobe `lobe` of clump number `index`: its centre and radius.
void leafClumpLobe(const LeafClump& c, uint32_t index, int lobe, uint32_t seed,
                   bromath::Vec3& center, float& radius);

/// How many leaves clump number `index` carries (deterministic).
int leafClumpCount(const LeafClump& c, uint32_t index, const LeafClumpOptions& opts);

/// Leaf number `leaf` (its index over all clumps' leaves) of clump number
/// `index`: its 16-float instance record, the layout `LeafPlacements`
/// documents (+Z tip, +Y card normal, +X side, white tint). False when the
/// leaf has no extent.
bool leafInClump(const LeafClump& c, uint32_t index, uint32_t leaf, const LeafClumpOptions& opts,
                 float out[16]);

/// Every clump's leaves, clump by clump. `segment` holds each leaf's clump
/// index; `branchRadius` the clump radius; `branchDepth` 0.
LeafPlacements placeLeavesInClumps(const std::vector<LeafClump>& clumps,
                                   const LeafClumpOptions& opts = {});

} // namespace bromesh
