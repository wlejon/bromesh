#pragma once

#include "bromesh/optimization/spatial_hash.h"
#include "bromesh/procedural/space_colonization.h"
#include "bromesh/procedural/vec_math.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace bromesh {

/// Capsule (line segment with radius). `tag` is caller-supplied — typically
/// the index of the source `BranchSegment`, used so a leaf scattered along
/// segment N can be excluded from collision with the same segment N.
struct Capsule {
    Vec3  a{};
    Vec3  b{};
    float radius = 0.0f;
    int   tag    = -1;
};

/// Sphere obstacle (e.g. reserve volume around a bloom). `tag == -1` matches
/// no segment, so spheres are never auto-excluded by per-segment loops.
struct Sphere {
    Vec3  center{};
    float radius = 0.0f;
    int   tag    = -1;
};

/// Static occupancy lookup over a set of capsules and spheres.
///
/// Build once per scatter / grow pass. Queries (`distance`, `nearest`,
/// `intersectsSphere`) are O(k) where k is the number of capsules in cells
/// overlapping the query — proportional to the local crowding, not the total
/// count.
///
/// Capsules longer than the chosen cell size are inserted at multiple cells
/// along their axis so range queries stay correct without inflating the
/// global cell size.
///
/// Tag-based exclusion: every query takes an optional `excludeTag`. When
/// non-negative, capsules / spheres whose `tag` equals it are skipped — this
/// lets `placeLeavesOnBranches` ignore the candidate's own segment without
/// doing per-call rebuilds.
class CapsuleField {
public:
    CapsuleField();
    /// `cellSize == 0` (default) auto-picks from the input radii / lengths.
    explicit CapsuleField(std::vector<Capsule> capsules,
                          std::vector<Sphere>  spheres = {},
                          float cellSize = 0.0f);

    bool   empty()        const { return caps_.empty() && spheres_.empty(); }
    size_t capsuleCount() const { return caps_.size(); }
    size_t sphereCount()  const { return spheres_.size(); }
    float  cellSize()     const { return cellSize_; }

    /// True if `p` lies on or inside the capsule/sphere surface (with
    /// optional `extraClearance` inflating every obstacle by that much).
    bool   contains(Vec3 p, int excludeTag = -1, float extraClearance = 0.0f) const;

    /// True if signed distance to the nearest surface is <= clearance.
    /// Equivalent to `distance(p, exclude) <= clearance`.
    bool   tooClose(Vec3 p, float clearance, int excludeTag = -1) const;

    /// Signed distance to the nearest obstacle surface. Negative inside,
    /// positive outside. Returns +infinity when the field is empty.
    float  distance(Vec3 p, int excludeTag = -1) const;

    struct Nearest {
        Vec3  point{};        // closest point on the surface
        Vec3  normal{};       // outward unit normal at that point
        float distance = 0.0f;// signed; negative inside
        int   tag      = -1;  // tag of the nearest obstacle, or -1 if empty
    };
    Nearest nearest(Vec3 p, int excludeTag = -1) const;

    /// True iff a sphere centered at `c` with `radius` overlaps any obstacle.
    bool   intersectsSphere(Vec3 c, float radius, int excludeTag = -1) const;

    /// Build a capsule list from a `BranchSegment` vector. Each segment with
    /// `radius > 0` becomes one capsule with `tag = segment index`. Segments
    /// whose `from == to` (synthetic roots produced by `spaceColonize`) are
    /// skipped.
    static std::vector<Capsule> capsulesFromSegments(
        const std::vector<BranchSegment>& segs,
        float radiusScale = 1.0f);

private:
    void buildIndex();

    std::vector<Capsule> caps_;
    std::vector<Sphere>  spheres_;
    SpatialHash3D        capHash_;   // capsule-axis sample points -> caps_ index
    SpatialHash3D        sphHash_;   // sphere centers -> spheres_ index
    float                cellSize_   = 1.0f;
    float                queryPad_   = 0.0f; // half-cell + max radius slack
};

} // namespace bromesh
