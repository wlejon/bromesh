#include "bromesh/procedural/obstacle_field.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace bromesh {

namespace {

// Closest point on segment [a,b] to query point p. Returns the parametric
// t in [0,1] such that closest = a + t*(b-a).
inline Vec3 closestOnSegment(Vec3 a, Vec3 b, Vec3 p, float& tOut) {
    Vec3 ab = b - a;
    float L2 = vdot(ab, ab);
    if (L2 < 1e-20f) { tOut = 0.0f; return a; }
    float t = vdot(p - a, ab) / L2;
    if (t < 0.0f) t = 0.0f;
    else if (t > 1.0f) t = 1.0f;
    tOut = t;
    return a + ab * t;
}

} // namespace

CapsuleField::CapsuleField() = default;

CapsuleField::CapsuleField(std::vector<Capsule> capsules,
                           std::vector<Sphere>  spheres,
                           float cellSize)
    : caps_(std::move(capsules)),
      spheres_(std::move(spheres)) {
    // Auto-size: 2 * (largest radius + half-length of the longest capsule).
    // Floor at 0.05 so empty / tiny inputs still produce a sane grid.
    float maxR     = 0.0f;
    float maxHalfL = 0.0f;
    for (const Capsule& c : caps_) {
        if (c.radius > maxR) maxR = c.radius;
        float h = vlen(c.b - c.a) * 0.5f;
        if (h > maxHalfL) maxHalfL = h;
    }
    for (const Sphere& s : spheres_) {
        if (s.radius > maxR) maxR = s.radius;
    }

    if (cellSize > 0.0f) {
        cellSize_ = cellSize;
    } else {
        cellSize_ = 2.0f * (maxR + maxHalfL);
        if (cellSize_ < 0.05f) cellSize_ = 0.05f;
    }
    queryPad_ = maxR;

    buildIndex();
}

void CapsuleField::buildIndex() {
    capHash_.reset(cellSize_);
    sphHash_.reset(cellSize_);

    // For each capsule, sample points along the axis no more than `cellSize_`
    // apart so a radius query at any point on or near the capsule is
    // guaranteed to find it via at least one sample's cell. We register the
    // SAME capsule id at each sample point; the nearest()/distance() pass
    // dedupes results internally.
    for (size_t i = 0; i < caps_.size(); ++i) {
        const Capsule& c = caps_[i];
        Vec3 ab = c.b - c.a;
        float L = vlen(ab);
        // Step a hair smaller than cellSize so adjacent cells overlap; this
        // closes the corner-of-cell gap when the capsule axis grazes a cell
        // boundary.
        float step = (cellSize_ > 0.0f) ? cellSize_ * 0.95f : 1.0f;
        if (L < step) {
            // Single sample at the midpoint covers the whole capsule.
            Vec3 mid = c.a + ab * 0.5f;
            capHash_.insert(mid, static_cast<int32_t>(i));
            continue;
        }
        int n = static_cast<int>(std::ceil(L / step));
        if (n < 1) n = 1;
        for (int s = 0; s <= n; ++s) {
            float t = (n == 0) ? 0.0f : (static_cast<float>(s) / static_cast<float>(n));
            Vec3 p = c.a + ab * t;
            capHash_.insert(p, static_cast<int32_t>(i));
        }
    }
    for (size_t i = 0; i < spheres_.size(); ++i) {
        sphHash_.insert(spheres_[i].center, static_cast<int32_t>(i));
    }
}

// Common kernel: walk every candidate capsule/sphere within a search radius
// and accumulate the running-best signed distance + nearest surface info.
//
// `searchRadius` only needs to contain the closest obstacle. We start with
// (cellSize + maxRadius) which is always large enough to find any obstacle
// whose surface could plausibly be the nearest — capsule axis samples are
// at most one cellSize apart along the axis, plus the capsule's radius
// reaches outward. If nothing comes back we expand to infinity (fall back
// to a brute scan), which covers the rare case of a tiny field built with a
// large user-supplied cellSize.
namespace {

struct BestHit {
    bool   any         = false;
    float  bestDist    = std::numeric_limits<float>::infinity();
    Vec3   bestPoint{};
    Vec3   bestNormal{0, 1, 0};
    int    bestTag     = -1;
};

inline void considerCapsule(const Capsule& c, Vec3 p, BestHit& h) {
    float t = 0.0f;
    Vec3 cp = closestOnSegment(c.a, c.b, p, t);
    Vec3 d  = p - cp;
    float L = vlen(d);
    Vec3 n;
    if (L > 1e-8f) {
        n = d * (1.0f / L);
    } else {
        // p sits on the axis: pick any perpendicular as the outward normal.
        Vec3 ab = c.b - c.a;
        float ablen = vlen(ab);
        Vec3 axis = (ablen > 1e-8f) ? ab * (1.0f / ablen) : Vec3{0, 1, 0};
        Vec3 alt = (std::fabs(axis.y) < 0.9f) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
        n = vnorm(vcross(axis, alt));
    }
    float signed_d = L - c.radius;
    if (signed_d < h.bestDist) {
        h.any        = true;
        h.bestDist   = signed_d;
        h.bestPoint  = cp + n * c.radius;
        h.bestNormal = n;
        h.bestTag    = c.tag;
    }
}

inline void considerSphere(const Sphere& s, Vec3 p, BestHit& h) {
    Vec3 d = p - s.center;
    float L = vlen(d);
    Vec3 n = (L > 1e-8f) ? d * (1.0f / L) : Vec3{0, 1, 0};
    float signed_d = L - s.radius;
    if (signed_d < h.bestDist) {
        h.any        = true;
        h.bestDist   = signed_d;
        h.bestPoint  = s.center + n * s.radius;
        h.bestNormal = n;
        h.bestTag    = s.tag;
    }
}

} // namespace

CapsuleField::Nearest CapsuleField::nearest(Vec3 p, int excludeTag) const {
    Nearest out;
    out.distance = std::numeric_limits<float>::infinity();
    if (caps_.empty() && spheres_.empty()) return out;

    BestHit h;

    auto query = [&](float radius) {
        std::vector<int32_t> ids;

        if (!caps_.empty()) {
            ids.clear();
            capHash_.radiusQuery(p, radius, ids);
            std::sort(ids.begin(), ids.end());
            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
            for (int32_t id : ids) {
                const Capsule& c = caps_[static_cast<size_t>(id)];
                if (excludeTag != -1 && c.tag == excludeTag) continue;
                considerCapsule(c, p, h);
            }
        }
        if (!spheres_.empty()) {
            ids.clear();
            sphHash_.radiusQuery(p, radius, ids);
            for (int32_t id : ids) {
                const Sphere& s = spheres_[static_cast<size_t>(id)];
                if (excludeTag != -1 && s.tag == excludeTag) continue;
                considerSphere(s, p, h);
            }
        }
    };

    // First pass: tight radius covering one cell + the largest capsule radius.
    query(cellSize_ + queryPad_);

    // If still nothing hit, the field is sparse relative to the grid — fall
    // back to a brute scan. Plant scenes are small (hundreds of segments) so
    // this is fine as a safety net and only kicks in for far-away points.
    if (!h.any) {
        for (const Capsule& c : caps_) {
            if (excludeTag != -1 && c.tag == excludeTag) continue;
            considerCapsule(c, p, h);
        }
        for (const Sphere& s : spheres_) {
            if (excludeTag != -1 && s.tag == excludeTag) continue;
            considerSphere(s, p, h);
        }
    }

    if (h.any) {
        out.point    = h.bestPoint;
        out.normal   = h.bestNormal;
        out.distance = h.bestDist;
        out.tag      = h.bestTag;
    }
    return out;
}

float CapsuleField::distance(Vec3 p, int excludeTag) const {
    return nearest(p, excludeTag).distance;
}

bool CapsuleField::contains(Vec3 p, int excludeTag, float extraClearance) const {
    return distance(p, excludeTag) <= extraClearance;
}

bool CapsuleField::tooClose(Vec3 p, float clearance, int excludeTag) const {
    return distance(p, excludeTag) <= clearance;
}

bool CapsuleField::intersectsSphere(Vec3 c, float radius, int excludeTag) const {
    return distance(c, excludeTag) < radius;
}

std::vector<Capsule> CapsuleField::capsulesFromSegments(
    const std::vector<BranchSegment>& segs,
    float radiusScale) {
    std::vector<Capsule> out;
    out.reserve(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        const BranchSegment& s = segs[i];
        if (s.radius <= 0.0f) continue;
        // Synthetic roots from spaceColonize have from == to; skip — they
        // can't act as a capsule and would degenerate to a sphere we don't
        // want.
        if (vdist2(s.from, s.to) < 1e-14f) continue;
        Capsule c;
        c.a = s.from;
        c.b = s.to;
        c.radius = s.radius * radiusScale;
        c.tag = static_cast<int>(i);
        out.push_back(c);
    }
    return out;
}

} // namespace bromesh
