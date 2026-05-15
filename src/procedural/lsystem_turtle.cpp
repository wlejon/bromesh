#include "bromesh/procedural/lsystem_turtle.h"

#include <bromath/bromath.h>

#include <cmath>

namespace bromesh {

using namespace bromath;

namespace {

constexpr float kDeg2Rad = bromath::DEG2RAD;

struct Pose {
    Vec3  position;
    Vec3  heading;
    Vec3  up;
    float radius;
    int   lastSeg;   // segment index that the next F should parent off
    int   rootSeg;   // synthetic root for this branch (-1 if none yet)
};

inline void reorthogonalize(Vec3& heading, Vec3& up) {
    heading = vnorm(heading);
    Vec3 right = vcross(heading, up);
    if (vdot(right, right) < 1e-12f) {
        // Heading collinear with up — pick any orthogonal as a fresh up.
        Vec3 alt = (std::fabs(heading.y) < 0.9f) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
        right = vnorm(vcross(heading, alt));
    } else {
        right = vnorm(right);
    }
    up = vnorm(vcross(right, heading));
}

inline float paramOr(const Module& m, size_t i, float fallback) {
    if (i < m.params.size()) return m.params[i];
    return fallback;
}

} // namespace

std::vector<BranchSegment> lsystemToBranches(
    const std::vector<Module>& modules,
    const TurtleOptions& opts) {
    std::vector<BranchSegment> segs;

    Pose cur;
    cur.position = opts.position;
    cur.heading  = vnormOr(opts.heading, {0, 1, 0});
    cur.up       = vnormOr(opts.up,      {0, 0, 1});
    reorthogonalize(cur.heading, cur.up);
    cur.radius   = opts.radius;
    cur.lastSeg  = -1;
    cur.rootSeg  = -1;

    std::vector<Pose> stack;

    auto emitRootIfNeeded = [&]() {
        if (cur.rootSeg >= 0) return;
        BranchSegment root;
        root.parent = -1;
        root.from = cur.position;
        root.to   = cur.position;
        root.radius = 0.0f;
        root.depth = 0;
        segs.push_back(root);
        cur.rootSeg = static_cast<int>(segs.size()) - 1;
        if (cur.lastSeg < 0) cur.lastSeg = cur.rootSeg;
    };

    for (const Module& m : modules) {
        switch (m.symbol) {
        case 'F': {
            emitRootIfNeeded();
            float len = paramOr(m, 0, opts.stepLength);
            float segR = paramOr(m, 1, cur.radius);
            Vec3 newPos = cur.position + cur.heading * len;
            BranchSegment s;
            s.parent = cur.lastSeg;
            s.from   = cur.position;
            s.to     = newPos;
            s.radius = segR;
            int parentDepth = (cur.lastSeg >= 0)
                ? segs[static_cast<size_t>(cur.lastSeg)].depth
                : 0;
            s.depth  = parentDepth + 1;
            segs.push_back(s);
            cur.lastSeg = static_cast<int>(segs.size()) - 1;
            cur.position = newPos;
            break;
        }
        case 'G':
        case 'f': {
            float len = paramOr(m, 0, opts.stepLength);
            cur.position = cur.position + cur.heading * len;
            break;
        }
        case '+':
        case '-':
        case '&':
        case '^':
        case '\\':
        case '/': {
            float ang = (!m.params.empty())
                ? m.params[0] * kDeg2Rad
                : opts.angle;
            Vec3 axis;
            switch (m.symbol) {
            case '+': axis = cur.up; break;
            case '-': axis = cur.up; ang = -ang; break;
            case '&': axis = vnorm(vcross(cur.heading, cur.up)); break;
            case '^': axis = vnorm(vcross(cur.heading, cur.up)); ang = -ang; break;
            case '\\': axis = cur.heading; break;
            case '/':  axis = cur.heading; ang = -ang; break;
            }
            Quat q = qaxisAngle(axis, ang);
            cur.heading = vnorm(qrotate(q, cur.heading));
            cur.up      = vnorm(qrotate(q, cur.up));
            reorthogonalize(cur.heading, cur.up);
            break;
        }
        case '|': {
            Quat q = qaxisAngle(cur.up, PI);
            cur.heading = vnorm(qrotate(q, cur.heading));
            reorthogonalize(cur.heading, cur.up);
            break;
        }
        case '!': {
            if (!m.params.empty()) cur.radius = m.params[0];
            break;
        }
        case '[': {
            stack.push_back(cur);
            break;
        }
        case ']': {
            if (!stack.empty()) {
                cur = stack.back();
                stack.pop_back();
            }
            break;
        }
        default:
            // Unknown symbols are silently ignored (annotations etc.).
            break;
        }
    }

    return segs;
}

} // namespace bromesh
