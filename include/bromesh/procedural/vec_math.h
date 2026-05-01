#pragma once

#include <cmath>
#include <cstddef>

namespace bromesh {

/// Minimal POD vec/quat types used by the procedural and sweep modules.
/// Kept local to bromesh to avoid pulling in a full math library; the
/// rest of bromesh stores positions as raw float arrays, and helpers in
/// this header bridge between the two styles.
struct Vec2 {
    float x = 0, y = 0;
};

struct Vec3 {
    float x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
};

struct Quat {
    float x = 0, y = 0, z = 0, w = 1;
    Quat() = default;
    Quat(float xx, float yy, float zz, float ww) : x(xx), y(yy), z(zz), w(ww) {}
};

inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 operator-(Vec3 a)         { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(Vec3 a, float s){ return {a.x*s, a.y*s, a.z*s}; }
inline Vec3 operator*(float s, Vec3 a){ return {a.x*s, a.y*s, a.z*s}; }
inline Vec3& operator+=(Vec3& a, Vec3 b) { a.x+=b.x; a.y+=b.y; a.z+=b.z; return a; }

inline float vdot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3  vcross(Vec3 a, Vec3 b) {
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
inline float vlen(Vec3 a) { return std::sqrt(vdot(a, a)); }
inline Vec3  vnorm(Vec3 a) {
    float L = vlen(a);
    return (L > 1e-20f) ? Vec3{a.x/L, a.y/L, a.z/L} : Vec3{0,0,0};
}
inline Vec3  vnormOr(Vec3 a, Vec3 fallback) {
    float L = vlen(a);
    return (L > 1e-20f) ? Vec3{a.x/L, a.y/L, a.z/L} : fallback;
}
inline float vdist(Vec3 a, Vec3 b) { return vlen(a - b); }
inline float vdist2(Vec3 a, Vec3 b) {
    Vec3 d = a - b;
    return vdot(d, d);
}

/// Build a quaternion that rotates unit vector `from` to unit vector `to`.
inline Quat quatFromTo(Vec3 from, Vec3 to) {
    float d = vdot(from, to);
    if (d > 0.999999f) return {0, 0, 0, 1};
    if (d < -0.999999f) {
        // 180-degree rotation around any axis perpendicular to `from`.
        Vec3 axis = vcross({1, 0, 0}, from);
        if (vdot(axis, axis) < 1e-12f) axis = vcross({0, 1, 0}, from);
        axis = vnorm(axis);
        return {axis.x, axis.y, axis.z, 0};
    }
    Vec3 c = vcross(from, to);
    float s = std::sqrt((1.0f + d) * 2.0f);
    float invs = 1.0f / s;
    return {c.x * invs, c.y * invs, c.z * invs, s * 0.5f};
}

inline Quat quatAxisAngle(Vec3 axis, float angle) {
    Vec3 a = vnorm(axis);
    float h = angle * 0.5f;
    float s = std::sin(h);
    return {a.x * s, a.y * s, a.z * s, std::cos(h)};
}

inline Quat quatMul(Quat a, Quat b) {
    return {
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    };
}

inline Vec3 quatRotate(Quat q, Vec3 v) {
    Vec3 u{q.x, q.y, q.z};
    float s = q.w;
    return u * (2.0f * vdot(u, v))
         + v * (s*s - vdot(u, u))
         + vcross(u, v) * (2.0f * s);
}

} // namespace bromesh
