#include "bromesh/analysis/convex_decomposition.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace bromesh {

namespace {

struct V3 {
    double x, y, z;
};

V3 sub(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
V3 cross(const V3& a, const V3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double len(const V3& a) { return std::sqrt(dot(a, a)); }

struct Face {
    std::array<uint32_t, 3> v;
    V3 n;
    double d;
    std::vector<uint32_t> outside;
    bool alive = true;
};

uint64_t edgeKey(uint32_t a, uint32_t b) { return (uint64_t(a) << 32) | b; }

class QuickHull {
public:
    QuickHull(const std::vector<V3>& pts, double eps) : P(pts), eps(eps) {}

    bool build() {
        std::array<uint32_t, 4> seed;
        if (!initialSimplex(seed)) return false;

        V3 centroid{0, 0, 0};
        for (uint32_t i : seed) {
            centroid.x += P[i].x * 0.25;
            centroid.y += P[i].y * 0.25;
            centroid.z += P[i].z * 0.25;
        }
        const std::array<std::array<int, 3>, 4> tris = {{{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}};
        std::vector<uint32_t> fresh;
        for (const auto& t : tris) {
            uint32_t a = seed[t[0]], b = seed[t[1]], c = seed[t[2]];
            V3 n = cross(sub(P[b], P[a]), sub(P[c], P[a]));
            if (dot(n, sub(centroid, P[a])) > 0) std::swap(b, c);
            fresh.push_back(addFace(a, b, c));
        }

        std::vector<uint32_t> all(P.size());
        for (uint32_t i = 0; i < P.size(); ++i) all[i] = i;
        assign(all, fresh);

        for (;;) {
            int pick = -1;
            for (size_t f = 0; f < F.size(); ++f)
                if (F[f].alive && !F[f].outside.empty()) { pick = static_cast<int>(f); break; }
            if (pick < 0) break;
            addPoint(static_cast<uint32_t>(pick));
        }
        return true;
    }

    MeshData toMesh() const {
        MeshData m;
        for (const Face& f : F) {
            if (!f.alive) continue;
            uint32_t base = static_cast<uint32_t>(m.vertexCount());
            for (uint32_t vi : f.v) {
                m.positions.push_back(static_cast<float>(P[vi].x));
                m.positions.push_back(static_cast<float>(P[vi].y));
                m.positions.push_back(static_cast<float>(P[vi].z));
                m.normals.push_back(static_cast<float>(f.n.x));
                m.normals.push_back(static_cast<float>(f.n.y));
                m.normals.push_back(static_cast<float>(f.n.z));
            }
            m.indices.push_back(base);
            m.indices.push_back(base + 1);
            m.indices.push_back(base + 2);
        }
        return m;
    }

private:
    const std::vector<V3>& P;
    double eps;
    std::vector<Face> F;
    std::unordered_map<uint64_t, uint32_t> edges;

    double dist(const Face& f, uint32_t p) const { return dot(f.n, P[p]) - f.d; }

    uint32_t addFace(uint32_t a, uint32_t b, uint32_t c) {
        Face f;
        f.v = {a, b, c};
        V3 n = cross(sub(P[b], P[a]), sub(P[c], P[a]));
        double l = len(n);
        f.n = l > 0 ? V3{n.x / l, n.y / l, n.z / l} : V3{0, 0, 0};
        f.d = dot(f.n, P[a]);
        uint32_t id = static_cast<uint32_t>(F.size());
        F.push_back(std::move(f));
        edges[edgeKey(a, b)] = id;
        edges[edgeKey(b, c)] = id;
        edges[edgeKey(c, a)] = id;
        return id;
    }

    void assign(const std::vector<uint32_t>& pts, const std::vector<uint32_t>& faces) {
        for (uint32_t p : pts) {
            double best = eps;
            int bestF = -1;
            for (uint32_t fi : faces) {
                double d = dist(F[fi], p);
                if (d > best) { best = d; bestF = static_cast<int>(fi); }
            }
            if (bestF >= 0) F[bestF].outside.push_back(p);
        }
    }

    void addPoint(uint32_t fi) {
        uint32_t eye = F[fi].outside[0];
        double far = dist(F[fi], eye);
        for (uint32_t p : F[fi].outside) {
            double d = dist(F[fi], p);
            if (d > far) { far = d; eye = p; }
        }

        std::vector<uint32_t> visible;
        std::vector<char> seen(F.size(), 0);
        std::vector<uint32_t> stack{fi};
        seen[fi] = 1;
        while (!stack.empty()) {
            uint32_t f = stack.back();
            stack.pop_back();
            visible.push_back(f);
            for (int e = 0; e < 3; ++e) {
                uint32_t a = F[f].v[e], b = F[f].v[(e + 1) % 3];
                auto it = edges.find(edgeKey(b, a));
                if (it == edges.end()) continue;
                uint32_t g = it->second;
                if (seen[g] || !F[g].alive) continue;
                if (dist(F[g], eye) > eps) {
                    seen[g] = 1;
                    stack.push_back(g);
                }
            }
        }

        std::vector<std::pair<uint32_t, uint32_t>> horizon;
        std::vector<char> isVisible(F.size(), 0);
        for (uint32_t f : visible) isVisible[f] = 1;
        for (uint32_t f : visible) {
            for (int e = 0; e < 3; ++e) {
                uint32_t a = F[f].v[e], b = F[f].v[(e + 1) % 3];
                auto it = edges.find(edgeKey(b, a));
                if (it == edges.end() || !isVisible[it->second]) horizon.emplace_back(a, b);
            }
        }

        std::vector<uint32_t> orphans;
        for (uint32_t f : visible) {
            Face& face = F[f];
            face.alive = false;
            for (uint32_t p : face.outside)
                if (p != eye) orphans.push_back(p);
            face.outside.clear();
            face.outside.shrink_to_fit();
            for (int e = 0; e < 3; ++e) {
                uint32_t a = face.v[e], b = face.v[(e + 1) % 3];
                auto it = edges.find(edgeKey(a, b));
                if (it != edges.end() && it->second == f) edges.erase(it);
            }
        }

        std::vector<uint32_t> fresh;
        fresh.reserve(horizon.size());
        for (const auto& [a, b] : horizon) fresh.push_back(addFace(a, b, eye));
        assign(orphans, fresh);
    }

    bool initialSimplex(std::array<uint32_t, 4>& s) const {
        const size_t n = P.size();
        std::array<uint32_t, 6> ext{};
        for (uint32_t i = 0; i < n; ++i) {
            if (P[i].x < P[ext[0]].x) ext[0] = i;
            if (P[i].x > P[ext[1]].x) ext[1] = i;
            if (P[i].y < P[ext[2]].y) ext[2] = i;
            if (P[i].y > P[ext[3]].y) ext[3] = i;
            if (P[i].z < P[ext[4]].z) ext[4] = i;
            if (P[i].z > P[ext[5]].z) ext[5] = i;
        }
        double best = -1;
        for (int i = 0; i < 6; ++i)
            for (int j = i + 1; j < 6; ++j) {
                double d = len(sub(P[ext[i]], P[ext[j]]));
                if (d > best) { best = d; s[0] = ext[i]; s[1] = ext[j]; }
            }
        if (best <= eps) return false;

        V3 axis = sub(P[s[1]], P[s[0]]);
        best = eps;
        bool found = false;
        for (uint32_t i = 0; i < n; ++i) {
            double d = len(cross(axis, sub(P[i], P[s[0]]))) / len(axis);
            if (d > best) { best = d; s[2] = i; found = true; }
        }
        if (!found) return false;

        V3 nrm = cross(sub(P[s[1]], P[s[0]]), sub(P[s[2]], P[s[0]]));
        double nl = len(nrm);
        best = eps;
        found = false;
        for (uint32_t i = 0; i < n; ++i) {
            double d = std::abs(dot(nrm, sub(P[i], P[s[0]]))) / nl;
            if (d > best) { best = d; s[3] = i; found = true; }
        }
        return found;
    }
};

}  // namespace

MeshData convexHull(const MeshData& mesh) {
    const size_t n = mesh.vertexCount();
    if (n < 4) return {};
    std::vector<V3> pts(n);
    V3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (size_t i = 0; i < n; ++i) {
        pts[i] = {mesh.positions[i * 3], mesh.positions[i * 3 + 1], mesh.positions[i * 3 + 2]};
        lo = {std::min(lo.x, pts[i].x), std::min(lo.y, pts[i].y), std::min(lo.z, pts[i].z)};
        hi = {std::max(hi.x, pts[i].x), std::max(hi.y, pts[i].y), std::max(hi.z, pts[i].z)};
    }
    const double scale = std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z});
    if (!(scale > 0) || !std::isfinite(scale)) return {};
    QuickHull qh(pts, scale * 1e-6 + 1e-12);
    if (!qh.build()) return {};
    return qh.toMesh();
}

}  // namespace bromesh
