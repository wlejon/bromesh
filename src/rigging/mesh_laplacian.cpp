#include "bromesh/rigging/mesh_laplacian.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>

namespace bromesh {

namespace {

struct V3 { double x, y, z; };
V3 operator-(V3 a, V3 b) { return { a.x-b.x, a.y-b.y, a.z-b.z }; }
double dot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
V3 cross(V3 a, V3 b) {
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
double norm(V3 a) { return std::sqrt(dot(a, a)); }

V3 vpos(const MeshData& m, int i) {
    return { (double)m.positions[i*3 + 0],
             (double)m.positions[i*3 + 1],
             (double)m.positions[i*3 + 2] };
}

// Cotangent of the angle at vertex opposite to edge (p,q) in triangle (opp,p,q).
// cot θ = cos θ / sin θ = (u·v) / |u x v|.
double cotAt(V3 opp, V3 p, V3 q) {
    V3 u = p - opp;
    V3 v = q - opp;
    double c = dot(u, v);
    double s = norm(cross(u, v));
    if (s < 1e-20) return 0.0;
    double ct = c / s;
    // Clamp to [0, big]: negative cotangents (obtuse corner) drop to 0.
    // Positive clamp prevents pathological triangles from blowing up the
    // condition number.
    if (ct < 0.0) ct = 0.0;
    if (ct > 1e6)  ct = 1e6;
    return ct;
}

// Mixed Voronoi area contribution of vertex v0 in triangle (v0, v1, v2).
// For non-obtuse triangles, returns the Voronoi region of v0. For obtuse
// triangles, returns triArea/2 if v0 is the obtuse corner else triArea/4
// (Meyer et al. 2002).
double mixedVoronoiArea(V3 v0, V3 v1, V3 v2) {
    V3 e01 = v1 - v0;
    V3 e02 = v2 - v0;
    V3 e12 = v2 - v1;
    double triArea = 0.5 * norm(cross(e01, e02));
    if (triArea < 1e-20) return 0.0;

    // Corner angles (law of cosines via dot products).
    double d00 = dot(e01, e01);      // |e01|^2
    double d11 = dot(e02, e02);      // |e02|^2
    double d22 = dot(e12, e12);      // |e12|^2

    // cos at v0 = (|e01|^2 + |e02|^2 - |e12|^2) / (2|e01||e02|)
    // Easier: check obtuseness by sign of dot products.
    double c0 = dot(e01, e02);
    double c1 = dot(v0 - v1, v2 - v1);
    double c2 = dot(v0 - v2, v1 - v2);

    bool obt0 = c0 < 0.0;
    bool obt1 = c1 < 0.0;
    bool obt2 = c2 < 0.0;
    bool anyObtuse = obt0 || obt1 || obt2;

    if (anyObtuse) {
        return obt0 ? triArea * 0.5 : triArea * 0.25;
    }
    // Voronoi area at v0:
    //   (1/8) [ |e01|^2 cot(angle at v2) + |e02|^2 cot(angle at v1) ]
    double cotAtV1 = cotAt(v1, v0, v2);
    double cotAtV2 = cotAt(v2, v0, v1);
    double vor = 0.125 * (d00 * cotAtV2 + d11 * cotAtV1);
    if (vor < 0.0) vor = 0.0;
    return vor;
}

} // namespace

void assembleCotangentLaplacian(const MeshData& mesh,
                                SparseCsr& L,
                                std::vector<double>& massDiag) {
    const int n = (int)mesh.vertexCount();
    const auto& idx = mesh.indices;

    // Accumulate edge weights into a map keyed by (min, max). Only one entry
    // per undirected edge; we re-emit symmetric rows at CSR build.
    // Using vector<map> is simpler and safe for modest meshes; could be
    // replaced by a sorted vector-of-pairs later if needed.
    std::vector<std::unordered_map<int, double>> edgeWeight(n);
    massDiag.assign(n, 0.0);

    for (size_t t = 0; t + 2 < idx.size(); t += 3) {
        int i0 = (int)idx[t + 0];
        int i1 = (int)idx[t + 1];
        int i2 = (int)idx[t + 2];
        V3 p0 = vpos(mesh, i0);
        V3 p1 = vpos(mesh, i1);
        V3 p2 = vpos(mesh, i2);

        // Cotangent at each corner contributes to the opposite edge weight.
        double cot0 = cotAt(p0, p1, p2); // opposite edge (1,2)
        double cot1 = cotAt(p1, p2, p0); // opposite edge (2,0)
        double cot2 = cotAt(p2, p0, p1); // opposite edge (0,1)

        auto add = [&](int a, int b, double w) {
            edgeWeight[a][b] += w;
            edgeWeight[b][a] += w;
        };
        // Standard factor 0.5 per cotangent contribution.
        add(i1, i2, 0.5 * cot0);
        add(i2, i0, 0.5 * cot1);
        add(i0, i1, 0.5 * cot2);

        // Mass contributions.
        massDiag[i0] += mixedVoronoiArea(p0, p1, p2);
        massDiag[i1] += mixedVoronoiArea(p1, p2, p0);
        massDiag[i2] += mixedVoronoiArea(p2, p0, p1);
    }

    // Avoid zero mass (isolated / degenerate verts).
    for (int i = 0; i < n; ++i) {
        if (massDiag[i] < 1e-20) massDiag[i] = 1e-20;
    }

    // Build CSR: L[i,i] = -sum_j w_ij ; L[i,j] = w_ij for j != i.
    // Sign convention: (Lf)[i] = sum_j w_ij (f[j] - f[i]) = (off-diag · f) + diag * f[i].
    // i.e. diagonal is -sum of row off-diagonals. This L is negative semi-definite.
    //
    // Downstream solves use -L (positive semi-definite). Expose L as-is;
    // callers negate as needed.
    L.rows = n;
    L.cols = n;
    L.rowStart.assign(n + 1, 0);
    // Count nnz per row.
    for (int i = 0; i < n; ++i) {
        L.rowStart[i + 1] = (int)edgeWeight[i].size() + 1; // +1 for diagonal
    }
    for (int i = 0; i < n; ++i) L.rowStart[i + 1] += L.rowStart[i];
    int nnz = L.rowStart[n];
    L.colIndex.assign(nnz, 0);
    L.values.assign(nnz, 0.0);

    for (int i = 0; i < n; ++i) {
        int base = L.rowStart[i];
        // Gather neighbors, sort by column for a canonical layout.
        std::vector<std::pair<int, double>> row;
        row.reserve(edgeWeight[i].size() + 1);
        double diag = 0.0;
        for (auto& kv : edgeWeight[i]) {
            row.emplace_back(kv.first, kv.second);
            diag -= kv.second;
        }
        row.emplace_back(i, diag);
        std::sort(row.begin(), row.end(),
                  [](const auto& a, const auto& b){ return a.first < b.first; });
        for (size_t k = 0; k < row.size(); ++k) {
            L.colIndex[base + k] = row[k].first;
            L.values[base + k]   = row[k].second;
        }
    }
}

void sparseSpmv(const SparseCsr& A,
                const std::vector<double>& x,
                std::vector<double>& y) {
    const int n = A.rows;
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        int b = A.rowStart[i];
        int e = A.rowStart[i + 1];
        for (int k = b; k < e; ++k) s += A.values[k] * x[A.colIndex[k]];
        y[i] = s;
    }
}

int sparseCg(const SparseCsr& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double tol,
             int maxIter) {
    const int n = A.rows;
    std::vector<double> r(n), z(n), p(n), Ap(n);
    // Jacobi preconditioner: inverse of diagonal.
    std::vector<double> invDiag(n, 1.0);
    for (int i = 0; i < n; ++i) {
        int bs = A.rowStart[i], es = A.rowStart[i + 1];
        double d = 0.0;
        for (int k = bs; k < es; ++k) if (A.colIndex[k] == i) { d = A.values[k]; break; }
        if (std::fabs(d) > 1e-20) invDiag[i] = 1.0 / d;
    }

    // r = b - A x
    sparseSpmv(A, x, Ap);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    double bnorm2 = 0.0;
    for (int i = 0; i < n; ++i) bnorm2 += b[i] * b[i];
    if (bnorm2 < 1e-30) return 0;

    // z = M^-1 r
    for (int i = 0; i < n; ++i) z[i] = invDiag[i] * r[i];
    p = z;
    double rz = 0.0;
    for (int i = 0; i < n; ++i) rz += r[i] * z[i];

    double tol2 = tol * tol * bnorm2;
    int it = 0;
    for (; it < maxIter; ++it) {
        sparseSpmv(A, p, Ap);
        double pAp = 0.0;
        for (int i = 0; i < n; ++i) pAp += p[i] * Ap[i];
        if (std::fabs(pAp) < 1e-30) break;
        double alpha = rz / pAp;
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        double rnorm2 = 0.0;
        for (int i = 0; i < n; ++i) rnorm2 += r[i] * r[i];
        if (rnorm2 < tol2) { ++it; break; }
        for (int i = 0; i < n; ++i) z[i] = invDiag[i] * r[i];
        double rzNew = 0.0;
        for (int i = 0; i < n; ++i) rzNew += r[i] * z[i];
        double beta = rzNew / rz;
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz = rzNew;
    }
    return it;
}

} // namespace bromesh
