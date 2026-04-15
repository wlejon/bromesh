#include "bromesh/rigging/bbw.h"
#include "bromesh/rigging/mesh_laplacian.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

#if BROMESH_HAS_OSQP
extern "C" {
#include "osqp.h"
}
#endif

namespace bromesh {

namespace {

using Vec3 = std::array<double, 3>;

void matInverse4(const float* m, double* inv) {
    double a[16];
    for (int i = 0; i < 16; ++i) a[i] = (double)m[i];
    inv[0]  =  a[5]*a[10]*a[15] - a[5]*a[11]*a[14] - a[9]*a[6]*a[15] + a[9]*a[7]*a[14] + a[13]*a[6]*a[11] - a[13]*a[7]*a[10];
    inv[4]  = -a[4]*a[10]*a[15] + a[4]*a[11]*a[14] + a[8]*a[6]*a[15] - a[8]*a[7]*a[14] - a[12]*a[6]*a[11] + a[12]*a[7]*a[10];
    inv[8]  =  a[4]*a[9]*a[15]  - a[4]*a[11]*a[13] - a[8]*a[5]*a[15] + a[8]*a[7]*a[13] + a[12]*a[5]*a[11] - a[12]*a[7]*a[9];
    inv[12] = -a[4]*a[9]*a[14]  + a[4]*a[10]*a[13] + a[8]*a[5]*a[14] - a[8]*a[6]*a[13] - a[12]*a[5]*a[10] + a[12]*a[6]*a[9];
    inv[1]  = -a[1]*a[10]*a[15] + a[1]*a[11]*a[14] + a[9]*a[2]*a[15] - a[9]*a[3]*a[14] - a[13]*a[2]*a[11] + a[13]*a[3]*a[10];
    inv[5]  =  a[0]*a[10]*a[15] - a[0]*a[11]*a[14] - a[8]*a[2]*a[15] + a[8]*a[3]*a[14] + a[12]*a[2]*a[11] - a[12]*a[3]*a[10];
    inv[9]  = -a[0]*a[9]*a[15]  + a[0]*a[11]*a[13] + a[8]*a[1]*a[15] - a[8]*a[3]*a[13] - a[12]*a[1]*a[11] + a[12]*a[3]*a[9];
    inv[13] =  a[0]*a[9]*a[14]  - a[0]*a[10]*a[13] - a[8]*a[1]*a[14] + a[8]*a[2]*a[13] + a[12]*a[1]*a[10] - a[12]*a[2]*a[9];
    inv[2]  =  a[1]*a[6]*a[15]  - a[1]*a[7]*a[14]  - a[5]*a[2]*a[15] + a[5]*a[3]*a[14] + a[13]*a[2]*a[7]  - a[13]*a[3]*a[6];
    inv[6]  = -a[0]*a[6]*a[15]  + a[0]*a[7]*a[14]  + a[4]*a[2]*a[15] - a[4]*a[3]*a[14] - a[12]*a[2]*a[7]  + a[12]*a[3]*a[6];
    inv[10] =  a[0]*a[5]*a[15]  - a[0]*a[7]*a[13]  - a[4]*a[1]*a[15] + a[4]*a[3]*a[13] + a[12]*a[1]*a[7]  - a[12]*a[3]*a[5];
    inv[14] = -a[0]*a[5]*a[14]  + a[0]*a[6]*a[13]  + a[4]*a[1]*a[14] - a[4]*a[2]*a[13] - a[12]*a[1]*a[6]  + a[12]*a[2]*a[5];
    inv[3]  = -a[1]*a[6]*a[11]  + a[1]*a[7]*a[10]  + a[5]*a[2]*a[11] - a[5]*a[3]*a[10] - a[9]*a[2]*a[7]   + a[9]*a[3]*a[6];
    inv[7]  =  a[0]*a[6]*a[11]  - a[0]*a[7]*a[10]  - a[4]*a[2]*a[11] + a[4]*a[3]*a[10] + a[8]*a[2]*a[7]   - a[8]*a[3]*a[6];
    inv[11] = -a[0]*a[5]*a[11]  + a[0]*a[7]*a[9]   + a[4]*a[1]*a[11] - a[4]*a[3]*a[9]  - a[8]*a[1]*a[7]   + a[8]*a[3]*a[5];
    inv[15] =  a[0]*a[5]*a[10]  - a[0]*a[6]*a[9]   - a[4]*a[1]*a[10] + a[4]*a[2]*a[9]  + a[8]*a[1]*a[6]   - a[8]*a[2]*a[5];
    double det = a[0]*inv[0] + a[1]*inv[4] + a[2]*inv[8] + a[3]*inv[12];
    if (std::fabs(det) < 1e-20) {
        for (int i = 0; i < 16; ++i) inv[i] = (i % 5 == 0) ? 1.0 : 0.0;
        return;
    }
    double id = 1.0 / det;
    for (int i = 0; i < 16; ++i) inv[i] *= id;
}

double distPointSegment(Vec3 p, Vec3 a, Vec3 b) {
    double dx = b[0]-a[0], dy = b[1]-a[1], dz = b[2]-a[2];
    double L2 = dx*dx + dy*dy + dz*dz;
    if (L2 < 1e-20) {
        double ex = p[0]-a[0], ey = p[1]-a[1], ez = p[2]-a[2];
        return std::sqrt(ex*ex + ey*ey + ez*ez);
    }
    double t = ((p[0]-a[0])*dx + (p[1]-a[1])*dy + (p[2]-a[2])*dz) / L2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double cx = a[0] + t*dx, cy = a[1] + t*dy, cz = a[2] + t*dz;
    double ex = p[0]-cx, ey = p[1]-cy, ez = p[2]-cz;
    return std::sqrt(ex*ex + ey*ey + ez*ez);
}

void extractBoneHeads(const Skeleton& skeleton,
                      std::vector<Vec3>& heads,
                      std::vector<Vec3>& tails) {
    const size_t n = skeleton.bones.size();
    heads.resize(n);
    tails.resize(n);
    for (size_t i = 0; i < n; ++i) {
        double w[16];
        matInverse4(skeleton.bones[i].inverseBind, w);
        heads[i] = { w[12], w[13], w[14] };
    }
    for (size_t i = 0; i < n; ++i) {
        int firstChild = -1;
        for (size_t j = 0; j < n; ++j) {
            if (skeleton.bones[j].parent == (int)i) { firstChild = (int)j; break; }
        }
        if (firstChild >= 0) {
            tails[i] = heads[firstChild];
        } else {
            double w[16];
            matInverse4(skeleton.bones[i].inverseBind, w);
            Vec3 yAxis = { w[4], w[5], w[6] };
            double len = 1.0;
            int p = skeleton.bones[i].parent;
            if (p >= 0) {
                double dx = heads[i][0]-heads[p][0];
                double dy = heads[i][1]-heads[p][1];
                double dz = heads[i][2]-heads[p][2];
                double d = std::sqrt(dx*dx+dy*dy+dz*dz);
                if (d > 1e-8) len = d;
            }
            tails[i] = { heads[i][0] + yAxis[0]*len,
                         heads[i][1] + yAxis[1]*len,
                         heads[i][2] + yAxis[2]*len };
        }
    }
}

// Compute Q = L * M^-1 * L as a sparse symmetric matrix (CSR). L must be
// sorted per-row. M is diagonal (mass[]). We emit a symmetric CSR (both
// halves) and let the OSQP conversion filter the upper triangle.
void buildBiharmonic(const SparseCsr& L,
                     const std::vector<double>& mass,
                     SparseCsr& Q) {
    const int n = L.rows;
    Q.rows = n;
    Q.cols = n;
    Q.rowStart.assign(n + 1, 0);

    // Two-pass: first count per-row nnz using a dense "seen" mask.
    std::vector<int> rowNnz(n, 0);
    std::vector<int> marker(n, -1);
    for (int i = 0; i < n; ++i) {
        int count = 0;
        int bi = L.rowStart[i], ei = L.rowStart[i + 1];
        for (int a = bi; a < ei; ++a) {
            int k = L.colIndex[a];
            int bk = L.rowStart[k], ek = L.rowStart[k + 1];
            for (int c = bk; c < ek; ++c) {
                int j = L.colIndex[c];
                if (marker[j] != i) { marker[j] = i; ++count; }
            }
        }
        rowNnz[i] = count;
    }
    for (int i = 0; i < n; ++i) Q.rowStart[i + 1] = Q.rowStart[i] + rowNnz[i];
    int nnz = Q.rowStart[n];
    Q.colIndex.assign(nnz, 0);
    Q.values.assign(nnz, 0.0);

    // Second pass: accumulate values with a dense scratch vector.
    std::vector<double> scratch(n, 0.0);
    std::vector<int> touched; touched.reserve(64);
    std::fill(marker.begin(), marker.end(), -1);
    for (int i = 0; i < n; ++i) {
        int bi = L.rowStart[i], ei = L.rowStart[i + 1];
        touched.clear();
        for (int a = bi; a < ei; ++a) {
            int k = L.colIndex[a];
            double Lik = L.values[a];
            double invMk = (mass[k] > 1e-30) ? (1.0 / mass[k]) : 0.0;
            double scale = Lik * invMk;
            int bk = L.rowStart[k], ek = L.rowStart[k + 1];
            for (int c = bk; c < ek; ++c) {
                int j = L.colIndex[c];
                double Lkj = L.values[c];
                if (marker[j] != i) { marker[j] = i; scratch[j] = 0.0; touched.push_back(j); }
                scratch[j] += scale * Lkj;
            }
        }
        // Sort touched by column for canonical row layout.
        std::sort(touched.begin(), touched.end());
        int base = Q.rowStart[i];
        for (size_t t = 0; t < touched.size(); ++t) {
            int j = touched[t];
            Q.colIndex[base + (int)t] = j;
            Q.values[base + (int)t]   = scratch[j];
        }
    }
}

} // namespace

#if BROMESH_HAS_OSQP

// Build an OSQP-compatible CSC of the upper triangle of a symmetric CSR.
// In the symmetric CSR (i,j) == (j,i), so CSC(upper) has columns where
// we collect only rows i ≤ j.
static void symCsrToUpperCsc(const SparseCsr& Asym,
                             std::vector<OSQPInt>& cp,
                             std::vector<OSQPInt>& ri,
                             std::vector<OSQPFloat>& vx) {
    int n = Asym.rows;
    cp.assign(n + 1, 0);
    // Count upper-triangular entries per column. Column j entries are rows i ≤ j.
    // Traverse symmetric CSR row i: for each (i, j) with i ≤ j, it lives in column j.
    for (int i = 0; i < n; ++i) {
        int b = Asym.rowStart[i], e = Asym.rowStart[i + 1];
        for (int k = b; k < e; ++k) {
            int j = Asym.colIndex[k];
            if (i <= j) cp[j + 1]++;
        }
    }
    for (int j = 0; j < n; ++j) cp[j + 1] += cp[j];
    OSQPInt nnz = cp[n];
    ri.assign((size_t)nnz, 0);
    vx.assign((size_t)nnz, 0.0);
    std::vector<OSQPInt> cursor(cp.begin(), cp.begin() + n);
    for (int i = 0; i < n; ++i) {
        int b = Asym.rowStart[i], e = Asym.rowStart[i + 1];
        for (int k = b; k < e; ++k) {
            int j = Asym.colIndex[k];
            if (i > j) continue;
            OSQPInt pos = cursor[j]++;
            ri[(size_t)pos] = i;
            vx[(size_t)pos] = Asym.values[k];
        }
    }
    // Rows within each column should be ascending — sort defensively.
    for (int j = 0; j < n; ++j) {
        OSQPInt b = cp[j], e = cp[j + 1];
        // Small columns (typ. <100) — insertion sort is fine.
        for (OSQPInt a = b + 1; a < e; ++a) {
            OSQPInt r = ri[(size_t)a]; double v = vx[(size_t)a];
            OSQPInt q = a - 1;
            while (q >= b && ri[(size_t)q] > r) {
                ri[(size_t)(q+1)] = ri[(size_t)q];
                vx[(size_t)(q+1)] = vx[(size_t)q];
                --q;
            }
            ri[(size_t)(q+1)] = r; vx[(size_t)(q+1)] = v;
        }
    }
}

#endif

SkinData bbwWeights(const MeshData& mesh,
                    const Skeleton& skeleton,
                    const BBWOptions& opts) {
    SkinData out;
    const size_t nVerts = mesh.vertexCount();
    const size_t nBones = skeleton.bones.size();
    const int K = std::max(1, opts.maxInfluences);

    out.boneCount = nBones;
    out.boneWeights.assign(nVerts * K, 0.0f);
    out.boneIndices.assign(nVerts * K, 0u);
    out.inverseBindMatrices.assign(nBones * 16, 0.0f);
    for (size_t i = 0; i < nBones; ++i) {
        std::memcpy(&out.inverseBindMatrices[i * 16],
                    skeleton.bones[i].inverseBind, sizeof(float) * 16);
    }
    if (nVerts == 0 || nBones == 0 || mesh.indices.empty()) return out;

#if !BROMESH_HAS_OSQP
    (void)opts;
    return out;
#else
    const int n = (int)nVerts;

    // Laplacian, mass, biharmonic Q.
    SparseCsr L; std::vector<double> mass;
    assembleCotangentLaplacian(mesh, L, mass);
    SparseCsr Q;
    buildBiharmonic(L, mass, Q);

    // Tiny diagonal regularization so Q is strictly PD (OSQP is robust to
    // semi-definite P but convergence is slower; a small shift is cheap).
    {
        double avgDiag = 0.0;
        for (int i = 0; i < n; ++i) {
            int b = Q.rowStart[i], e = Q.rowStart[i + 1];
            for (int k = b; k < e; ++k) if (Q.colIndex[k] == i) { avgDiag += Q.values[k]; break; }
        }
        avgDiag = (n > 0) ? (avgDiag / n) : 1.0;
        double reg = avgDiag * 1e-6 + 1e-12;
        for (int i = 0; i < n; ++i) {
            int b = Q.rowStart[i], e = Q.rowStart[i + 1];
            for (int k = b; k < e; ++k) if (Q.colIndex[k] == i) { Q.values[k] += reg; break; }
        }
    }

    // Bone segments.
    std::vector<Vec3> heads, tails;
    extractBoneHeads(skeleton, heads, tails);

    // Anchor vertices per bone: the K nearest verts to the bone's segment.
    int anchorsPer = std::max(1, opts.anchorsPerBone);
    std::vector<std::vector<int>> anchorsByBone(nBones);
    for (size_t b = 0; b < nBones; ++b) {
        std::vector<std::pair<double, int>> dv;
        dv.reserve(nVerts);
        for (size_t v = 0; v < nVerts; ++v) {
            Vec3 p = { mesh.positions[v*3+0],
                       mesh.positions[v*3+1],
                       mesh.positions[v*3+2] };
            dv.emplace_back(distPointSegment(p, heads[b], tails[b]), (int)v);
        }
        int take = std::min((int)dv.size(), anchorsPer);
        std::partial_sort(dv.begin(), dv.begin() + take, dv.end(),
                          [](auto& a, auto& c){ return a.first < c.first; });
        anchorsByBone[b].reserve(take);
        for (int i = 0; i < take; ++i) anchorsByBone[b].push_back(dv[i].second);
    }

    // isAnchorOfSomeBone[v] = bone index if v is anchor of some bone (first
    // wins — conflicts are rare with distinct bones), else -1. Used to pin
    // other bones' solves to zero at this vertex.
    std::vector<int> anchorOwner(nVerts, -1);
    for (size_t b = 0; b < nBones; ++b) {
        for (int v : anchorsByBone[b]) {
            if (anchorOwner[v] < 0) anchorOwner[v] = (int)b;
        }
    }

    // OSQP data shared across bones: P (upper tri of Q) and A (identity).
    std::vector<OSQPInt> Pp, Pi;
    std::vector<OSQPFloat> Px;
    symCsrToUpperCsc(Q, Pp, Pi, Px);

    // A = identity in CSC: one entry per column.
    std::vector<OSQPInt> Ap(n + 1, 0), Ai(n, 0);
    std::vector<OSQPFloat> Ax(n, 1.0);
    for (int j = 0; j <= n; ++j) Ap[j] = j;
    for (int j = 0; j < n; ++j) Ai[j] = j;

    std::vector<OSQPFloat> q(n, 0.0);
    std::vector<OSQPFloat> l(n, 0.0), u(n, 1.0);

    OSQPCscMatrix Pm, Am;
    OSQPCscMatrix_set_data(&Pm, n, n, (OSQPInt)Px.size(), Px.data(), Pi.data(), Pp.data());
    OSQPCscMatrix_set_data(&Am, n, n, (OSQPInt)Ax.size(), Ax.data(), Ai.data(), Ap.data());

    OSQPSettings* settings = OSQPSettings_new();
    osqp_set_default_settings(settings);
    settings->verbose = 0;
    settings->eps_abs = opts.eps;
    settings->eps_rel = opts.eps;
    settings->max_iter = opts.maxIter;
    settings->polishing = 1;

    OSQPSolver* solver = nullptr;
    OSQPInt exit_setup = osqp_setup(&solver, &Pm, q.data(), &Am,
                                    l.data(), u.data(), (OSQPInt)n, (OSQPInt)n, settings);
    if (exit_setup != 0 || !solver) {
        if (solver) osqp_cleanup(solver);
        OSQPSettings_free(settings);
        return out;
    }

    std::vector<std::vector<std::pair<int, float>>> perVertBone(nVerts);
    for (auto& pv : perVertBone) pv.reserve((size_t)K + 2);

    for (size_t b = 0; b < nBones; ++b) {
        // Build bound vectors for this bone.
        for (int v = 0; v < n; ++v) { l[v] = 0.0; u[v] = 1.0; }
        for (int v : anchorsByBone[b]) { l[v] = 1.0; u[v] = 1.0; }
        for (size_t ob = 0; ob < nBones; ++ob) {
            if (ob == b) continue;
            for (int v : anchorsByBone[ob]) { l[v] = 0.0; u[v] = 0.0; }
        }
        osqp_update_data_vec(solver, nullptr, l.data(), u.data());
        OSQPInt exit_solve = osqp_solve(solver);
        if (exit_solve != 0 || !solver->solution || !solver->solution->x) continue;

        const OSQPFloat* w = solver->solution->x;
        const double minW = (double)opts.minWeight;
        for (size_t v = 0; v < nVerts; ++v) {
            double wv = w[v];
            if (wv < minW) continue;
            if (wv > 1.0) wv = 1.0;
            perVertBone[v].emplace_back((int)b, (float)wv);
        }
    }

    osqp_cleanup(solver);
    OSQPSettings_free(settings);

    // Top-K + renormalize per vertex.
    for (size_t v = 0; v < nVerts; ++v) {
        auto& list = perVertBone[v];
        std::sort(list.begin(), list.end(),
                  [](const auto& a, const auto& c){ return a.second > c.second; });
        if ((int)list.size() > K) list.resize((size_t)K);
        float sum = 0.0f;
        for (auto& pr : list) sum += pr.second;
        if (sum < 1e-20f) {
            // Fallback: nearest-bone gets full weight.
            double bestD = std::numeric_limits<double>::infinity();
            int bestB = 0;
            Vec3 p = { mesh.positions[v*3+0],
                       mesh.positions[v*3+1],
                       mesh.positions[v*3+2] };
            for (size_t b = 0; b < nBones; ++b) {
                double d = distPointSegment(p, heads[b], tails[b]);
                if (d < bestD) { bestD = d; bestB = (int)b; }
            }
            out.boneWeights[v * K] = 1.0f;
            out.boneIndices[v * K] = (uint32_t)bestB;
            continue;
        }
        float inv = 1.0f / sum;
        for (int k = 0; k < (int)list.size(); ++k) {
            out.boneWeights[v * K + k] = list[k].second * inv;
            out.boneIndices[v * K + k] = (uint32_t)list[k].first;
        }
    }

    return out;
#endif
}

} // namespace bromesh
