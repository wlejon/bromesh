#include "bromesh/rigging/bone_heat.h"
#include "bromesh/rigging/mesh_laplacian.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

namespace bromesh {

namespace {

using Vec3 = std::array<double, 3>;

// --- duplicated from voxel_bind for bone segment extraction ----------------
// Keeping a local copy rather than carving a shared helper: the 4x4 inverse
// is tiny and not worth a header for.
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

// Distance from point p to segment (a, b). Returns 0 for zero-length seg.
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
            // Leaf: short stub along +Y of bind pose.
            double w[16];
            matInverse4(skeleton.bones[i].inverseBind, w);
            Vec3 yAxis = { w[4], w[5], w[6] };
            // Length: use parent distance if available, else unit.
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

} // namespace

SkinData boneHeatWeights(const MeshData& mesh,
                         const Skeleton& skeleton,
                         const BoneHeatOptions& opts) {
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

    // Cotangent Laplacian + mass matrix.
    SparseCsr L;
    std::vector<double> mass;
    assembleCotangentLaplacian(mesh, L, mass);

    // Bone segments.
    std::vector<Vec3> heads, tails;
    extractBoneHeads(skeleton, heads, tails);

    // Per-vertex nearest-bone distance (for heat weight p_j and for the
    // "visible" indicator). We use the simple "nearest bone wins" rule: if
    // bone j is the nearest, the vertex gets a strong heat source tying
    // w_j toward 1. Baran & Popović also introduce a visibility test via
    // raycast; we skip that — for reasonably-fit skeletons, the nearest
    // bone is almost always the right one, and raycasting adds runtime
    // and a BVH dependency. Callers with pathological topology should use
    // voxelBindWeights instead.
    std::vector<int> nearestBone(nVerts, 0);
    std::vector<double> nearestDist(nVerts, 0.0);
    for (size_t v = 0; v < nVerts; ++v) {
        Vec3 p = { mesh.positions[v*3+0],
                   mesh.positions[v*3+1],
                   mesh.positions[v*3+2] };
        double bestD = std::numeric_limits<double>::infinity();
        int    bestB = 0;
        for (size_t b = 0; b < nBones; ++b) {
            double d = distPointSegment(p, heads[b], tails[b]);
            if (d < bestD) { bestD = d; bestB = (int)b; }
        }
        nearestBone[v] = bestB;
        nearestDist[v] = bestD;
    }

    // Reference distance for heat-weight scaling: median nearest-bone distance.
    // Gives the algorithm a scale-invariant baseline.
    double refDist = 0.0;
    {
        std::vector<double> tmp = nearestDist;
        if (!tmp.empty()) {
            std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
            refDist = tmp[tmp.size()/2];
        }
        if (refDist < 1e-8) refDist = 1.0;
    }

    // Assemble the system matrix A = -L + diag(H * mass_i) * mass_i?
    // We use the standard formulation: the equilibrium PDE solved over
    // the mass-weighted Laplacian.
    //
    //   -L w + M H w = M H p
    //
    // where M is the mass matrix (diagonal), L is the cotangent Laplacian
    // (our L has negative diagonal — so -L is SPD), H is the per-vertex
    // heat weight chosen to pin w near p_j at vertex v if bone j is v's
    // nearest bone. Canonical choice: H_v = c / d_v^2 where d_v is the
    // distance from vertex v to its nearest bone.
    //
    // In CSR form we precompute the row pattern of -L, then add the extra
    // diagonal MH term at solve time per bone. Since H doesn't depend on
    // the bone (only on the nearest-bone distance per vertex), the system
    // matrix is SHARED across bones — one CG matrix, many right-hand sides.
    SparseCsr Asys = L;
    for (auto& v : Asys.values) v = -v;
    // Add diagonal term M * H per vertex.
    std::vector<double> heatDiag(nVerts, 0.0);
    double strength = (double)opts.heatStrength;
    for (size_t v = 0; v < nVerts; ++v) {
        double d = nearestDist[v];
        // Avoid blowup near-zero distances; also avoid underflow for far-away
        // verts by anchoring to refDist.
        double dsq = d * d + 1e-8 * refDist * refDist;
        heatDiag[v] = strength * (refDist * refDist) / dsq;
    }
    for (int i = 0; i < Asys.rows; ++i) {
        int bs = Asys.rowStart[i], es = Asys.rowStart[i + 1];
        for (int k = bs; k < es; ++k) {
            if (Asys.colIndex[k] == i) {
                Asys.values[k] += mass[i] * heatDiag[i];
                break;
            }
        }
    }

    // Per-bone solve: w_j with rhs b_j[v] = mass[v] * H[v] * (nearestBone[v]==j ? 1 : 0).
    // Stack onto per-vertex top-K list.
    std::vector<std::vector<std::pair<int, float>>> perVertBone(nVerts);
    for (auto& pv : perVertBone) pv.reserve((size_t)K + 2);

    std::vector<double> rhs(nVerts, 0.0), wSol(nVerts, 0.0);
    for (size_t b = 0; b < nBones; ++b) {
        // Build rhs.
        bool anyAssigned = false;
        for (size_t v = 0; v < nVerts; ++v) {
            if (nearestBone[v] == (int)b) {
                rhs[v] = mass[v] * heatDiag[v];
                anyAssigned = true;
            } else {
                rhs[v] = 0.0;
            }
        }
        if (!anyAssigned) continue;
        std::fill(wSol.begin(), wSol.end(), 0.0);
        sparseCg(Asys, rhs, wSol, opts.solverTol, opts.solverMaxIter);

        // Insert into per-vertex top-K.
        const double minW = (double)opts.minWeight;
        for (size_t v = 0; v < nVerts; ++v) {
            double w = wSol[v];
            if (w < minW) continue;
            // Keep list sorted descending; trim to K.
            auto& list = perVertBone[v];
            list.emplace_back((int)b, (float)w);
        }
    }

    // Finalize per vertex: sort desc, keep top K, renormalize.
    for (size_t v = 0; v < nVerts; ++v) {
        auto& list = perVertBone[v];
        std::sort(list.begin(), list.end(),
                  [](const auto& a, const auto& c){ return a.second > c.second; });
        if ((int)list.size() > K) list.resize((size_t)K);
        float sum = 0.0f;
        for (auto& pr : list) sum += pr.second;
        if (sum < 1e-20f) {
            // Fallback: assign full weight to nearest bone.
            out.boneWeights[v * K] = 1.0f;
            out.boneIndices[v * K] = (uint32_t)nearestBone[v];
            continue;
        }
        float inv = 1.0f / sum;
        for (int k = 0; k < (int)list.size(); ++k) {
            out.boneWeights[v * K + k] = list[k].second * inv;
            out.boneIndices[v * K + k] = (uint32_t)list[k].first;
        }
    }

    return out;
}

} // namespace bromesh
