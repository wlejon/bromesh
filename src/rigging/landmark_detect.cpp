#include "bromesh/rigging/landmark_detect.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace bromesh {

namespace {

using Vec3 = std::array<float, 3>;

Vec3 normalize(const float v[3]) {
    float l = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (l < 1e-12f) return { 0, 1, 0 };
    return { v[0]/l, v[1]/l, v[2]/l };
}

Vec3 cross(Vec3 a, Vec3 b) {
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}

// World-space point at local (right, up, fwd) coordinates in the chosen frame.
Vec3 unproject(Vec3 right, Vec3 up, Vec3 fwd, float r, float u, float f) {
    return { right[0]*r + up[0]*u + fwd[0]*f,
             right[1]*r + up[1]*u + fwd[1]*f,
             right[2]*r + up[2]*u + fwd[2]*f };
}

void setLm(Landmarks& lm, const char* name, Vec3 p) {
    lm.set(name, p[0], p[1], p[2]);
}

} // namespace

Landmarks detectHumanoidLandmarks(const MeshData& mesh,
                                  const LandmarkDetectOptions& opts) {
    Landmarks out;
    if (mesh.vertexCount() == 0) return out;

    // Orthonormal frame: up, forward, right = cross(up, forward).
    Vec3 up  = normalize(opts.upAxis);
    Vec3 fwd = normalize(opts.forwardAxis);
    // Re-orthogonalise forward against up.
    float upDotFwd = up[0]*fwd[0] + up[1]*fwd[1] + up[2]*fwd[2];
    fwd = { fwd[0] - up[0]*upDotFwd,
            fwd[1] - up[1]*upDotFwd,
            fwd[2] - up[2]*upDotFwd };
    fwd = normalize(fwd.data());
    Vec3 right = cross(up, fwd);

    const size_t V = mesh.vertexCount();
    std::vector<float> R(V), U(V), F(V);
    for (size_t i = 0; i < V; ++i) {
        float x = mesh.positions[i*3+0];
        float y = mesh.positions[i*3+1];
        float z = mesh.positions[i*3+2];
        R[i] = x*right[0] + y*right[1] + z*right[2];
        U[i] = x*up[0]    + y*up[1]    + z*up[2];
        F[i] = x*fwd[0]   + y*fwd[1]   + z*fwd[2];
    }

    float uMin = U[0], uMax = U[0];
    float rMin = R[0], rMax = R[0];
    for (size_t i = 1; i < V; ++i) {
        uMin = std::min(uMin, U[i]); uMax = std::max(uMax, U[i]);
        rMin = std::min(rMin, R[i]); rMax = std::max(rMax, R[i]);
    }
    const float H = uMax - uMin;                    // body height
    const float W = std::max(rMax - rMin, 1e-6f);   // arm-span
    const float rMidSym = 0.5f * (rMin + rMax);     // symmetry plane

    // -- Crown: highest-up vertex near the midline (filter out arms) ----------
    float crownU = -std::numeric_limits<float>::infinity();
    float crownR = 0.0f, crownF = 0.0f;
    for (size_t i = 0; i < V; ++i) {
        if (std::fabs(R[i] - rMidSym) > 0.15f * W) continue;
        if (U[i] > crownU) { crownU = U[i]; crownR = R[i]; crownF = F[i]; }
    }

    // -- Ankles: centroid of the bottom 5% of vertices on each side. Using a
    // single min-U vertex made ankle.F land on whichever toe-tip or heel was
    // lowest, which in turn placed the hip column at the front/back surface
    // of the body (not its centerline). The centroid is stable against foot
    // shape — toes/heels cancel out, leaving the foot's geometric center.
    auto pickAnkle = [&](bool leftSide, float& oR, float& oU, float& oF) {
        std::vector<size_t> cand;
        cand.reserve(V / 4);
        for (size_t i = 0; i < V; ++i) {
            float dr = R[i] - rMidSym;
            if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
            if (U[i] > uMin + 0.25f * H) continue;
            cand.push_back(i);
        }
        if (cand.empty()) { oR = 0; oU = 0; oF = 0; return; }
        std::sort(cand.begin(), cand.end(),
            [&](size_t a, size_t b){ return U[a] < U[b]; });
        size_t n = std::max((size_t)1, cand.size() / 20);
        double sR = 0, sU = 0, sF = 0;
        for (size_t i = 0; i < n; ++i) {
            sR += R[cand[i]]; sU += U[cand[i]]; sF += F[cand[i]];
        }
        oR = (float)(sR / (double)n);
        oU = (float)(sU / (double)n);
        oF = (float)(sF / (double)n);
    };
    float ankLR, ankLU, ankLF, ankRR, ankRU, ankRF;
    pickAnkle(true,  ankLR, ankLU, ankLF);
    pickAnkle(false, ankRR, ankRU, ankRF);

    // -- Wrists: extreme side in the upper half --------------------------------
    auto pickWrist = [&](bool leftSide, float& oR, float& oU, float& oF) {
        float bestR = leftSide ? std::numeric_limits<float>::infinity()
                               : -std::numeric_limits<float>::infinity();
        oR = 0; oU = 0; oF = 0;
        for (size_t i = 0; i < V; ++i) {
            if (U[i] < uMin + 0.4f * H) continue; // upper body only
            if (leftSide ? R[i] < bestR : R[i] > bestR) {
                bestR = R[i]; oR = R[i]; oU = U[i]; oF = F[i];
            }
        }
    };
    float wrLR, wrLU, wrLF, wrRR, wrRU, wrRF;
    pickWrist(true,  wrLR, wrLU, wrLF);
    pickWrist(false, wrRR, wrRU, wrRF);

    // -- Shoulders: torso half-width measured at the wrist's up-level ---------
    // Scan a thin band at shoulder height, collect vertices that are clearly
    // not arm-extrema (|r - midSym| < half the wrist offset), take the widest.
    auto shoulderAt = [&](bool leftSide, float shoulderU,
                          float& oR, float& oU, float& oF) {
        float band = 0.05f * H;
        float wristR = leftSide ? wrLR : wrRR;
        float limit  = 0.5f * std::fabs(wristR - rMidSym);
        float best = 0.0f;
        float bestF = 0.0f;
        float count = 0.0f;
        for (size_t i = 0; i < V; ++i) {
            if (std::fabs(U[i] - shoulderU) > band) continue;
            float dr = R[i] - rMidSym;
            if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
            if (std::fabs(dr) > limit) continue; // skip arm vertices
            if (std::fabs(dr) > std::fabs(best)) { best = dr; bestF = F[i]; }
            count += 1.0f;
        }
        if (count < 1.0f) {
            // Fallback: 25% of the way from midline to wrist.
            best  = 0.25f * (wristR - rMidSym);
            bestF = leftSide ? wrLF : wrRF;
        }
        oR = rMidSym + best;
        oU = shoulderU;
        oF = bestF;
    };
    // T-pose shoulder height ≈ wrist up-coord (arms horizontal).
    float shoulderU = 0.5f * (wrLU + wrRU);
    float shLR, shLU, shLF, shRR, shRU, shRF;
    shoulderAt(true,  shoulderU, shLR, shLU, shLF);
    shoulderAt(false, shoulderU, shRR, shRU, shRF);

    // -- Hips: ~half-way up the body, at the leg columns' horizontal position.
    // Legs in a T-pose hang straight down, so hip_R ≈ ankle_R.
    float hipU = uMin + 0.5f * H;
    float hipLR = ankLR, hipLF = ankLF;
    float hipRR = ankRR, hipRF = ankRF;

    // -- Pelvis, chest, neck_base, elbows, knees, toes, crown -----------------
    Vec3 crown    = unproject(right, up, fwd, crownR, crownU, crownF);
    Vec3 ankleL   = unproject(right, up, fwd, ankLR, ankLU, ankLF);
    Vec3 ankleR   = unproject(right, up, fwd, ankRR, ankRU, ankRF);
    Vec3 wristL   = unproject(right, up, fwd, wrLR,  wrLU,  wrLF);
    Vec3 wristR   = unproject(right, up, fwd, wrRR,  wrRU,  wrRF);
    Vec3 shoulderL= unproject(right, up, fwd, shLR,  shLU,  shLF);
    Vec3 shoulderR= unproject(right, up, fwd, shRR,  shRU,  shRF);
    Vec3 hipL     = unproject(right, up, fwd, hipLR, hipU,  hipLF);
    Vec3 hipR     = unproject(right, up, fwd, hipRR, hipU,  hipRF);

    Vec3 pelvis = { 0.5f*(hipL[0]+hipR[0]),
                    0.5f*(hipL[1]+hipR[1]),
                    0.5f*(hipL[2]+hipR[2]) };
    Vec3 shoulderMid = { 0.5f*(shoulderL[0]+shoulderR[0]),
                         0.5f*(shoulderL[1]+shoulderR[1]),
                         0.5f*(shoulderL[2]+shoulderR[2]) };
    // Chest: 60% up from pelvis to shoulders.
    Vec3 chest = { pelvis[0] + 0.60f*(shoulderMid[0]-pelvis[0]),
                   pelvis[1] + 0.60f*(shoulderMid[1]-pelvis[1]),
                   pelvis[2] + 0.60f*(shoulderMid[2]-pelvis[2]) };
    // Neck base: just above shoulder midpoint, toward crown.
    Vec3 neckBase = { shoulderMid[0] + 0.15f*(crown[0]-shoulderMid[0]),
                      shoulderMid[1] + 0.15f*(crown[1]-shoulderMid[1]),
                      shoulderMid[2] + 0.15f*(crown[2]-shoulderMid[2]) };

    auto lerp = [](Vec3 a, Vec3 b, float t) {
        return Vec3{ a[0]+t*(b[0]-a[0]), a[1]+t*(b[1]-a[1]), a[2]+t*(b[2]-a[2]) };
    };
    Vec3 elbowL = lerp(shoulderL, wristL, 0.5f);
    Vec3 elbowR = lerp(shoulderR, wristR, 0.5f);
    Vec3 kneeL  = lerp(hipL, ankleL, 0.5f);
    Vec3 kneeR  = lerp(hipR, ankleR, 0.5f);
    float footLen = opts.footLengthFrac * H;
    Vec3 toeL = { ankleL[0]+fwd[0]*footLen,
                  ankleL[1]+fwd[1]*footLen,
                  ankleL[2]+fwd[2]*footLen };
    Vec3 toeR = { ankleR[0]+fwd[0]*footLen,
                  ankleR[1]+fwd[1]*footLen,
                  ankleR[2]+fwd[2]*footLen };

    setLm(out, "pelvis",     pelvis);
    setLm(out, "chest",      chest);
    setLm(out, "neck_base",  neckBase);
    setLm(out, "crown",      crown);
    setLm(out, "shoulder_L", shoulderL);
    setLm(out, "shoulder_R", shoulderR);
    setLm(out, "elbow_L",    elbowL);
    setLm(out, "elbow_R",    elbowR);
    setLm(out, "wrist_L",    wristL);
    setLm(out, "wrist_R",    wristR);
    setLm(out, "hip_L",      hipL);
    setLm(out, "hip_R",      hipR);
    setLm(out, "knee_L",     kneeL);
    setLm(out, "knee_R",     kneeR);
    setLm(out, "ankle_L",    ankleL);
    setLm(out, "ankle_R",    ankleR);
    setLm(out, "toe_L",      toeL);
    setLm(out, "toe_R",      toeR);
    return out;
}

Landmarks detectQuadrupedLandmarks(const MeshData& mesh,
                                   const LandmarkDetectOptions& opts) {
    Landmarks out;
    if (mesh.vertexCount() == 0) return out;

    Vec3 up  = normalize(opts.upAxis);
    Vec3 fwd = normalize(opts.forwardAxis);
    float upDotFwd = up[0]*fwd[0] + up[1]*fwd[1] + up[2]*fwd[2];
    fwd = { fwd[0] - up[0]*upDotFwd,
            fwd[1] - up[1]*upDotFwd,
            fwd[2] - up[2]*upDotFwd };
    fwd = normalize(fwd.data());
    Vec3 right = cross(up, fwd);

    const size_t V = mesh.vertexCount();
    std::vector<float> R(V), U(V), F(V);
    for (size_t i = 0; i < V; ++i) {
        float x = mesh.positions[i*3+0];
        float y = mesh.positions[i*3+1];
        float z = mesh.positions[i*3+2];
        R[i] = x*right[0] + y*right[1] + z*right[2];
        U[i] = x*up[0]    + y*up[1]    + z*up[2];
        F[i] = x*fwd[0]   + y*fwd[1]   + z*fwd[2];
    }

    float uMin = U[0], uMax = U[0];
    float rMin = R[0], rMax = R[0];
    float fMin = F[0], fMax = F[0];
    for (size_t i = 1; i < V; ++i) {
        uMin = std::min(uMin, U[i]); uMax = std::max(uMax, U[i]);
        rMin = std::min(rMin, R[i]); rMax = std::max(rMax, R[i]);
        fMin = std::min(fMin, F[i]); fMax = std::max(fMax, F[i]);
    }
    const float H = uMax - uMin;
    const float L = std::max(fMax - fMin, 1e-6f);
    const float W = std::max(rMax - rMin, 1e-6f);
    (void)L; (void)W;
    const float rMidSym = 0.5f*(rMin + rMax);
    const float fMid    = 0.5f*(fMin + fMax);

    // Muzzle / tail tip: fwd extrema in the upper half of the body.
    float muR=0, muU=0, muF=-std::numeric_limits<float>::infinity();
    float ttR=0, ttU=0, ttF= std::numeric_limits<float>::infinity();
    for (size_t i = 0; i < V; ++i) {
        if (U[i] < uMin + 0.4f*H) continue;
        if (F[i] > muF) { muF=F[i]; muR=R[i]; muU=U[i]; }
        if (F[i] < ttF) { ttF=F[i]; ttR=R[i]; ttU=U[i]; }
    }

    // Crown: centroid of top-5% up vertices in the front half. More stable
    // than argmax on discrete meshes where the "top" is a whole face.
    float uThresh = uMax - 0.05f*H;
    float crR = rMidSym, crU = uMax, crF = muF;
    {
        double sumR=0, sumU=0, sumF=0; int n=0;
        for (size_t i = 0; i < V; ++i) {
            if (U[i] < uThresh) continue;
            if (F[i] < fMid)    continue;
            sumR += R[i]; sumU += U[i]; sumF += F[i]; ++n;
        }
        if (n > 0) {
            crR = float(sumR / n);
            crU = float(sumU / n);
            crF = float(sumF / n);
        }
    }

    // Four paws: argmin-up in each (left/right × front/back) quadrant, filtered
    // to the lower third of the body.
    auto pickPaw = [&](bool leftSide, bool front,
                       float& oR, float& oU, float& oF) {
        float bestU = std::numeric_limits<float>::infinity();
        oR = 0; oU = 0; oF = 0;
        for (size_t i = 0; i < V; ++i) {
            if (U[i] > uMin + 0.3f*H) continue;
            float dr = R[i] - rMidSym;
            if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
            float df = F[i] - fMid;
            if (front ? !(df > 0) : !(df < 0)) continue;
            if (U[i] < bestU) { bestU=U[i]; oR=R[i]; oU=U[i]; oF=F[i]; }
        }
    };
    float fpLR, fpLU, fpLF, fpRR, fpRU, fpRF;
    float hpLR, hpLU, hpLF, hpRR, hpRU, hpRF;
    pickPaw(true,  true,  fpLR, fpLU, fpLF);
    pickPaw(false, true,  fpRR, fpRU, fpRF);
    pickPaw(true,  false, hpLR, hpLU, hpLF);
    pickPaw(false, false, hpRR, hpRU, hpRF);

    // Shoulder / hip U: ~75% of body height. In a natural stance the legs
    // hang vertically, so the leg column's R,F carries straight up to the
    // torso attachment; we don't need to hunt for it in the mesh.
    const float shoulderHipU = uMin + 0.75f*H;

    Vec3 muzzle    = unproject(right, up, fwd, muR,  muU,  muF);
    Vec3 tailTip   = unproject(right, up, fwd, ttR,  ttU,  ttF);
    Vec3 crown     = unproject(right, up, fwd, crR,  crU,  crF);
    Vec3 fpawL     = unproject(right, up, fwd, fpLR, fpLU, fpLF);
    Vec3 fpawR     = unproject(right, up, fwd, fpRR, fpRU, fpRF);
    Vec3 hpawL     = unproject(right, up, fwd, hpLR, hpLU, hpLF);
    Vec3 hpawR     = unproject(right, up, fwd, hpRR, hpRU, hpRF);
    Vec3 shoulderL = unproject(right, up, fwd, fpLR, shoulderHipU, fpLF);
    Vec3 shoulderR = unproject(right, up, fwd, fpRR, shoulderHipU, fpRF);
    Vec3 hipL      = unproject(right, up, fwd, hpLR, shoulderHipU, hpLF);
    Vec3 hipR      = unproject(right, up, fwd, hpRR, shoulderHipU, hpRF);

    auto mid  = [](Vec3 a, Vec3 b) {
        return Vec3{0.5f*(a[0]+b[0]), 0.5f*(a[1]+b[1]), 0.5f*(a[2]+b[2])};
    };
    auto lerp = [](Vec3 a, Vec3 b, float t) {
        return Vec3{a[0]+t*(b[0]-a[0]), a[1]+t*(b[1]-a[1]), a[2]+t*(b[2]-a[2])};
    };

    Vec3 chest    = mid(shoulderL, shoulderR);
    Vec3 pelvis   = mid(hipL, hipR);
    Vec3 neckBase = { chest[0] + 0.20f*(crown[0]-chest[0]),
                      chest[1] + 0.30f*(crown[1]-chest[1]),
                      chest[2] + 0.40f*(crown[2]-chest[2]) };

    Vec3 felbowL = lerp(shoulderL, fpawL, 0.5f);
    Vec3 felbowR = lerp(shoulderR, fpawR, 0.5f);
    Vec3 hkneeL  = lerp(hipL, hpawL, 0.5f);
    Vec3 hkneeR  = lerp(hipR, hpawR, 0.5f);

    // Tail base: 25% past the pelvis along the chest→pelvis direction. This
    // reliably lands at the rump without needing a torso/tail separation
    // heuristic that might trip on messy topology.
    Vec3 tailBase = { pelvis[0] + 0.25f*(pelvis[0]-chest[0]),
                      pelvis[1] + 0.25f*(pelvis[1]-chest[1]),
                      pelvis[2] + 0.25f*(pelvis[2]-chest[2]) };

    setLm(out, "pelvis",      pelvis);
    setLm(out, "chest",       chest);
    setLm(out, "neck_base",   neckBase);
    setLm(out, "crown",       crown);
    setLm(out, "muzzle",      muzzle);
    setLm(out, "tail_base",   tailBase);
    setLm(out, "tail_tip",    tailTip);
    setLm(out, "fshoulder_L", shoulderL);
    setLm(out, "fshoulder_R", shoulderR);
    setLm(out, "felbow_L",    felbowL);
    setLm(out, "felbow_R",    felbowR);
    setLm(out, "fpaw_L",      fpawL);
    setLm(out, "fpaw_R",      fpawR);
    setLm(out, "hip_L",       hipL);
    setLm(out, "hip_R",       hipR);
    setLm(out, "hknee_L",     hkneeL);
    setLm(out, "hknee_R",     hkneeR);
    setLm(out, "hpaw_L",      hpawL);
    setLm(out, "hpaw_R",      hpawR);
    return out;
}

} // namespace bromesh
