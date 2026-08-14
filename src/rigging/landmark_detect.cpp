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

    // -- Crown: centroid of highest vertices near midline --------------------
    float uThresh = uMax - 0.02f * H;
    double crownSumR = 0, crownSumU = 0, crownSumF = 0;
    int crownN = 0;
    for (size_t i = 0; i < V; ++i) {
        if (std::fabs(R[i] - rMidSym) > 0.15f * W) continue;
        if (U[i] >= uThresh) {
            crownSumR += R[i]; crownSumU += U[i]; crownSumF += F[i];
            ++crownN;
        }
    }
    float crownR = (crownN > 0) ? (float)(crownSumR / crownN) : rMidSym;
    float crownU = (crownN > 0) ? (float)(crownSumU / crownN) : uMax;
    float crownF = (crownN > 0) ? (float)(crownSumF / crownN) : 0.0f;

    // -- Ankles: centroid of the bottom vertices on each side.
    auto pickAnkle = [&](bool leftSide, float& oR, float& oU, float& oF) {
        float ankThresh = uMin + 0.05f * H;
        double sR = 0, sU = 0, sF = 0;
        int n = 0;
        for (size_t i = 0; i < V; ++i) {
            float dr = R[i] - rMidSym;
            if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
            if (U[i] <= ankThresh) {
                sR += R[i]; sU += U[i]; sF += F[i];
                ++n;
            }
        }
        if (n == 0) {
            // fallback: find single minimum
            float bestU = std::numeric_limits<float>::infinity();
            for (size_t i = 0; i < V; ++i) {
                float dr = R[i] - rMidSym;
                if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
                if (U[i] < bestU) { bestU = U[i]; sR = R[i]; sU = U[i]; sF = F[i]; n = 1; }
            }
        }
        oR = (n > 0) ? (float)(sR / n) : (leftSide ? rMidSym - 0.2f * W : rMidSym + 0.2f * W);
        oU = (n > 0) ? (float)(sU / n) : uMin;
        oF = (n > 0) ? (float)(sF / n) : 0.0f;
    };
    float ankLR, ankLU, ankLF, ankRR, ankRU, ankRF;
    pickAnkle(true,  ankLR, ankLU, ankLF);
    pickAnkle(false, ankRR, ankRU, ankRF);

    // -- Wrists: centroid of outermost vertices in the upper half -------------
    auto pickWrist = [&](bool leftSide, float& oR, float& oU, float& oF) {
        float extremeR = leftSide ? rMin : rMax;
        float rThresh = leftSide ? (extremeR + 0.05f * W) : (extremeR - 0.05f * W);
        double sR = 0, sU = 0, sF = 0;
        int n = 0;
        for (size_t i = 0; i < V; ++i) {
            if (U[i] < uMin + 0.4f * H) continue;
            if (leftSide ? (R[i] <= rThresh) : (R[i] >= rThresh)) {
                sR += R[i]; sU += U[i]; sF += F[i];
                ++n;
            }
        }
        if (n == 0) {
            float bestR = leftSide ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity();
            for (size_t i = 0; i < V; ++i) {
                if (U[i] < uMin + 0.4f * H) continue;
                if (leftSide ? R[i] < bestR : R[i] > bestR) {
                    bestR = R[i]; sR = R[i]; sU = U[i]; sF = F[i]; n = 1;
                }
            }
        }
        oR = (n > 0) ? (float)(sR / n) : (leftSide ? rMin : rMax);
        oU = (n > 0) ? (float)(sU / n) : (uMin + 0.7f * H);
        oF = (n > 0) ? (float)(sF / n) : 0.0f;
    };
    float wrLR, wrLU, wrLF, wrRR, wrRU, wrRF;
    pickWrist(true,  wrLR, wrLU, wrLF);
    pickWrist(false, wrRR, wrRU, wrRF);

    // -- Shoulders: torso half-width measured at the wrist's up-level ---------
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
            best  = 0.25f * (wristR - rMidSym);
            bestF = leftSide ? wrLF : wrRF;
        }
        oR = rMidSym + best;
        oU = shoulderU;
        oF = bestF;
    };
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
    // Chest: upper torso, near shoulder midpoint
    Vec3 chest = { pelvis[0] + 0.90f*(shoulderMid[0]-pelvis[0]),
                   pelvis[1] + 0.90f*(shoulderMid[1]-pelvis[1]),
                   pelvis[2] + 0.90f*(shoulderMid[2]-pelvis[2]) };
    // Neck base: just above shoulder midpoint, toward crown.
    Vec3 neckBase = { shoulderMid[0] + 0.25f*(crown[0]-shoulderMid[0]),
                      shoulderMid[1] + 0.25f*(crown[1]-shoulderMid[1]),
                      shoulderMid[2] + 0.25f*(crown[2]-shoulderMid[2]) };

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
    (void)W;
    const float rMidSym = 0.5f*(rMin + rMax);
    const float fMid    = 0.5f*(fMin + fMax);

    // Muzzle / tail tip: fwd extrema in the upper half of the body.
    float muThresh = fMax - 0.05f * L;
    float ttThresh = fMin + 0.05f * L;
    double muSumR = 0, muSumU = 0, muSumF = 0; int muN = 0;
    double ttSumR = 0, ttSumU = 0, ttSumF = 0; int ttN = 0;
    for (size_t i = 0; i < V; ++i) {
        if (U[i] < uMin + 0.4f * H) continue;
        if (F[i] >= muThresh) { muSumR += R[i]; muSumU += U[i]; muSumF += F[i]; ++muN; }
        if (F[i] <= ttThresh) { ttSumR += R[i]; ttSumU += U[i]; ttSumF += F[i]; ++ttN; }
    }
    float muR = (muN > 0) ? (float)(muSumR / muN) : rMidSym;
    float muU = (muN > 0) ? (float)(muSumU / muN) : (uMin + 0.7f * H);
    float muF = (muN > 0) ? (float)(muSumF / muN) : fMax;
    float ttR = (ttN > 0) ? (float)(ttSumR / ttN) : rMidSym;
    float ttU = (ttN > 0) ? (float)(ttSumU / ttN) : (uMin + 0.7f * H);
    float ttF = (ttN > 0) ? (float)(ttSumF / ttN) : fMin;

    // Crown: centroid of top-5% up vertices in the front half.
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

    // Four paws: centroid of lowest vertices in each quadrant
    auto pickPaw = [&](bool leftSide, bool front,
                       float& oR, float& oU, float& oF) {
        float threshU = uMin + 0.05f * H;
        double sR = 0, sU = 0, sF = 0;
        int n = 0;
        for (size_t i = 0; i < V; ++i) {
            float dr = R[i] - rMidSym;
            if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
            float df = F[i] - fMid;
            if (front ? !(df > 0) : !(df < 0)) continue;
            if (U[i] <= threshU) {
                sR += R[i]; sU += U[i]; sF += F[i];
                ++n;
            }
        }
        if (n == 0) {
            float bestU = std::numeric_limits<float>::infinity();
            for (size_t i = 0; i < V; ++i) {
                if (U[i] > uMin + 0.3f*H) continue;
                float dr = R[i] - rMidSym;
                if (leftSide ? !(dr < 0) : !(dr > 0)) continue;
                float df = F[i] - fMid;
                if (front ? !(df > 0) : !(df < 0)) continue;
                if (U[i] < bestU) { bestU=U[i]; sR=R[i]; sU=U[i]; sF=F[i]; n = 1; }
            }
        }
        oR = (n > 0) ? (float)(sR / n) : (leftSide ? rMidSym - 0.2f * W : rMidSym + 0.2f * W);
        oU = (n > 0) ? (float)(sU / n) : uMin;
        oF = (n > 0) ? (float)(sF / n) : (front ? fMid + 0.2f * L : fMid - 0.2f * L);
    };
    float fpLR, fpLU, fpLF, fpRR, fpRU, fpRF;
    float hpLR, hpLU, hpLF, hpRR, hpRU, hpRF;
    pickPaw(true,  true,  fpLR, fpLU, fpLF);
    pickPaw(false, true,  fpRR, fpRU, fpRF);
    pickPaw(true,  false, hpLR, hpLU, hpLF);
    pickPaw(false, false, hpRR, hpRU, hpRF);

    // Shoulder / hip U: ~75% of body height.
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

    // Tail base: 25% past the pelvis along the chest→pelvis direction.
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
