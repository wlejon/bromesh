#include "bromesh/procedural/leaf_cluster.h"

#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/normals.h"
#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace bromesh {

using namespace bromath;

namespace {

// Append a sub-mesh to destination MeshData, remapping triangle indices.
void appendMesh(MeshData& dest, const MeshData& src) {
    if (src.empty()) return;
    const uint32_t baseIndex = static_cast<uint32_t>(dest.vertexCount());
    const size_t srcVCount = src.vertexCount();
    const bool hasNormals = src.hasNormals();
    const bool hasUVs = src.hasUVs();
    const bool hasColors = src.hasColors();

    dest.positions.insert(dest.positions.end(), src.positions.begin(), src.positions.end());

    if (hasNormals) {
        dest.normals.insert(dest.normals.end(), src.normals.begin(), src.normals.end());
    } else {
        // Default outward +Y normals if source lacked them
        dest.normals.reserve(dest.normals.size() + srcVCount * 3);
        for (size_t i = 0; i < srcVCount; ++i) {
            dest.normals.push_back(0.0f);
            dest.normals.push_back(1.0f);
            dest.normals.push_back(0.0f);
        }
    }

    if (hasUVs) {
        dest.uvs.insert(dest.uvs.end(), src.uvs.begin(), src.uvs.end());
    } else {
        dest.uvs.reserve(dest.uvs.size() + srcVCount * 2);
        for (size_t i = 0; i < srcVCount; ++i) {
            dest.uvs.push_back(0.0f);
            dest.uvs.push_back(0.0f);
        }
    }

    if (hasColors) {
        dest.colors.insert(dest.colors.end(), src.colors.begin(), src.colors.end());
    } else {
        dest.colors.reserve(dest.colors.size() + srcVCount * 4);
        for (size_t i = 0; i < srcVCount; ++i) {
            dest.colors.push_back(0.0f);
            dest.colors.push_back(0.0f);
            dest.colors.push_back(0.0f);
            dest.colors.push_back(1.0f);
        }
    }

    dest.indices.reserve(dest.indices.size() + src.indices.size());
    for (uint32_t idx : src.indices) {
        dest.indices.push_back(idx + baseIndex);
    }
}

// Build a tapered cylinder / tube along local +Z from (0,0,0) to (0,0,length).
MeshData buildMicroTwig(float length, float radius, int sides = 6, int rings = 4) {
    MeshData m;
    if (length <= 1e-5f || radius <= 1e-6f) return m;

    sides = std::max(3, sides);
    rings = std::max(2, rings);

    const int vcount = rings * (sides + 1);
    m.positions.reserve(vcount * 3);
    m.normals.reserve(vcount * 3);
    m.uvs.reserve(vcount * 2);
    m.colors.reserve(vcount * 4);

    for (int j = 0; j < rings; ++j) {
        float t = static_cast<float>(j) / static_cast<float>(rings - 1);
        float z = t * length;
        // Natural taper: base is full radius, tip narrows to 50%
        float r = radius * (1.0f - 0.5f * t);
        // Wind bend: 0.0 at base attachment, up to 0.35 at twig tip
        float wb = t * 0.35f;

        for (int i = 0; i <= sides; ++i) {
            float u = static_cast<float>(i) / static_cast<float>(sides);
            float theta = u * TWO_PI;
            float cosT = std::cos(theta);
            float sinT = std::sin(theta);

            m.positions.push_back(r * cosT);
            m.positions.push_back(r * sinT);
            m.positions.push_back(z);

            m.normals.push_back(cosT);
            m.normals.push_back(sinT);
            m.normals.push_back(0.0f);

            m.uvs.push_back(u);
            m.uvs.push_back(t);

            m.colors.push_back(wb);
            m.colors.push_back(0.0f);
            m.colors.push_back(0.0f);
            m.colors.push_back(1.0f);
        }
    }

    const int rowStride = sides + 1;
    for (int j = 0; j < rings - 1; ++j) {
        for (int i = 0; i < sides; ++i) {
            uint32_t a = static_cast<uint32_t>(j * rowStride + i);
            uint32_t b = static_cast<uint32_t>(j * rowStride + i + 1);
            uint32_t c = static_cast<uint32_t>((j + 1) * rowStride + i);
            uint32_t d = static_cast<uint32_t>((j + 1) * rowStride + i + 1);

            m.indices.push_back(a);
            m.indices.push_back(c);
            m.indices.push_back(b);

            m.indices.push_back(b);
            m.indices.push_back(c);
            m.indices.push_back(d);
        }
    }

    // End caps: bottom cap at z=0, top cap at z=length
    // Base cap (facing -Z)
    uint32_t baseCenterIdx = static_cast<uint32_t>(m.vertexCount());
    m.positions.push_back(0.0f);
    m.positions.push_back(0.0f);
    m.positions.push_back(0.0f);
    m.normals.push_back(0.0f);
    m.normals.push_back(0.0f);
    m.normals.push_back(-1.0f);
    m.uvs.push_back(0.5f);
    m.uvs.push_back(0.0f);
    m.colors.push_back(0.0f);
    m.colors.push_back(0.0f);
    m.colors.push_back(0.0f);
    m.colors.push_back(1.0f);

    for (int i = 0; i < sides; ++i) {
        uint32_t a = static_cast<uint32_t>(i);
        uint32_t b = static_cast<uint32_t>(i + 1);
        m.indices.push_back(baseCenterIdx);
        m.indices.push_back(a);
        m.indices.push_back(b);
    }

    // Tip cap (facing +Z)
    uint32_t tipCenterIdx = static_cast<uint32_t>(m.vertexCount());
    m.positions.push_back(0.0f);
    m.positions.push_back(0.0f);
    m.positions.push_back(length);
    m.normals.push_back(0.0f);
    m.normals.push_back(0.0f);
    m.normals.push_back(1.0f);
    m.uvs.push_back(0.5f);
    m.uvs.push_back(1.0f);
    m.colors.push_back(0.35f);
    m.colors.push_back(0.0f);
    m.colors.push_back(0.0f);
    m.colors.push_back(1.0f);

    int topRowBase = (rings - 1) * rowStride;
    for (int i = 0; i < sides; ++i) {
        uint32_t a = static_cast<uint32_t>(topRowBase + i);
        uint32_t b = static_cast<uint32_t>(topRowBase + i + 1);
        m.indices.push_back(tipCenterIdx);
        m.indices.push_back(b);
        m.indices.push_back(a);
    }

    return m;
}

// Pick any unit vector perpendicular to t (assumed unit length).
Vec3 perpendicularUnit(Vec3 t) {
    Vec3 worldUp{0, 1, 0};
    Vec3 c = vcross(t, worldUp);
    if (vdot(c, c) < 1e-8f) {
        c = vcross(t, Vec3{1, 0, 0});
    }
    return vnorm(c);
}

// Build a small tapered petiole (leaf stalk) connecting twig attachment to leaf card base.
MeshData buildPetiole(Vec3 from, Vec3 to, float radius, float baseBend, float tipBend, int sides = 4) {
    MeshData m;
    Vec3 axis = to - from;
    float len = vlen(axis);
    if (len <= 1e-5f || radius <= 1e-6f) return m;

    Vec3 dir = axis * (1.0f / len);
    Vec3 uAxis = perpendicularUnit(dir);
    Vec3 vAxis = vnorm(vcross(dir, uAxis));

    sides = std::max(3, sides);
    const int rings = 2;
    const int vcount = rings * (sides + 1);
    m.positions.reserve(vcount * 3);
    m.normals.reserve(vcount * 3);
    m.uvs.reserve(vcount * 2);
    m.colors.reserve(vcount * 4);

    for (int j = 0; j < rings; ++j) {
        float tj = static_cast<float>(j) / static_cast<float>(rings - 1);
        Vec3 center = from + axis * tj;
        float r = radius * (1.0f - 0.25f * tj);
        float wb = baseBend + tj * (tipBend - baseBend);

        for (int i = 0; i <= sides; ++i) {
            float u = static_cast<float>(i) / static_cast<float>(sides);
            float theta = u * TWO_PI;
            float cosT = std::cos(theta);
            float sinT = std::sin(theta);
            Vec3 radial = uAxis * cosT + vAxis * sinT;
            Vec3 p = center + radial * r;

            m.positions.push_back(p.x);
            m.positions.push_back(p.y);
            m.positions.push_back(p.z);

            m.normals.push_back(radial.x);
            m.normals.push_back(radial.y);
            m.normals.push_back(radial.z);

            m.uvs.push_back(u);
            m.uvs.push_back(tj);

            m.colors.push_back(wb);
            m.colors.push_back(0.0f);
            m.colors.push_back(0.0f);
            m.colors.push_back(1.0f);
        }
    }

    const int rowStride = sides + 1;
    for (int i = 0; i < sides; ++i) {
        uint32_t a = static_cast<uint32_t>(i);
        uint32_t b = static_cast<uint32_t>(i + 1);
        uint32_t c = static_cast<uint32_t>(rowStride + i);
        uint32_t d = static_cast<uint32_t>(rowStride + i + 1);

        m.indices.push_back(a);
        m.indices.push_back(c);
        m.indices.push_back(b);

        m.indices.push_back(b);
        m.indices.push_back(c);
        m.indices.push_back(d);
    }

    return m;
}

struct LeafSlot {
    float t;           // Fractional position along twig [0, 1]
    float phi;         // Azimuth around twig axis
    float spreadAngle; // Outward fan angle from +Z
};

std::vector<LeafSlot> computeSlots(Phyllotaxy phyllotaxy, int count, float spread) {
    std::vector<LeafSlot> slots;
    if (count <= 0) return slots;
    slots.reserve(static_cast<size_t>(count));

    switch (phyllotaxy) {
    case Phyllotaxy::Alternate: {
        // Distichous / alternating sides (0, PI, 0, PI) along twig
        for (int k = 0; k < count; ++k) {
            float t = (count == 1) ? 0.5f : (0.15f + 0.80f * (static_cast<float>(k) / static_cast<float>(count - 1)));
            float phi = (k % 2 == 0) ? 0.0f : PI;
            slots.push_back({t, phi, spread});
        }
        break;
    }
    case Phyllotaxy::Opposite: {
        // Decussate pairs (opposite each other, rotated 90° per node)
        int numPairs = (count + 1) / 2;
        for (int p = 0; p < numPairs; ++p) {
            float t = (numPairs == 1) ? 0.5f : (0.15f + 0.80f * (static_cast<float>(p) / static_cast<float>(numPairs - 1)));
            float phiBase = static_cast<float>(p) * (PI * 0.5f);
            slots.push_back({t, phiBase, spread});
            if (static_cast<int>(slots.size()) < count) {
                slots.push_back({t, phiBase + PI, spread});
            }
        }
        break;
    }
    case Phyllotaxy::Spiral: {
        // Golden angle rosette (~137.507764°)
        constexpr float kGoldenAngle = 2.39996322972865332f;
        for (int k = 0; k < count; ++k) {
            float t = (count == 1) ? 0.5f : (0.10f + 0.85f * (static_cast<float>(k) / static_cast<float>(count - 1)));
            float phi = static_cast<float>(k) * kGoldenAngle;
            slots.push_back({t, phi, spread});
        }
        break;
    }
    case Phyllotaxy::Fascicle: {
        // Pine needle bundle radiating from basal sheath
        for (int k = 0; k < count; ++k) {
            float phi = static_cast<float>(k) * (TWO_PI / static_cast<float>(count));
            slots.push_back({0.0f, phi, spread});
        }
        break;
    }
    case Phyllotaxy::CompoundPinnate: {
        // Central rachis: lateral pairs + 1 terminal leaflet at tip
        if (count == 1) {
            slots.push_back({1.0f, 0.0f, 0.0f});
        } else {
            int lateralCount = count - 1;
            int numPairs = (lateralCount + 1) / 2;
            for (int p = 0; p < numPairs; ++p) {
                float t = (numPairs == 1) ? 0.45f : (0.15f + 0.75f * (static_cast<float>(p) / static_cast<float>(numPairs - 1)));
                // Left leaflet (+X side)
                slots.push_back({t, PI * 0.5f, spread});
                // Right leaflet (-X side)
                if (static_cast<int>(slots.size()) < count - 1) {
                    slots.push_back({t, -PI * 0.5f, spread});
                }
            }
            // Terminal leaflet at tip
            slots.push_back({1.0f, 0.0f, 0.0f});
        }
        break;
    }
    }

    return slots;
}

} // namespace

MeshData leafCluster(Phyllotaxy phyllotaxy, const LeafClusterOptions& opts) {
    MeshData out;
    if (opts.count <= 0 && !opts.includeTwigMesh) {
        return out;
    }

    // 1. Build micro-twig / rachis mesh if enabled
    if (opts.includeTwigMesh && opts.twigLength > 1e-5f && opts.twigRadius > 1e-6f) {
        MeshData twig = buildMicroTwig(opts.twigLength, opts.twigRadius);
        appendMesh(out, twig);
    }

    if (opts.count <= 0) {
        return out;
    }

    // 2. Base template leaf card
    LeafCardOptions cardOpts;
    cardOpts.width = opts.leafWidth;
    cardOpts.length = opts.leafLength;
    cardOpts.bend = opts.leafBend;
    cardOpts.curl = opts.leafCurl;
    cardOpts.cup = opts.leafCup;
    cardOpts.fullUV = opts.fullUV;
    cardOpts.shapedSilhouette = opts.shapedSilhouette;
    cardOpts.stemOffset = true;
    cardOpts.widthSegments = 4;
    cardOpts.lengthSegments = 8;

    MeshData baseCard = leafCard(opts.shape, cardOpts);
    if (baseCard.empty()) return out;

    const size_t cardVCount = baseCard.vertexCount();
    const size_t cardICount = baseCard.indices.size();
    const bool cardHasNormals = baseCard.hasNormals();
    const bool cardHasUVs = baseCard.hasUVs();
    const bool cardHasColors = baseCard.hasColors();

    std::vector<LeafSlot> slots = computeSlots(phyllotaxy, opts.count, opts.spread);
    const Vec3 worldUp{0.0f, 1.0f, 0.0f};

    for (const LeafSlot& slot : slots) {
        float t = slot.t;
        Vec3 attachPoint{0.0f, 0.0f, t * opts.twigLength};

        // Outward fan direction from twig axis
        Vec3 dir;
        if (slot.spreadAngle <= 1e-4f) {
            // Terminal leaflet pointing along +Z
            dir = Vec3{0.0f, 0.0f, 1.0f};
            if (opts.droop > 0.0f) {
                dir = vnorm(dir * (1.0f - opts.droop * 0.4f) + Vec3{0.0f, -1.0f, 0.0f} * (opts.droop * 0.4f));
            }
        } else {
            Vec3 radial{std::cos(slot.phi), std::sin(slot.phi), 0.0f};
            dir = Vec3{radial.x * std::sin(slot.spreadAngle),
                       radial.y * std::sin(slot.spreadAngle),
                       std::cos(slot.spreadAngle)};

            // Gravity droop sag
            if (opts.droop > 0.0f) {
                dir = vnorm(dir * (1.0f - opts.droop) + Vec3{0.0f, -1.0f, 0.0f} * opts.droop);
            }
            // Phototropic lift toward sky
            if (opts.upBias > 0.0f) {
                dir = vnormOr(dir * (1.0f - opts.upBias * 0.35f) + worldUp * (opts.upBias * 0.35f), dir);
            }
        }

        // Petiole geometry
        Vec3 leafBase = attachPoint;
        float baseBend = t * 0.35f;
        float petioleTipBend = baseBend + 0.25f;

        if (opts.petioleLength > 1e-4f) {
            leafBase = attachPoint + dir * opts.petioleLength;
            float petioleRadius = opts.twigRadius * 0.65f;
            MeshData petiole = buildPetiole(attachPoint, leafBase, petioleRadius, baseBend, petioleTipBend);
            appendMesh(out, petiole);
        }

        // Leaf orientation coordinate frame:
        // F = forward along leaf length
        // N = normal / upward adaxial surface (biased toward +Y)
        // side = lateral width axis
        Vec3 F = dir;
        Vec3 side = vcross(F, worldUp);
        if (vdot(side, side) < 1e-6f) {
            side = vcross(F, Vec3{1.0f, 0.0f, 0.0f});
        }
        side = vnorm(side);
        Vec3 N = vnorm(vcross(side, F));

        // When upBias is not 100%, blend with radial normal for natural 3D fullness
        if (opts.upBias < 0.99f && slot.spreadAngle > 1e-4f) {
            Vec3 radialN{std::cos(slot.phi), std::sin(slot.phi), 0.0f};
            Vec3 radialNormal = vnormOr(vcross(vcross(F, radialN), F), N);
            N = vnormOr(N * opts.upBias + radialNormal * (1.0f - opts.upBias), N);
            side = vnorm(vcross(F, N));
        }

        // Transform base card into cluster space and append
        uint32_t vOffset = static_cast<uint32_t>(out.vertexCount());
        out.positions.reserve(out.positions.size() + cardVCount * 3);
        if (cardHasNormals) out.normals.reserve(out.normals.size() + cardVCount * 3);
        if (cardHasUVs) out.uvs.reserve(out.uvs.size() + cardVCount * 2);
        if (cardHasColors) out.colors.reserve(out.colors.size() + cardVCount * 4);

        for (size_t v = 0; v < cardVCount; ++v) {
            float lx = baseCard.positions[v * 3 + 0];
            float ly = baseCard.positions[v * 3 + 1];
            float lz = baseCard.positions[v * 3 + 2];

            Vec3 pos = leafBase + side * lx + N * ly + F * lz;
            out.positions.push_back(pos.x);
            out.positions.push_back(pos.y);
            out.positions.push_back(pos.z);

            if (cardHasNormals) {
                float nx = baseCard.normals[v * 3 + 0];
                float ny = baseCard.normals[v * 3 + 1];
                float nz = baseCard.normals[v * 3 + 2];
                Vec3 rawNorm = side * nx + N * ny + F * nz;
                Vec3 norm = vnormOr(rawNorm, N);
                out.normals.push_back(norm.x);
                out.normals.push_back(norm.y);
                out.normals.push_back(norm.z);
            }

            if (cardHasUVs) {
                out.uvs.push_back(baseCard.uvs[v * 2 + 0]);
                out.uvs.push_back(baseCard.uvs[v * 2 + 1]);
            }

            if (cardHasColors) {
                float localBend = baseCard.colors[v * 4 + 0];
                float finalBend = petioleTipBend + localBend * (1.0f - petioleTipBend);
                finalBend = std::clamp(finalBend, 0.0f, 1.0f);
                out.colors.push_back(finalBend);
                out.colors.push_back(0.0f);
                out.colors.push_back(0.0f);
                out.colors.push_back(1.0f);
            }
        }

        out.indices.reserve(out.indices.size() + cardICount);
        for (size_t i = 0; i < cardICount; ++i) {
            out.indices.push_back(baseCard.indices[i] + vOffset);
        }
    }

    return out;
}

} // namespace bromesh
