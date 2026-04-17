#include "bromesh/manipulation/polygon.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#if BROMESH_HAS_MANIFOLD
#include <manifold/polygon.h>
#endif

namespace bromesh {

namespace {

// Pick a 2D basis (u, v) on the plane defined by unit-length `n`. `u` is
// chosen as `n` × ref where ref is the world axis least parallel to `n`;
// `v` = n × u. The result: any 3D point p on the plane maps to
// (dot(p-origin, u), dot(p-origin, v)).
void planeBasis(const float n[3], float u[3], float v[3]) {
    float ref[3];
    const float ax = std::fabs(n[0]);
    const float ay = std::fabs(n[1]);
    const float az = std::fabs(n[2]);
    if (ax <= ay && ax <= az)      { ref[0] = 1; ref[1] = 0; ref[2] = 0; }
    else if (ay <= az)             { ref[0] = 0; ref[1] = 1; ref[2] = 0; }
    else                           { ref[0] = 0; ref[1] = 0; ref[2] = 1; }
    u[0] = n[1] * ref[2] - n[2] * ref[1];
    u[1] = n[2] * ref[0] - n[0] * ref[2];
    u[2] = n[0] * ref[1] - n[1] * ref[0];
    const float lu = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    if (lu > 1e-20f) {
        u[0] /= lu; u[1] /= lu; u[2] /= lu;
    } else {
        u[0] = 1; u[1] = 0; u[2] = 0;
    }
    v[0] = n[1]*u[2] - n[2]*u[1];
    v[1] = n[2]*u[0] - n[0]*u[2];
    v[2] = n[0]*u[1] - n[1]*u[0];
}

// Signed area of a 2D polygon (shoelace). Positive = CCW.
double signedArea(const std::vector<float>& xy) {
    const size_t n = xy.size() / 2;
    if (n < 3) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const size_t j = (i + 1) % n;
        sum += double(xy[i*2 + 0]) * double(xy[j*2 + 1]);
        sum -= double(xy[j*2 + 0]) * double(xy[i*2 + 1]);
    }
    return 0.5 * sum;
}

} // anonymous namespace


MeshData triangulatePolygon2D(
    const std::vector<float>& outer,
    const std::vector<std::vector<float>>& holes,
    float z)
{
#if BROMESH_HAS_MANIFOLD
    const size_t outerN = outer.size() / 2;
    if (outer.size() < 6 || (outer.size() % 2) != 0) return {};

    // Build manifold's Polygons (first contour is outer, rest are holes).
    manifold::Polygons polys;
    polys.reserve(1 + holes.size());
    manifold::SimplePolygon simpleOuter;
    simpleOuter.reserve(outerN);
    for (size_t i = 0; i < outerN; ++i) {
        simpleOuter.push_back({double(outer[i*2 + 0]), double(outer[i*2 + 1])});
    }
    polys.push_back(std::move(simpleOuter));
    for (const auto& h : holes) {
        if (h.size() < 6 || (h.size() % 2) != 0) continue;
        const size_t hn = h.size() / 2;
        manifold::SimplePolygon sh;
        sh.reserve(hn);
        for (size_t i = 0; i < hn; ++i) {
            sh.push_back({double(h[i*2 + 0]), double(h[i*2 + 1])});
        }
        polys.push_back(std::move(sh));
    }

    std::vector<manifold::ivec3> tris;
    try {
        tris = manifold::Triangulate(polys);
    } catch (...) {
        return {};
    }
    if (tris.empty()) return {};

    // Flatten the polygons (outer + holes, in manifold's index order) back
    // into our vertex buffer. manifold::Triangulate numbers vertices
    // sequentially in the input's traversal order, so we do the same here.
    MeshData out;
    size_t totalVerts = outerN;
    for (const auto& h : holes) {
        if (h.size() >= 6 && (h.size() % 2) == 0) totalVerts += h.size() / 2;
    }
    out.positions.reserve(totalVerts * 3);
    out.normals.reserve(totalVerts * 3);

    // Front-face direction is +z if the outer contour is CCW, else -z.
    const double area = signedArea(outer);
    const float nz = (area >= 0.0) ? 1.0f : -1.0f;

    for (size_t i = 0; i < outerN; ++i) {
        out.positions.push_back(outer[i*2 + 0]);
        out.positions.push_back(outer[i*2 + 1]);
        out.positions.push_back(z);
        out.normals.push_back(0.0f);
        out.normals.push_back(0.0f);
        out.normals.push_back(nz);
    }
    for (const auto& h : holes) {
        if (h.size() < 6 || (h.size() % 2) != 0) continue;
        const size_t hn = h.size() / 2;
        for (size_t i = 0; i < hn; ++i) {
            out.positions.push_back(h[i*2 + 0]);
            out.positions.push_back(h[i*2 + 1]);
            out.positions.push_back(z);
            out.normals.push_back(0.0f);
            out.normals.push_back(0.0f);
            out.normals.push_back(nz);
        }
    }

    out.indices.reserve(tris.size() * 3);
    for (const auto& t : tris) {
        out.indices.push_back(static_cast<uint32_t>(t.x));
        out.indices.push_back(static_cast<uint32_t>(t.y));
        out.indices.push_back(static_cast<uint32_t>(t.z));
    }
    return out;
#else
    (void)outer; (void)holes; (void)z;
    return {};
#endif
}


MeshData triangulatePolygon3D(
    const std::vector<float>& outer,
    const std::vector<std::vector<float>>& holes,
    const float normal[3])
{
#if BROMESH_HAS_MANIFOLD
    if (outer.size() < 9 || (outer.size() % 3) != 0) return {};
    const size_t outerN = outer.size() / 3;

    // Project 3D input onto a 2D basis on the plane. The origin is the
    // first outer vertex — everything is relative to it to keep projected
    // coordinates near the origin regardless of how far the polygon is from
    // world zero.
    float u[3], v[3];
    planeBasis(normal, u, v);
    const float ox = outer[0], oy = outer[1], oz = outer[2];

    auto projectContour = [&](const std::vector<float>& src,
                              std::vector<float>& out2d) {
        const size_t n = src.size() / 3;
        out2d.resize(n * 2);
        for (size_t i = 0; i < n; ++i) {
            const float dx = src[i*3 + 0] - ox;
            const float dy = src[i*3 + 1] - oy;
            const float dz = src[i*3 + 2] - oz;
            out2d[i*2 + 0] = dx * u[0] + dy * u[1] + dz * u[2];
            out2d[i*2 + 1] = dx * v[0] + dy * v[1] + dz * v[2];
        }
    };

    std::vector<float> outer2d;
    projectContour(outer, outer2d);
    std::vector<std::vector<float>> holes2d;
    holes2d.reserve(holes.size());
    for (const auto& h : holes) {
        if (h.size() < 9 || (h.size() % 3) != 0) continue;
        std::vector<float> h2d;
        projectContour(h, h2d);
        holes2d.push_back(std::move(h2d));
    }

    // Manifold's Triangulate expects the outer contour CCW and holes CW.
    // Our projection basis may orient the caller's CCW-in-world input as CW
    // in (u, v) — in which case we reverse both contours (swapping tri
    // indices post-triangulation to restore world winding). Without this,
    // non-convex inputs come back with mixed tri winding.
    const bool reversed = signedArea(outer2d) < 0.0;
    if (reversed) {
        std::reverse(outer2d.begin(), outer2d.end());
        // `outer2d` is a flat [x0,y0, x1,y1, ...]; the reverse above flipped
        // into [..., y1,x1, y0,x0]. Re-pair the coords.
        for (size_t i = 0; i + 1 < outer2d.size(); i += 2) {
            std::swap(outer2d[i], outer2d[i+1]);
        }
        for (auto& h2d : holes2d) {
            std::reverse(h2d.begin(), h2d.end());
            for (size_t i = 0; i + 1 < h2d.size(); i += 2) {
                std::swap(h2d[i], h2d[i+1]);
            }
        }
    }

    // Run through the 2D core, then rewrite positions with the original
    // 3D coordinates and swap normals to the caller-supplied vector.
    MeshData mesh2d = triangulatePolygon2D(outer2d, holes2d, 0.0f);
    if (mesh2d.indices.empty()) return {};

    // If we reversed the 2D input, manifold's tri indices point into the
    // reversed vertex order — remap back to original 3D index order. For
    // the outer contour, index i in the reversed layout corresponds to
    // (outerN - 1 - i) in the original. Holes follow the outer block.
    if (reversed) {
        std::vector<uint32_t> segStart = { 0, static_cast<uint32_t>(outerN) };
        std::vector<uint32_t> segLen   = { static_cast<uint32_t>(outerN) };
        for (const auto& h : holes) {
            if (h.size() < 9 || (h.size() % 3) != 0) continue;
            const uint32_t hn = static_cast<uint32_t>(h.size() / 3);
            segLen.push_back(hn);
            segStart.push_back(segStart.back() + hn);
        }
        auto remap = [&](uint32_t idx) -> uint32_t {
            for (size_t s = 0; s < segLen.size(); ++s) {
                const uint32_t base = segStart[s];
                const uint32_t len  = segLen[s];
                if (idx >= base && idx < base + len) {
                    return base + (len - 1 - (idx - base));
                }
            }
            return idx;
        };
        for (uint32_t& i : mesh2d.indices) i = remap(i);
        // Reversing the input also flipped tri winding; swap to restore it.
        for (size_t t = 0; t + 2 < mesh2d.indices.size(); t += 3) {
            std::swap(mesh2d.indices[t + 1], mesh2d.indices[t + 2]);
        }
    }

    const size_t vCount = mesh2d.vertexCount();

    // Assemble a flat 3D input that matches manifold's outer+holes ordering.
    std::vector<float> pos3d;
    pos3d.reserve(vCount * 3);
    for (size_t i = 0; i < outerN; ++i) {
        pos3d.push_back(outer[i*3 + 0]);
        pos3d.push_back(outer[i*3 + 1]);
        pos3d.push_back(outer[i*3 + 2]);
    }
    for (const auto& h : holes) {
        if (h.size() < 9 || (h.size() % 3) != 0) continue;
        const size_t hn = h.size() / 3;
        for (size_t i = 0; i < hn; ++i) {
            pos3d.push_back(h[i*3 + 0]);
            pos3d.push_back(h[i*3 + 1]);
            pos3d.push_back(h[i*3 + 2]);
        }
    }
    if (pos3d.size() != vCount * 3) return {};   // holes were filtered mid-build

    MeshData out;
    out.positions = std::move(pos3d);
    out.normals.resize(vCount * 3);
    for (size_t i = 0; i < vCount; ++i) {
        out.normals[i*3 + 0] = normal[0];
        out.normals[i*3 + 1] = normal[1];
        out.normals[i*3 + 2] = normal[2];
    }
    out.indices = std::move(mesh2d.indices);
    return out;
#else
    (void)outer; (void)holes; (void)normal;
    return {};
#endif
}

} // namespace bromesh
