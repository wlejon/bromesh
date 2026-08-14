#include "test_framework.h"
#include <cmath>

#if BROMESH_HAS_MANIFOLD

static float rawVolume(const bromesh::MeshData& mesh) {
    double sum = 0.0;
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];
        double ax = mesh.positions[i0*3], ay = mesh.positions[i0*3+1], az = mesh.positions[i0*3+2];
        double bx = mesh.positions[i1*3], by = mesh.positions[i1*3+1], bz = mesh.positions[i1*3+2];
        double cx = mesh.positions[i2*3], cy = mesh.positions[i2*3+1], cz = mesh.positions[i2*3+2];
        sum += ax*(by*cz-bz*cy) + ay*(bz*cx-bx*cz) + az*(bx*cy-by*cx);
    }
    return static_cast<float>(std::fabs(sum) / 6.0);
}

TEST(boolean_union) {
    // Two overlapping spheres
    auto a = bromesh::sphere(1.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 1.0f;
    }

    auto result = bromesh::booleanUnion(a, b);
    ASSERT(!result.empty(), "bool_union: should produce non-empty result");
    ASSERT(result.triangleCount() > 0, "bool_union: should have triangles");
    float volA = rawVolume(a);
    float volResult = rawVolume(result);
    ASSERT(volResult > volA * 0.9f, "bool_union: result volume should exceed single sphere");
}

TEST(boolean_difference) {
    auto a = bromesh::sphere(2.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 1.5f;
    }

    auto result = bromesh::booleanDifference(a, b);
    ASSERT(!result.empty(), "bool_diff: should produce non-empty result");
    ASSERT(result.triangleCount() > 0, "bool_diff: should have triangles");
    // "Didn't just concatenate the inputs" is a volume property, not a vertex
    // count one: booleanOp runs crease-normal splitting on its output, which
    // duplicates vertices along the cut and can push the result past
    // a.vertexCount() + b.vertexCount(). Concatenating A and B would sum their
    // volumes; a real difference carves B's overlap out of A instead.
    float volA = rawVolume(a);
    float volResult = rawVolume(result);
    ASSERT(volResult > 0.0f, "bool_diff: result should enclose positive volume");
    ASSERT(volResult < volA - 0.01f,
           "bool_diff: should remove B's overlap from A, not concatenate");
    // BBox should fit within A's bbox
    auto bboxA = bromesh::computeBBox(a);
    auto bboxR = bromesh::computeBBox(result);
    ASSERT(bboxR.min.x >= bboxA.min.x - 0.01f && bboxR.max.x <= bboxA.max.x + 0.01f,
           "bool_diff: result should fit within A's bbox on X");
}

TEST(boolean_intersection) {
    auto a = bromesh::sphere(1.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 0.5f;
    }

    auto result = bromesh::booleanIntersection(a, b);
    ASSERT(!result.empty(), "bool_isect: should produce non-empty result");
    float volA = rawVolume(a);
    float volResult = rawVolume(result);
    ASSERT(volResult < volA * 1.01f, "bool_isect: intersection volume should be less than full sphere");
    ASSERT(volResult > 0.01f, "bool_isect: intersection should have meaningful volume");
}

TEST(boolean_no_overlap) {
    auto a = bromesh::sphere(1.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 5.0f;
    }

    auto result = bromesh::booleanUnion(a, b);
    ASSERT(!result.empty(), "bool_nooverlap: should produce result");
    float volA = rawVolume(a);
    float volB = rawVolume(b);
    float volResult = rawVolume(result);
    float expected = volA + volB;
    ASSERT(std::fabs(volResult - expected) < expected * 0.15f,
           "bool_nooverlap: union volume should be sum of parts");
}

TEST(split_by_plane) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);

    auto [top, bottom] = bromesh::splitByPlane(mesh, 0, 1, 0, 0);
    ASSERT(!top.empty(), "split_plane: top half should be non-empty");
    ASSERT(!bottom.empty(), "split_plane: bottom half should be non-empty");

    float topVol = rawVolume(top);
    float bottomVol = rawVolume(bottom);

    // Both halves should be roughly equal for a centered sphere split at Y=0
    ASSERT(std::fabs(topVol - bottomVol) < (topVol + bottomVol) * 0.2f,
           "split_plane: halves should be roughly equal");
}

TEST(boolean_box_minus_sphere) {
    auto cube = bromesh::box(1.5f, 1.5f, 1.5f);
    auto ball = bromesh::sphere(1.0f, 16, 12);

    auto result = bromesh::booleanDifference(cube, ball);
    ASSERT(!result.empty(), "bool_box_sphere: should produce result");
    ASSERT(result.triangleCount() > 0, "bool_box_sphere: should have triangles");
    // Result should be more complex than the original cube (added sphere boundary)
    ASSERT(result.triangleCount() > cube.triangleCount(),
           "bool_box_sphere: result should have more triangles than original cube");
    // BBox should fit within cube's bbox
    auto bboxC = bromesh::computeBBox(cube);
    auto bboxR = bromesh::computeBBox(result);
    ASSERT(bboxR.min.x >= bboxC.min.x - 0.01f && bboxR.max.x <= bboxC.max.x + 0.01f,
           "bool_box_sphere: result should fit within cube bbox");
}

TEST(boolean_disjoint_cubes_exact_volume) {
    // Cube A: size 1x1x1 at [-2.5, -1.5]^3 (center -2.0, volume 1.0)
    auto a = bromesh::box(0.5f, 0.5f, 0.5f);
    bromesh::translateMesh(a, -2.0f, -2.0f, -2.0f);

    // Cube B: size 1x1x1 at [1.5, 2.5]^3 (center +2.0, volume 1.0)
    auto b = bromesh::box(0.5f, 0.5f, 0.5f);
    bromesh::translateMesh(b, 2.0f, 2.0f, 2.0f);

    auto un = bromesh::booleanUnion(a, b);
    ASSERT(std::fabs(rawVolume(un) - 2.0f) < 1e-3f, "disjoint union volume == 2.0");

    auto isect = bromesh::booleanIntersection(a, b);
    ASSERT(isect.empty() || isect.triangleCount() == 0 || rawVolume(isect) < 1e-3f,
           "disjoint intersection is empty or volume == 0.0");

    auto diff = bromesh::booleanDifference(a, b);
    ASSERT(std::fabs(rawVolume(diff) - 1.0f) < 1e-3f, "disjoint difference volume == 1.0");
}

TEST(boolean_concentric_cubes_exact_volume) {
    // Cube A: size 2x2x2 centered at origin (volume 8.0)
    auto a = bromesh::box(1.0f, 1.0f, 1.0f);
    // Cube B: size 1x1x1 centered at origin (volume 1.0)
    auto b = bromesh::box(0.5f, 0.5f, 0.5f);

    auto un = bromesh::booleanUnion(a, b);
    ASSERT(std::fabs(rawVolume(un) - 8.0f) < 1e-3f, "concentric union volume == 8.0");

    auto isect = bromesh::booleanIntersection(a, b);
    ASSERT(std::fabs(rawVolume(isect) - 1.0f) < 1e-3f, "concentric intersection volume == 1.0");

    auto diff = bromesh::booleanDifference(a, b);
    ASSERT(std::fabs(rawVolume(diff) - 7.0f) < 1e-3f, "concentric difference volume == 7.0");
}

TEST(boolean_overlapping_cubes_exact_volume) {
    // Box A at [-1, 1]^3 (size 2x2x2, volume 8.0)
    auto a = bromesh::box(1.0f, 1.0f, 1.0f);
    // Box B at [0, 2] x [-1, 1] x [-1, 1] (size 2x2x2, shifted +1 in X, volume 8.0)
    auto b = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::translateMesh(b, 1.0f, 0.0f, 0.0f);

    // Overlap is [0, 1] x [-1, 1] x [-1, 1] (volume 1 x 2 x 2 = 4.0)
    auto un = bromesh::booleanUnion(a, b);
    ASSERT(std::fabs(rawVolume(un) - 12.0f) < 1e-3f, "overlapping union volume == 12.0");

    auto isect = bromesh::booleanIntersection(a, b);
    ASSERT(std::fabs(rawVolume(isect) - 4.0f) < 1e-3f, "overlapping intersection volume == 4.0");

    auto diff = bromesh::booleanDifference(a, b);
    ASSERT(std::fabs(rawVolume(diff) - 4.0f) < 1e-3f, "overlapping difference volume == 4.0");
}

TEST(boolean_sphere_sphere_analytic_lens) {
    // Two spheres of radius R=1.0 with center distance d=1.0
    // Analytic lens volume V_lens = pi * (4R + d) * (2R - d)^2 / 12 = 5*pi/12 ~= 1.308997
    // Analytic sphere volume V_s = 4/3 * pi * R^3 ~= 4.188790
    auto a = bromesh::sphere(1.0f, 48, 32);
    auto b = bromesh::sphere(1.0f, 48, 32);
    bromesh::translateMesh(b, 1.0f, 0.0f, 0.0f);

    const float expectedLens = 1.3089969f;
    const float expectedSphere = 4.1887902f;
    const float expectedUnion = 2.0f * expectedSphere - expectedLens; // ~= 7.0685835
    const float expectedDiff = expectedSphere - expectedLens;        // ~= 2.8797933

    auto isect = bromesh::booleanIntersection(a, b);
    float volIsect = rawVolume(isect);
    ASSERT(std::fabs(volIsect - expectedLens) < 0.05f * expectedLens,
           "sphere lens intersection volume matches 5pi/12 within 5%");

    auto un = bromesh::booleanUnion(a, b);
    float volUnion = rawVolume(un);
    ASSERT(std::fabs(volUnion - expectedUnion) < 0.05f * expectedUnion,
           "sphere lens union volume matches 2*Vs - Vlens within 5%");

    auto diff = bromesh::booleanDifference(a, b);
    float volDiff = rawVolume(diff);
    ASSERT(std::fabs(volDiff - expectedDiff) < 0.05f * expectedDiff,
           "sphere lens difference volume matches Vs - Vlens within 5%");
}

TEST(polygon2d_square) {
    // A unit CCW square must triangulate into 2 triangles, 4 verts, area 1.
    std::vector<float> outer = {
        0.0f, 0.0f,
        1.0f, 0.0f,
        1.0f, 1.0f,
        0.0f, 1.0f,
    };
    auto m = bromesh::triangulatePolygon2D(outer);
    ASSERT(m.vertexCount() == 4,     "polygon2d_square: 4 verts");
    ASSERT(m.triangleCount() == 2,   "polygon2d_square: 2 tris");
    ASSERT(m.indices.size() == 6,    "polygon2d_square: 6 index entries");
    ASSERT(m.normals.size() == 12,   "polygon2d_square: normals populated");
    // All vertices in the z=0 plane.
    for (size_t i = 0; i < 4; ++i) {
        ASSERT(std::fabs(m.positions[i*3 + 2]) < 1e-6f, "polygon2d_square: z==0");
    }
    // All normals should be +Z (CCW outer → front-face toward +Z).
    for (size_t i = 0; i < 4; ++i) {
        ASSERT(std::fabs(m.normals[i*3 + 2] - 1.0f) < 1e-6f,
               "polygon2d_square: +Z normals");
    }
}

TEST(polygon2d_cw_flips_normal) {
    // Same square, wound CW — normal should flip to -Z.
    std::vector<float> outer = {
        0.0f, 0.0f,
        0.0f, 1.0f,
        1.0f, 1.0f,
        1.0f, 0.0f,
    };
    auto m = bromesh::triangulatePolygon2D(outer);
    ASSERT(m.vertexCount() == 4, "polygon2d_cw: 4 verts");
    for (size_t i = 0; i < 4; ++i) {
        ASSERT(std::fabs(m.normals[i*3 + 2] + 1.0f) < 1e-6f,
               "polygon2d_cw: -Z normals");
    }
}

TEST(polygon2d_with_hole) {
    // 2x2 outer square (CCW) with a 1x1 hole in the middle (CW).
    std::vector<float> outer = {
        0.0f, 0.0f,
        2.0f, 0.0f,
        2.0f, 2.0f,
        0.0f, 2.0f,
    };
    std::vector<std::vector<float>> holes = {{
        0.5f, 0.5f,
        0.5f, 1.5f,
        1.5f, 1.5f,
        1.5f, 0.5f,
    }};
    auto m = bromesh::triangulatePolygon2D(outer, holes);
    ASSERT(m.vertexCount() == 8, "polygon2d_hole: 4+4 verts");
    // Any valid triangulation of a square with a square hole has 8 tris
    // (outer perimeter forms a ring around the hole).
    ASSERT(m.triangleCount() == 8, "polygon2d_hole: 8 tris");
    // Spot-check area by summing triangle areas — should equal 3 (4 - 1).
    double area = 0.0;
    for (size_t t = 0; t < m.indices.size(); t += 3) {
        const uint32_t a = m.indices[t+0];
        const uint32_t b = m.indices[t+1];
        const uint32_t c = m.indices[t+2];
        const double ax = m.positions[a*3+0], ay = m.positions[a*3+1];
        const double bx = m.positions[b*3+0], by = m.positions[b*3+1];
        const double cx = m.positions[c*3+0], cy = m.positions[c*3+1];
        area += 0.5 * std::fabs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay));
    }
    ASSERT(std::fabs(area - 3.0) < 1e-4, "polygon2d_hole: area == 3");
}

TEST(polygon2d_degenerate_rejected) {
    // Fewer than 3 vertices → empty mesh.
    auto m1 = bromesh::triangulatePolygon2D({0.0f, 0.0f, 1.0f, 1.0f});
    ASSERT(m1.empty(), "polygon2d_degenerate: 2-pt outer rejected");
    // Empty outer.
    auto m2 = bromesh::triangulatePolygon2D({});
    ASSERT(m2.empty(), "polygon2d_degenerate: empty outer rejected");
}

TEST(polygon2d_concave_L) {
    // L-shaped polygon — concave, 6 verts, must triangulate without leaving
    // the L's interior.
    std::vector<float> outer = {
        0.0f, 0.0f,
        2.0f, 0.0f,
        2.0f, 1.0f,
        1.0f, 1.0f,
        1.0f, 2.0f,
        0.0f, 2.0f,
    };
    auto m = bromesh::triangulatePolygon2D(outer);
    ASSERT(m.vertexCount() == 6, "polygon2d_L: 6 verts");
    ASSERT(m.triangleCount() == 4, "polygon2d_L: 4 tris (n-2 for simple polygon)");
    double area = 0.0;
    for (size_t t = 0; t < m.indices.size(); t += 3) {
        const uint32_t a = m.indices[t+0];
        const uint32_t b = m.indices[t+1];
        const uint32_t c = m.indices[t+2];
        const double ax = m.positions[a*3+0], ay = m.positions[a*3+1];
        const double bx = m.positions[b*3+0], by = m.positions[b*3+1];
        const double cx = m.positions[c*3+0], cy = m.positions[c*3+1];
        area += 0.5 * std::fabs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay));
    }
    ASSERT(std::fabs(area - 3.0) < 1e-4, "polygon2d_L: area == 3");
}

TEST(polygon3d_xz_plane) {
    // Unit square on the XZ plane (normal = +Y). Output positions must match
    // the input verbatim.
    std::vector<float> outer = {
        0.0f, 5.0f, 0.0f,
        1.0f, 5.0f, 0.0f,
        1.0f, 5.0f, 1.0f,
        0.0f, 5.0f, 1.0f,
    };
    const float n[3] = {0.0f, 1.0f, 0.0f};
    auto m = bromesh::triangulatePolygon3D(outer, {}, n);
    ASSERT(m.vertexCount() == 4, "polygon3d_xz: 4 verts");
    ASSERT(m.triangleCount() == 2, "polygon3d_xz: 2 tris");
    for (size_t i = 0; i < 4; ++i) {
        ASSERT(std::fabs(m.positions[i*3 + 0] - outer[i*3 + 0]) < 1e-6f,
               "polygon3d_xz: x preserved");
        ASSERT(std::fabs(m.positions[i*3 + 1] - outer[i*3 + 1]) < 1e-6f,
               "polygon3d_xz: y preserved");
        ASSERT(std::fabs(m.positions[i*3 + 2] - outer[i*3 + 2]) < 1e-6f,
               "polygon3d_xz: z preserved");
        ASSERT(std::fabs(m.normals[i*3 + 1] - 1.0f) < 1e-6f,
               "polygon3d_xz: +Y normals");
    }
}

TEST(polygon3d_tilted_plane) {
    // Square on a plane normal to (1,1,1)/sqrt(3). Verify area matches
    // (hand-computed: 4 verts → 2 right triangles → area 1).
    const float s = 0.7071067811865476f;  // 1/sqrt(2)
    // Points on a plane through the origin with normal (1,1,1)/sqrt(3):
    // easiest: pick two orthonormal tangent vectors and emit a unit square.
    const float u[3] = { s, -s, 0 };
    const float v[3] = { s,  s, -2.0f * s };
    // Normalize v. (s,s,-2s) has length sqrt(s^2+s^2+4s^2) = sqrt(6)s. Divide.
    const float vl = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    const float vv[3] = { v[0]/vl, v[1]/vl, v[2]/vl };
    auto pt = [&](float a, float b) {
        return std::vector<float>{
            a*u[0] + b*vv[0],
            a*u[1] + b*vv[1],
            a*u[2] + b*vv[2],
        };
    };
    const auto p0 = pt(0, 0), p1 = pt(1, 0), p2 = pt(1, 1), p3 = pt(0, 1);
    std::vector<float> outer = {
        p0[0], p0[1], p0[2],
        p1[0], p1[1], p1[2],
        p2[0], p2[1], p2[2],
        p3[0], p3[1], p3[2],
    };
    const float inv3 = 0.5773502691896258f;
    const float n[3] = {inv3, inv3, inv3};
    auto m = bromesh::triangulatePolygon3D(outer, {}, n);
    ASSERT(m.triangleCount() == 2, "polygon3d_tilted: 2 tris");
    double area = 0.0;
    for (size_t t = 0; t < m.indices.size(); t += 3) {
        const uint32_t ai = m.indices[t+0];
        const uint32_t bi = m.indices[t+1];
        const uint32_t ci = m.indices[t+2];
        const double ax = m.positions[ai*3+0], ay = m.positions[ai*3+1], az = m.positions[ai*3+2];
        const double bx = m.positions[bi*3+0], by = m.positions[bi*3+1], bz = m.positions[bi*3+2];
        const double cx = m.positions[ci*3+0], cy = m.positions[ci*3+1], cz = m.positions[ci*3+2];
        const double ex = bx-ax, ey = by-ay, ez = bz-az;
        const double fx = cx-ax, fy = cy-ay, fz = cz-az;
        const double crx = ey*fz - ez*fy;
        const double cry = ez*fx - ex*fz;
        const double crz = ex*fy - ey*fx;
        area += 0.5 * std::sqrt(crx*crx + cry*cry + crz*crz);
    }
    ASSERT(std::fabs(area - 1.0) < 1e-4, "polygon3d_tilted: area == 1");
}

#endif // BROMESH_HAS_MANIFOLD
