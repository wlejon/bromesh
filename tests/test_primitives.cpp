#include "test_framework.h"
#include <cmath>

static bool checkOutwardWinding(const bromesh::MeshData& m, float cx, float cy, float cz) {
    size_t triCount = m.triangleCount();
    for (size_t t = 0; t < triCount; t++) {
        uint32_t i0 = m.indices[t * 3 + 0];
        uint32_t i1 = m.indices[t * 3 + 1];
        uint32_t i2 = m.indices[t * 3 + 2];

        float ax = m.positions[i0 * 3 + 0], ay = m.positions[i0 * 3 + 1], az = m.positions[i0 * 3 + 2];
        float bx = m.positions[i1 * 3 + 0], by = m.positions[i1 * 3 + 1], bz = m.positions[i1 * 3 + 2];
        float ex = m.positions[i2 * 3 + 0], ey = m.positions[i2 * 3 + 1], ez = m.positions[i2 * 3 + 2];

        float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
        float e2x = ex - ax, e2y = ey - ay, e2z = ez - az;

        float nx = e1y * e2z - e1z * e2y;
        float ny = e1z * e2x - e1x * e2z;
        float nz = e1x * e2y - e1y * e2x;

        float dx = (ax + bx + ex) / 3.0f - cx;
        float dy = (ay + by + ey) / 3.0f - cy;
        float dz = (az + bz + ez) / 3.0f - cz;

        float dot = nx * dx + ny * dy + nz * dz;
        if (dot < -1e-5f) return false;
    }
    return true;
}

static bool checkTorusWinding(const bromesh::MeshData& m, float majorRadius) {
    size_t triCount = m.triangleCount();
    for (size_t t = 0; t < triCount; t++) {
        uint32_t i0 = m.indices[t * 3 + 0];
        uint32_t i1 = m.indices[t * 3 + 1];
        uint32_t i2 = m.indices[t * 3 + 2];

        float ax = m.positions[i0 * 3 + 0], ay = m.positions[i0 * 3 + 1], az = m.positions[i0 * 3 + 2];
        float bx = m.positions[i1 * 3 + 0], by = m.positions[i1 * 3 + 1], bz = m.positions[i1 * 3 + 2];
        float ex = m.positions[i2 * 3 + 0], ey = m.positions[i2 * 3 + 1], ez = m.positions[i2 * 3 + 2];

        float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
        float e2x = ex - ax, e2y = ey - ay, e2z = ez - az;

        float fnx = e1y * e2z - e1z * e2y;
        float fny = e1z * e2x - e1x * e2z;
        float fnz = e1x * e2y - e1y * e2x;

        float cx = (ax + bx + ex) / 3.0f;
        float cy = (ay + by + ey) / 3.0f;
        float cz = (az + bz + ez) / 3.0f;

        float len = std::sqrt(cx * cx + cz * cz);
        float tcx = (cx / len) * majorRadius;
        float tcz = (cz / len) * majorRadius;

        float dx = cx - tcx, dy = cy, dz = cz - tcz;
        float dot = fnx * dx + fny * dy + fnz * dz;
        if (dot < -1e-5f) return false;
    }
    return true;
}

TEST(stubs_link) {
    // Just verify all stub functions link without crashing
    float field[8] = {-1,-1,-1,-1, 1,1,1,1};
    auto mc = bromesh::marchingCubes(field, 2, 2, 2);
    ASSERT(!mc.empty(), "marching cubes with crossing returns non-empty");

    auto sn = bromesh::surfaceNets(field, 2, 2, 2);
    ASSERT(!sn.empty(), "surface nets with crossing returns non-empty");

    auto dc = bromesh::dualContour(field, 2, 2, 2);
    ASSERT(!dc.empty(), "dual contour with crossing returns non-empty");

    uint8_t voxels[8] = {1,1,1,1, 0,0,0,0};
    auto gm = bromesh::greedyMesh(voxels, 2, 2, 2);
    ASSERT(!gm.empty(), "greedy mesh returns non-empty for partial volume");

    auto b = bromesh::box(1, 1, 1);
    ASSERT(!b.empty(), "box returns non-empty mesh");
    ASSERT(b.vertexCount() == 24, "box has 24 vertices");
    ASSERT(b.triangleCount() == 12, "box has 12 triangles");
    ASSERT(b.hasNormals(), "box has normals");
    ASSERT(b.hasUVs(), "box has UVs");

    auto s = bromesh::sphere(1);
    ASSERT(!s.empty(), "sphere returns non-empty mesh");
    ASSERT(s.hasNormals(), "sphere has normals");
    ASSERT(s.hasUVs(), "sphere has UVs");

    auto cyl = bromesh::cylinder(1, 2, 16);
    ASSERT(!cyl.empty(), "cylinder returns non-empty mesh");
    ASSERT(cyl.hasNormals(), "cylinder has normals");

    auto cap = bromesh::capsule(1, 1, 16, 8);
    ASSERT(!cap.empty(), "capsule returns non-empty mesh");
    ASSERT(cap.hasNormals(), "capsule has normals");

    auto pl = bromesh::plane(1, 1, 2, 2);
    ASSERT(!pl.empty(), "plane returns non-empty mesh");
    ASSERT(pl.vertexCount() == 9, "plane 2x2 has 9 vertices");
    ASSERT(pl.triangleCount() == 8, "plane 2x2 has 8 triangles");

    auto tor = bromesh::torus(2, 0.5f, 24, 12);
    ASSERT(!tor.empty(), "torus returns non-empty mesh");
    ASSERT(tor.hasNormals(), "torus has normals");

    float heights[9] = {0,0,0, 0,1,0, 0,0,0};
    auto hm = bromesh::heightmapGrid(heights, 3, 3, 1.0f);
    ASSERT(!hm.empty(), "heightmap returns non-empty mesh");
    ASSERT(hm.vertexCount() == 9, "heightmap 3x3 has 9 vertices");
    ASSERT(hm.hasNormals(), "heightmap has normals");
}


#if BROMESH_HAS_PAR_SHAPES
TEST(par_icosahedron) {
    auto mesh = bromesh::icosahedron();
    ASSERT(!mesh.empty(), "par_ico: should be non-empty");
    ASSERT(mesh.triangleCount() == 20, "par_ico: icosahedron has 20 faces");
    ASSERT(mesh.vertexCount() == 12, "par_ico: icosahedron has 12 vertices");
}

TEST(par_dodecahedron) {
    auto mesh = bromesh::dodecahedron();
    ASSERT(!mesh.empty(), "par_dodec: should be non-empty");
    ASSERT(mesh.vertexCount() > 0, "par_dodec: should have vertices");
    ASSERT(mesh.triangleCount() > 0, "par_dodec: should have triangles");
}

TEST(par_octahedron) {
    auto mesh = bromesh::octahedron();
    ASSERT(!mesh.empty(), "par_oct: should be non-empty");
    ASSERT(mesh.triangleCount() == 8, "par_oct: octahedron has 8 faces");
    ASSERT(mesh.vertexCount() == 6, "par_oct: octahedron has 6 vertices");
}

TEST(par_tetrahedron) {
    auto mesh = bromesh::tetrahedron();
    ASSERT(!mesh.empty(), "par_tet: should be non-empty");
    ASSERT(mesh.triangleCount() == 4, "par_tet: tetrahedron has 4 faces");
    ASSERT(mesh.vertexCount() == 4, "par_tet: tetrahedron has 4 vertices");
}

TEST(par_geodesic_sphere) {
    auto mesh = bromesh::geodesicSphere(2.0f, 2);
    ASSERT(!mesh.empty(), "par_geod: should be non-empty");
    ASSERT(mesh.vertexCount() > 12, "par_geod: subdivided should have more than icosahedron");
    ASSERT(mesh.hasNormals(), "par_geod: should have normals");

    // Verify all vertices are approximately on the sphere surface
    bool allOnSurface = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float x = mesh.positions[v * 3 + 0];
        float y = mesh.positions[v * 3 + 1];
        float z = mesh.positions[v * 3 + 2];
        float dist = std::sqrt(x * x + y * y + z * z);
        if (std::fabs(dist - 2.0f) > 0.1f) { allOnSurface = false; break; }
    }
    ASSERT(allOnSurface, "par_geod: vertices should be on sphere surface");
}

TEST(par_geodesic_sphere_topology) {
    // A geodesic sphere at subdivision n has EXACTLY 10*4^n + 2 vertices and
    // 20*4^n triangles. par_shapes builds it as an unwelded soup and welds the
    // bit-identical duplicates afterwards; its stock weld epsilon exceeded the
    // edge length past subdivision 6, merging DISTINCT vertices and dropping
    // the collapsed triangles — a planet-scale sphere came back with 11% of
    // its triangles, sliver facets and holes. The epsilon now scales with the
    // subdivision order; this pins the exact counts at the resolutions that
    // used to break (7, 8) and one that never did (5).
    for (int n : {5, 7, 8}) {
        auto mesh = bromesh::geodesicSphere(1.0f, n);
        size_t expV = 10u * ((size_t)1 << (2 * n)) + 2u;
        size_t expT = 20u * ((size_t)1 << (2 * n));
        ASSERT(mesh.vertexCount() == expV,
               "par_geod_topo: vertex count must be exactly 10*4^n + 2");
        ASSERT(mesh.triangleCount() == expT,
               "par_geod_topo: triangle count must be exactly 20*4^n");
    }
}

TEST(par_cone) {
    const float r = 1.0f, h = 2.0f;
    auto mesh = bromesh::cone(r, h, 16, 4);
    ASSERT(!mesh.empty(), "par_cone: should be non-empty");
    ASSERT(mesh.vertexCount() > 10, "par_cone: should have vertices");
    ASSERT(mesh.triangleCount() > 10, "par_cone: should have triangles");

    // Canonical layout: base disc at Y=0 with radius `r`, apex at Y=`h`.
    float minX = 1e9f, maxX = -1e9f, minY = 1e9f, maxY = -1e9f, minZ = 1e9f, maxZ = -1e9f;
    for (size_t i = 0; i < mesh.vertexCount(); ++i) {
        float x = mesh.positions[i * 3 + 0];
        float y = mesh.positions[i * 3 + 1];
        float z = mesh.positions[i * 3 + 2];
        if (x < minX) minX = x; if (x > maxX) maxX = x;
        if (y < minY) minY = y; if (y > maxY) maxY = y;
        if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
    }
    ASSERT(std::fabs(minY) < 1e-4f, "par_cone: base sits at Y=0");
    ASSERT(std::fabs(maxY - h) < 1e-4f, "par_cone: apex sits at Y=height");
    ASSERT(std::fabs(maxX - r) < 1e-4f && std::fabs(minX + r) < 1e-4f,
           "par_cone: X spans [-radius, radius]");
    ASSERT(std::fabs(maxZ - r) < 1e-4f && std::fabs(minZ + r) < 1e-4f,
           "par_cone: Z spans [-radius, radius]");

    // Capped variant: triangle count grows by `slices` (one per fan slice).
    auto open   = bromesh::cone(r, h, 16, 4, false);
    auto capped = bromesh::cone(r, h, 16, 4, true);
    ASSERT(capped.triangleCount() == open.triangleCount() + 16,
           "par_cone: capBase adds `slices` fan triangles");
    ASSERT(capped.vertexCount() == open.vertexCount() + 17,
           "par_cone: capBase adds 1 center + `slices` ring vertices");

    // Lateral surface must face outward. Sample any side vertex (not the
    // apex, which has |x|+|z|≈0): its normal should point away from the
    // central Y axis. Triangle winding should also be outward — cross of
    // edges should point away from the centroid's lateral component.
    bool anyOutwardN = false;
    for (size_t i = 0; i < open.vertexCount(); ++i) {
        float x = open.positions[i * 3 + 0];
        float z = open.positions[i * 3 + 2];
        float r2 = x*x + z*z;
        if (r2 < 1e-3f) continue;
        float nx = open.normals[i * 3 + 0];
        float nz = open.normals[i * 3 + 2];
        float dot = nx * x + nz * z;
        ASSERT(dot >= -1e-4f, "par_cone: side normal not inward");
        if (dot > 1e-3f) anyOutwardN = true;
    }
    ASSERT(anyOutwardN, "par_cone: at least one side normal points outward");

    bool anyOutwardW = false;
    for (size_t t = 0; t < open.triangleCount(); ++t) {
        uint32_t ai = open.indices[t * 3 + 0];
        uint32_t bi = open.indices[t * 3 + 1];
        uint32_t ci = open.indices[t * 3 + 2];
        float ax = open.positions[ai * 3 + 0], ay = open.positions[ai * 3 + 1], az = open.positions[ai * 3 + 2];
        float bx = open.positions[bi * 3 + 0], by = open.positions[bi * 3 + 1], bz = open.positions[bi * 3 + 2];
        float cx = open.positions[ci * 3 + 0], cy = open.positions[ci * 3 + 1], cz = open.positions[ci * 3 + 2];
        float ux = bx - ax, uy = by - ay, uz = bz - az;
        float vx = cx - ax, vy = cy - ay, vz = cz - az;
        float crx = uy * vz - uz * vy;
        float crz = ux * vy - uy * vx;
        float ccx = (ax + bx + cx) / 3.0f;
        float ccz = (az + bz + cz) / 3.0f;
        float r2c = ccx * ccx + ccz * ccz;
        if (r2c < 1e-3f) continue;
        float dot = crx * ccx + crz * ccz;
        if (dot > 0.0f) anyOutwardW = true;
    }
    ASSERT(anyOutwardW, "par_cone: triangle winding is outward");

    // Seam vertices (duplicate positions emitted by par_shapes for UV
    // unwrapping) must share the same normal on the lateral surface, so
    // the seam line is invisible under smooth shading. The apex is a
    // singularity — all (slices+1) apex copies share position (0,h,0)
    // and intentionally retain par_shapes' per-face normals to keep the
    // silhouette pointed.
    for (size_t i = 0; i < open.vertexCount(); ++i) {
        float xi = open.positions[i * 3 + 0];
        float zi = open.positions[i * 3 + 2];
        float ri = std::sqrt(xi * xi + zi * zi);
        if (ri < 1e-4f) continue;
        for (size_t j = i + 1; j < open.vertexCount(); ++j) {
            float dx = open.positions[j * 3 + 0] - xi;
            float dy = open.positions[j * 3 + 1] - open.positions[i * 3 + 1];
            float dz = open.positions[j * 3 + 2] - zi;
            if (dx*dx + dy*dy + dz*dz > 1e-8f) continue;
            float nd0 = open.normals[i * 3 + 0] - open.normals[j * 3 + 0];
            float nd1 = open.normals[i * 3 + 1] - open.normals[j * 3 + 1];
            float nd2 = open.normals[i * 3 + 2] - open.normals[j * 3 + 2];
            ASSERT(std::sqrt(nd0*nd0 + nd1*nd1 + nd2*nd2) < 1e-4f,
                   "par_cone: lateral seam normals match across duplicates");
        }
    }
}

TEST(par_disc) {
    auto mesh = bromesh::disc(1.5f, 16);
    ASSERT(!mesh.empty(), "par_disc: should be non-empty");
    ASSERT(mesh.vertexCount() > 0, "par_disc: should have vertices");
    // Disc lies in the XZ plane (Y is the normal — bro's Y-up convention).
    float maxAbsY = 0.0f, maxAbsX = 0.0f, maxAbsZ = 0.0f;
    for (size_t i = 0; i < mesh.vertexCount(); ++i) {
        maxAbsX = std::max(maxAbsX, std::fabs(mesh.positions[i * 3 + 0]));
        maxAbsY = std::max(maxAbsY, std::fabs(mesh.positions[i * 3 + 1]));
        maxAbsZ = std::max(maxAbsZ, std::fabs(mesh.positions[i * 3 + 2]));
    }
    ASSERT(maxAbsY < 1e-5f, "par_disc: every vertex sits in the XZ plane");
    ASSERT(std::fabs(maxAbsX - 1.5f) < 1e-4f, "par_disc: X extent matches radius");
    ASSERT(std::fabs(maxAbsZ - 1.5f) < 1e-4f, "par_disc: Z extent matches radius");
    if (mesh.hasNormals()) {
        for (size_t i = 0; i < mesh.vertexCount(); ++i) {
            ASSERT(mesh.normals[i * 3 + 1] > 0.99f,
                   "par_disc: normals point +Y");
        }
    }
}

TEST(par_rock) {
    auto mesh = bromesh::rock(1.0f, 42, 2);
    ASSERT(!mesh.empty(), "par_rock: should be non-empty");
    ASSERT(mesh.vertexCount() > 20, "par_rock: should have reasonable vertex count");
    ASSERT(mesh.hasNormals(), "par_rock: should have normals");

    // Different seeds should produce different shapes
    auto mesh2 = bromesh::rock(1.0f, 99, 2);
    ASSERT(!mesh2.empty(), "par_rock2: should be non-empty");
    // Vertices should differ
    bool differ = false;
    size_t checkCount = std::min(mesh.vertexCount(), mesh2.vertexCount());
    for (size_t i = 0; i < checkCount * 3 && !differ; ++i) {
        if (std::fabs(mesh.positions[i] - mesh2.positions[i]) > 0.001f) differ = true;
    }
    ASSERT(differ, "par_rock: different seeds should produce different shapes");
}

TEST(par_trefoil_knot) {
    auto mesh = bromesh::trefoilKnot(1.0f, 32, 8);
    ASSERT(!mesh.empty(), "par_trefoil: should be non-empty");
    ASSERT(mesh.vertexCount() > 50, "par_trefoil: should have many vertices");
    ASSERT(mesh.hasNormals(), "par_trefoil: should have normals");
}

TEST(par_klein_bottle) {
    auto mesh = bromesh::kleinBottle(16, 8);
    ASSERT(!mesh.empty(), "par_klein: should be non-empty");
    ASSERT(mesh.vertexCount() > 50, "par_klein: should have many vertices");
    ASSERT(mesh.hasNormals(), "par_klein: should have normals");
}

#endif // BROMESH_HAS_PAR_SHAPES

TEST(winding_order_box) {
    auto m = bromesh::box(2, 3, 4);
    ASSERT(checkOutwardWinding(m, 0, 0, 0), "box: all face normals should point outward");
}

TEST(winding_order_sphere) {
    auto m = bromesh::sphere(2.0f, 16, 12);
    ASSERT(checkOutwardWinding(m, 0, 0, 0), "sphere: all face normals should point outward");
}

TEST(winding_order_cylinder) {
    auto m = bromesh::cylinder(1.5f, 2.0f, 24);
    ASSERT(checkOutwardWinding(m, 0, 0, 0), "cylinder: all face normals should point outward");
}

TEST(winding_order_capsule) {
    auto m = bromesh::capsule(1.0f, 1.5f, 16, 8);
    ASSERT(checkOutwardWinding(m, 0, 0, 0), "capsule: all face normals should point outward");
}

TEST(winding_order_torus) {
    auto m = bromesh::torus(3.0f, 0.5f, 24, 12);
    ASSERT(checkTorusWinding(m, 3.0f), "torus: all face normals should point outward from tube");
}

TEST(winding_order_sphere_high_res) {
    auto m = bromesh::sphere(5.0f, 64, 48);
    ASSERT(checkOutwardWinding(m, 0, 0, 0), "sphere_high_res: all face normals should point outward");
}

TEST(winding_order_capsule_tall) {
    auto m = bromesh::capsule(0.5f, 4.0f, 32, 16);
    ASSERT(checkOutwardWinding(m, 0, 0, 0), "capsule_tall: all face normals should point outward");
}

