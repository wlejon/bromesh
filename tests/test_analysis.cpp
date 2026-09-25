#include "test_framework.h"
#include <cmath>

TEST(bbox_basics) {
    bromath::AABB3 b;
    b.min = {-1, -2, -3};
    b.max = { 1,  2,  3};
    auto c = bromath::acenter(b);
    auto e = bromath::aextent(b);
    auto h = bromath::ahalfExtent(b);
    ASSERT(std::fabs(c.x) < 0.001f, "center X");
    ASSERT(std::fabs(c.y) < 0.001f, "center Y");
    ASSERT(std::fabs(c.z) < 0.001f, "center Z");
    ASSERT(std::fabs(e.x - 2.0f) < 0.001f, "extent X (full size)");
    ASSERT(std::fabs(e.y - 4.0f) < 0.001f, "extent Y (full size)");
    ASSERT(std::fabs(e.z - 6.0f) < 0.001f, "extent Z (full size)");
    ASSERT(std::fabs(h.x - 1.0f) < 0.001f, "halfExtent X");
    ASSERT(std::fabs(h.y - 2.0f) < 0.001f, "halfExtent Y");
    ASSERT(std::fabs(h.z - 3.0f) < 0.001f, "halfExtent Z");
}

TEST(compute_bbox_box) {
    auto b = bromesh::box(1, 1, 1);
    auto bbox = bromesh::computeBBox(b);
    ASSERT(std::fabs(bbox.min.x - (-1.0f)) < 0.001f, "bbox min x == -1");
    ASSERT(std::fabs(bbox.min.y - (-1.0f)) < 0.001f, "bbox min y == -1");
    ASSERT(std::fabs(bbox.min.z - (-1.0f)) < 0.001f, "bbox min z == -1");
    ASSERT(std::fabs(bbox.max.x - 1.0f) < 0.001f, "bbox max x == 1");
    ASSERT(std::fabs(bbox.max.y - 1.0f) < 0.001f, "bbox max y == 1");
    ASSERT(std::fabs(bbox.max.z - 1.0f) < 0.001f, "bbox max z == 1");
}

TEST(is_manifold_box) {
    auto b = bromesh::box(1, 1, 1);
    ASSERT(bromesh::isManifold(b), "box should be manifold");
}

TEST(compute_volume_box) {
    auto b = bromesh::box(1, 1, 1);
    float vol = bromesh::computeVolume(b);
    ASSERT(std::fabs(vol - 8.0f) < 0.1f, "box(1,1,1) volume should be ~8.0");
}

TEST(surface_sample_basic) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    bromesh::computeNormals(mesh);

    auto samples = bromesh::sampleSurface(mesh, 100, 42);
    ASSERT(samples.vertexCount() == 100, "sample: 100 points");
    ASSERT(samples.hasNormals(), "sample: has normals");
    ASSERT(samples.indices.empty(), "sample: no indices (point cloud)");
}

TEST(surface_sample_on_sphere) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    bromesh::computeNormals(mesh);

    auto samples = bromesh::sampleSurface(mesh, 200, 42);

    // All points should be approximately on the sphere surface (radius ~= 2)
    bool onSurface = true;
    for (size_t v = 0; v < samples.vertexCount(); ++v) {
        float x = samples.positions[v*3+0];
        float y = samples.positions[v*3+1];
        float z = samples.positions[v*3+2];
        float r = std::sqrt(x*x + y*y + z*z);
        if (std::fabs(r - 2.0f) > 0.15f) { onSurface = false; break; }
    }
    ASSERT(onSurface, "sample_sphere: all points near radius 2");
}

TEST(surface_sample_with_uvs) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    // box already has UVs

    auto samples = bromesh::sampleSurface(mesh, 50, 42);
    ASSERT(samples.hasUVs(), "sample_uvs: should have UVs when source does");
}

TEST(surface_sample_deterministic) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);

    auto s1 = bromesh::sampleSurface(mesh, 50, 123);
    auto s2 = bromesh::sampleSurface(mesh, 50, 123);

    ASSERT(s1.positions == s2.positions, "sample_deterministic: same seed = same result");
}

TEST(surface_area_box) {
    // Box with half-extents 1 -> side length 2 -> 6 faces * 4 = 24
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    float area = bromesh::computeSurfaceArea(mesh);
    ASSERT(std::fabs(area - 24.0f) < 0.1f, "surface_area: unit box ~= 24");
}

TEST(triangle_areas_count) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto areas = bromesh::computeTriangleAreas(mesh);
    ASSERT(areas.size() == mesh.triangleCount(), "tri_areas: one per triangle");

    // All areas should be positive
    bool allPositive = true;
    for (float a : areas) {
        if (a <= 0.0f) { allPositive = false; break; }
    }
    ASSERT(allPositive, "tri_areas: all positive");
}

TEST(raycast_sphere_center) {
    auto mesh = bromesh::sphere(2.0f, 16, 12);
    bromesh::computeNormals(mesh);

    // Shoot ray from outside along -Z toward center
    float origin[3] = {0, 0, 5};
    float dir[3] = {0, 0, -1};
    auto hit = bromesh::raycast(mesh, origin, dir);

    ASSERT(hit.hit, "raycast_sphere: should hit");
    ASSERT(std::fabs(hit.position[2] - 2.0f) < 0.3f, "raycast_sphere: hit near +Z pole");
    ASSERT(hit.distance > 0, "raycast_sphere: positive distance");
    // Normal should be roughly aligned with +Z or -Z (face normal may point inward depending on winding)
    ASSERT(std::fabs(hit.normal[2]) > 0.5f, "raycast_sphere: normal has strong Z component");
}

TEST(raycast_miss) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);

    // Shoot ray that misses entirely
    float origin[3] = {10, 10, 10};
    float dir[3] = {1, 0, 0}; // pointing away
    auto hit = bromesh::raycast(mesh, origin, dir);

    ASSERT(!hit.hit, "raycast_miss: should not hit");
}

TEST(raycast_max_distance) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);

    float origin[3] = {0, 0, 10};
    float dir[3] = {0, 0, -1};

    // Max distance too short
    auto hit = bromesh::raycast(mesh, origin, dir, 5.0f);
    ASSERT(!hit.hit, "raycast_maxdist: too far, should miss");

    // Max distance sufficient
    auto hit2 = bromesh::raycast(mesh, origin, dir, 20.0f);
    ASSERT(hit2.hit, "raycast_maxdist: close enough, should hit");
}

TEST(raycast_all_through_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);

    // Ray through center of box should hit 2 faces (entry + exit)
    float origin[3] = {0, 0, 5};
    float dir[3] = {0, 0, -1};
    auto hits = bromesh::raycastAll(mesh, origin, dir);

    ASSERT(hits.size() >= 2, "raycast_all: at least 2 hits through box");
    // Should be sorted by distance
    if (hits.size() >= 2) {
        ASSERT(hits[0].distance <= hits[1].distance, "raycast_all: sorted by distance");
    }
}

TEST(raycast_test_fast) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);

    float origin[3] = {0, 0, 5};
    float dir[3] = {0, 0, -1};
    ASSERT(bromesh::raycastTest(mesh, origin, dir), "raycast_test: should hit box");

    float dir2[3] = {0, 0, 1}; // away
    ASSERT(!bromesh::raycastTest(mesh, origin, dir2), "raycast_test: should miss");
}

TEST(bvh_build_empty) {
    bromesh::MeshData empty;
    auto bvh = bromesh::MeshBVH::build(empty);
    ASSERT(bvh.empty(), "bvh_empty: empty mesh → empty BVH");
    ASSERT(bvh.nodeCount() == 0, "bvh_empty: 0 nodes");
}

TEST(bvh_build_box) {
    auto mesh = bromesh::box(1.0f, 2.0f, 3.0f);
    auto bvh = bromesh::MeshBVH::build(mesh);
    ASSERT(!bvh.empty(), "bvh_box: non-empty");
    ASSERT(bvh.triangleCount() == mesh.triangleCount(), "bvh_box: indexes every tri");
    auto bb = bvh.bounds();
    ASSERT(std::fabs(bb.min.x - -1.0f) < 1e-4f, "bvh_box: bounds minX");
    ASSERT(std::fabs(bb.max.y -  2.0f) < 1e-4f, "bvh_box: bounds maxY");
    ASSERT(std::fabs(bb.max.z -  3.0f) < 1e-4f, "bvh_box: bounds maxZ");
}

TEST(bvh_raycast_matches_brute_force) {
    // Use a dense sphere so the BVH actually has multiple levels.
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    bromesh::computeNormals(mesh);
    auto bvh = bromesh::MeshBVH::build(mesh);

    // Shoot a handful of rays and cross-check with the brute-force raycast.
    struct Ray { float o[3], d[3]; };
    Ray rays[] = {
        {{0, 0,  5}, { 0,  0, -1}},
        {{5, 0,  0}, {-1,  0,  0}},
        {{0, 5,  0}, { 0, -1,  0}},
        {{3, 3,  3}, {-1, -1, -1}},
        {{0, 0, 10}, { 0,  0, -1}},
    };
    for (const auto& r : rays) {
        auto a = bromesh::raycast(mesh, r.o, r.d);
        auto b = bvh.raycast(mesh, r.o, r.d);
        ASSERT(a.hit == b.hit, "bvh_raycast_match: hit flag agrees");
        if (a.hit && b.hit) {
            ASSERT(std::fabs(a.distance - b.distance) < 1e-3f,
                   "bvh_raycast_match: distance agrees");
            ASSERT(a.triangleIndex == b.triangleIndex ||
                   std::fabs(a.distance - b.distance) < 1e-3f,
                   "bvh_raycast_match: same or equidistant triangle");
        }
    }
}

TEST(bvh_raycast_miss) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto bvh = bromesh::MeshBVH::build(mesh);

    float origin[3] = {10, 10, 10};
    float dir[3]    = {1, 0, 0};
    auto hit = bvh.raycast(mesh, origin, dir);
    ASSERT(!hit.hit, "bvh_miss: ray pointing away");
}

TEST(bvh_raycast_max_distance) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto bvh = bromesh::MeshBVH::build(mesh);

    float origin[3] = {0, 0, 10};
    float dir[3]    = {0, 0, -1};
    auto tooShort = bvh.raycast(mesh, origin, dir, 5.0f);
    ASSERT(!tooShort.hit, "bvh_maxdist: short ray misses");
    auto ok = bvh.raycast(mesh, origin, dir, 20.0f);
    ASSERT(ok.hit, "bvh_maxdist: long ray hits");
}

TEST(bvh_raycast_axis_aligned_ray) {
    // Pure axis-aligned rays stress the slab test (one or more invDir entries
    // would be infinite). Make sure it handles them.
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto bvh = bromesh::MeshBVH::build(mesh);

    float origin[3] = {0, 0, 5};
    float dir[3]    = {0, 0, -1};
    auto hit = bvh.raycast(mesh, origin, dir);
    ASSERT(hit.hit, "bvh_axis: axis-aligned ray hits box");
    ASSERT(std::fabs(hit.position[2] - 1.0f) < 1e-3f, "bvh_axis: hits front face Z=1");
}

TEST(bvh_raycast_test_fast) {
    auto mesh = bromesh::sphere(1.5f, 24, 16);
    auto bvh = bromesh::MeshBVH::build(mesh);

    float origin[3] = {0, 0, 5};
    float dir[3]    = {0, 0, -1};
    ASSERT(bvh.raycastTest(mesh, origin, dir), "bvh_test: should hit");
    float dir2[3]   = {0, 0, 1};
    ASSERT(!bvh.raycastTest(mesh, origin, dir2), "bvh_test: should miss");
}

TEST(bvh_handles_degenerate_triangles) {
    // Three collinear vertices → degenerate triangle. BVH should still build
    // without infinite recursion and raycast should not hit.
    bromesh::MeshData m;
    m.positions = {
        0, 0, 0,
        1, 0, 0,
        2, 0, 0,
        // A real triangle to make raycast meaningful.
        0, 0, -1,
        1, 0, -1,
        0, 1, -1,
    };
    m.indices = { 0, 1, 2,  3, 4, 5 };
    auto bvh = bromesh::MeshBVH::build(m);
    ASSERT(!bvh.empty(), "bvh_degen: still builds");

    float o[3] = {0.2f, 0.2f, 5};
    float d[3] = {0, 0, -1};
    auto hit = bvh.raycast(m, o, d);
    ASSERT(hit.hit, "bvh_degen: hits the real triangle");
    ASSERT(hit.triangleIndex == 1, "bvh_degen: hits triangle 1, not the degenerate one");
}

TEST(closest_point_on_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    // Point above the box
    float point[3] = {0, 5, 0};
    auto cp = bromesh::closestPoint(mesh, point);

    ASSERT(cp.hit, "closest_point: should find a point");
    ASSERT(std::fabs(cp.position[1] - 1.0f) < 0.01f, "closest_point: on top face Y=1");
    ASSERT(cp.distance > 3.5f && cp.distance < 4.5f, "closest_point: distance ~4");
}

TEST(closest_point_inside) {
    auto mesh = bromesh::sphere(2.0f, 16, 12);

    // Point at center
    float point[3] = {0, 0, 0};
    auto cp = bromesh::closestPoint(mesh, point);

    ASSERT(cp.hit, "closest_inside: should find a point");
    // Distance should be approximately the radius
    ASSERT(std::fabs(cp.distance - 2.0f) < 0.3f, "closest_inside: distance ~= radius");
}

TEST(self_intersect_clean_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);

    ASSERT(!bromesh::hasSelfIntersections(mesh), "clean_box: no self-intersections");
    auto pairs = bromesh::findSelfIntersections(mesh);
    ASSERT(pairs.empty(), "clean_box: no intersection pairs");
}

TEST(self_intersect_clean_sphere) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    ASSERT(!bromesh::hasSelfIntersections(mesh), "clean_sphere: no self-intersections");
}

TEST(self_intersect_created) {
    // Create a self-intersecting mesh by merging two overlapping boxes
    auto box1 = bromesh::box(1.0f, 1.0f, 1.0f);
    auto box2 = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::translateMesh(box2, 0.5f, 0.5f, 0.5f);

    auto merged = bromesh::mergeMeshes({box1, box2});
    ASSERT(bromesh::hasSelfIntersections(merged),
           "overlapping_boxes: should have self-intersections");

    auto pairs = bromesh::findSelfIntersections(merged);
    ASSERT(!pairs.empty(), "overlapping_boxes: should find intersection pairs");
}

TEST(meshes_intersect_overlapping) {
    auto box1 = bromesh::box(1.0f, 1.0f, 1.0f);
    auto box2 = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::translateMesh(box2, 0.5f, 0.5f, 0.5f);

    ASSERT(bromesh::meshesIntersect(box1, box2), "overlap: boxes should intersect");
}

TEST(meshes_intersect_separated) {
    auto box1 = bromesh::box(1.0f, 1.0f, 1.0f);
    auto box2 = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::translateMesh(box2, 10.0f, 0.0f, 0.0f);

    ASSERT(!bromesh::meshesIntersect(box1, box2), "separated: boxes should not intersect");
}

TEST(bvh_closest_point_matches_brute_force) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    bromesh::computeNormals(mesh);
    ASSERT(mesh.triangleCount() >= 500, "mesh has 500+ triangles");

    bromesh::MeshBVH bvh = bromesh::MeshBVH::build(mesh);
    ASSERT(!bvh.empty(), "BVH built");

    const size_t triCount = mesh.triangleCount();
    for (int i = 0; i < 50; ++i) {
        // Sample interior of a triangle to ensure a unique closest triangle
        uint32_t tri = static_cast<uint32_t>((i * 31 + 7) % triCount);
        uint32_t i0 = mesh.indices[tri * 3 + 0];
        uint32_t i1 = mesh.indices[tri * 3 + 1];
        uint32_t i2 = mesh.indices[tri * 3 + 2];

        float p0[3] = { mesh.positions[i0*3+0], mesh.positions[i0*3+1], mesh.positions[i0*3+2] };
        float p1[3] = { mesh.positions[i1*3+0], mesh.positions[i1*3+1], mesh.positions[i1*3+2] };
        float p2[3] = { mesh.positions[i2*3+0], mesh.positions[i2*3+1], mesh.positions[i2*3+2] };

        // Face normal
        float e1[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
        float e2[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
        float fn[3] = {
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0]
        };
        float fnLen = std::sqrt(fn[0]*fn[0] + fn[1]*fn[1] + fn[2]*fn[2]);
        if (fnLen > 1e-8f) { fn[0] /= fnLen; fn[1] /= fnLen; fn[2] /= fnLen; }

        // Barycentric weights in triangle interior, well away from edges
        float u = 0.25f, v = 0.35f, w = 0.40f;
        float surf[3] = {
            u * p0[0] + v * p1[0] + w * p2[0],
            u * p0[1] + v * p1[1] + w * p2[1],
            u * p0[2] + v * p1[2] + w * p2[2]
        };

        // Offset along face normal: inside (-0.02), near (0.01), outside (0.5)
        float d = (i % 3 == 0) ? -0.02f : ((i % 3 == 1) ? 0.01f : 0.5f);
        float pt[3] = {
            surf[0] + fn[0] * d,
            surf[1] + fn[1] * d,
            surf[2] + fn[2] * d
        };

        bromesh::RayHit brute = bromesh::closestPoint(mesh, pt);
        bromesh::RayHit fast = bvh.closestPoint(mesh, pt);

        ASSERT(brute.hit == fast.hit, "hit status matches");
        ASSERT(std::fabs(brute.distance - fast.distance) < 1e-4f, "distance matches");
        ASSERT(std::fabs(brute.position[0] - fast.position[0]) < 1e-4f, "position X matches");
        ASSERT(std::fabs(brute.position[1] - fast.position[1]) < 1e-4f, "position Y matches");
        ASSERT(std::fabs(brute.position[2] - fast.position[2]) < 1e-4f, "position Z matches");
        ASSERT(std::fabs(brute.normal[0] - fast.normal[0]) < 1e-4f, "normal X matches");
        ASSERT(std::fabs(brute.normal[1] - fast.normal[1]) < 1e-4f, "normal Y matches");
        ASSERT(std::fabs(brute.normal[2] - fast.normal[2]) < 1e-4f, "normal Z matches");
        ASSERT(brute.triangleIndex == fast.triangleIndex, "triangleIndex matches");
        ASSERT(std::fabs(brute.baryU - fast.baryU) < 1e-4f, "baryU matches");
        ASSERT(std::fabs(brute.baryV - fast.baryV) < 1e-4f, "baryV matches");
        ASSERT(std::fabs(brute.baryW - fast.baryW) < 1e-4f, "baryW matches");
    }
}

static bromesh::MeshData pointCloud(const std::vector<float>& xyz) {
    bromesh::MeshData m;
    m.positions = xyz;
    return m;
}

static bool hullIsOutward(const bromesh::MeshData& h) {
    float cx = 0, cy = 0, cz = 0;
    size_t n = h.vertexCount();
    for (size_t i = 0; i < n; ++i) { cx += h.positions[i * 3]; cy += h.positions[i * 3 + 1]; cz += h.positions[i * 3 + 2]; }
    cx /= n; cy /= n; cz /= n;
    for (size_t t = 0; t < h.triangleCount(); ++t) {
        const float* a = &h.positions[h.indices[t * 3] * 3];
        const float* b = &h.positions[h.indices[t * 3 + 1] * 3];
        const float* c = &h.positions[h.indices[t * 3 + 2] * 3];
        float ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
        float vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
        float nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
        if (nx * (a[0] - cx) + ny * (a[1] - cy) + nz * (a[2] - cz) <= 0) return false;
    }
    return true;
}

TEST(convex_hull_chamfered_box_point_cloud) {
    const float hx = 1.0f, hy = 0.6f, hz = 0.8f, c = 0.15f;
    std::vector<float> pts;
    for (float sx : {-1.f, 1.f}) for (float sy : {-1.f, 1.f}) for (float sz : {-1.f, 1.f}) {
        pts.insert(pts.end(), {sx * (hx - c), sy * hy, sz * (hz - c)});
        pts.insert(pts.end(), {sx * hx, sy * (hy - c), sz * (hz - c)});
        pts.insert(pts.end(), {sx * (hx - c), sy * (hy - c), sz * hz});
    }
    pts.insert(pts.end(), {0.1f, 0.2f, -0.1f, -0.3f, 0.0f, 0.4f});
    auto h = bromesh::convexHull(pointCloud(pts));
    ASSERT(h.triangleCount() == 44, "24 extreme points on a chamfered box hull into 44 triangles");
    ASSERT(h.vertexCount() == h.triangleCount() * 3, "hull is flat-shaded, three vertices per triangle");
    ASSERT(h.normals.size() == h.positions.size(), "hull carries normals");
    ASSERT(hullIsOutward(h), "every hull triangle winds outward");
    float box = 8 * hx * hy * hz;
    float removed = 4 * c * c * ((hx - c) + (hy - c) + (hz - c)) + 20 * c * c * c / 3;
    float vol = bromesh::computeVolume(h);
    ASSERT(std::fabs(vol - (box - removed)) < 1e-3f, "chamfered box hull volume is exact");
}

TEST(convex_hull_ignores_interior_and_topology) {
    auto sphere = bromesh::sphere(1.0f, 24, 16);
    auto h = bromesh::convexHull(sphere);
    ASSERT(!h.empty(), "sphere hull is non-empty");
    ASSERT(hullIsOutward(h), "sphere hull winds outward");
    float vs = bromesh::computeVolume(h);
    ASSERT(vs > 3.8f && vs < 4.19f, "sphere hull volume approaches 4/3 pi");
    auto b = bromesh::box(1, 2, 3);
    auto hb = bromesh::convexHull(b);
    ASSERT(std::fabs(bromesh::computeVolume(hb) - 48.0f) < 0.01f, "box hull keeps the box volume");
}

TEST(convex_hull_degenerate_inputs) {
    ASSERT(bromesh::convexHull(pointCloud({0, 0, 0, 1, 0, 0, 0, 1, 0})).empty(), "three points have no hull");
    ASSERT(bromesh::convexHull(pointCloud({0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0.5f, 0.5f, 0})).empty(), "coplanar points have no hull");
    ASSERT(bromesh::convexHull(pointCloud({0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3})).empty(), "collinear points have no hull");
    ASSERT(bromesh::convexHull(bromesh::MeshData{}).empty(), "empty input has no hull");
}


