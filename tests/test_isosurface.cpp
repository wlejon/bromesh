#include "test_framework.h"
#include <cmath>
#include <vector>

static void fillSphereField(float* field, int N, float radius) {
    float c = (N - 1) * 0.5f;
    for (int z = 0; z < N; ++z)
        for (int y = 0; y < N; ++y)
            for (int x = 0; x < N; ++x) {
                float dx = x - c, dy = y - c, dz = z - c;
                field[z * N * N + y * N + x] = std::sqrt(dx*dx + dy*dy + dz*dz) - radius;
            }
}

TEST(marching_cubes_sphere) {
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 5.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = dist - radius;
            }
        }
    }

    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!mesh.empty(), "sphere mesh should be non-empty");
    ASSERT(mesh.hasNormals(), "sphere mesh should have normals");
    ASSERT(mesh.vertexCount() >= 100, "sphere should have at least 100 vertices");
    ASSERT(mesh.vertexCount() <= 10000, "sphere should have at most 10000 vertices");
    ASSERT(mesh.triangleCount() >= 30, "sphere should have at least 30 triangles");

    // Verify all vertices lie within 0.6 cellSize of analytic sphere surface
    float cellSize = 1.0f;
    bool allVertsNearSurface = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float dx = mesh.positions[v * 3 + 0] - cx;
        float dy = mesh.positions[v * 3 + 1] - cy;
        float dz = mesh.positions[v * 3 + 2] - cz;
        float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (std::fabs(dist - radius) >= 0.6f * cellSize) {
            allVertsNearSurface = false;
            break;
        }
    }
    ASSERT(allVertsNearSurface, "marching cubes: all vertices within 0.6 cellSize of sphere surface");

    // Verify mesh volume matches analytic sphere volume V = 4/3 * pi * R^3 within 5%
    float analyticVolume = (4.0f / 3.0f) * 3.14159265f * radius * radius * radius;
    float meshVolume = bromesh::computeVolume(mesh);
    ASSERT(std::fabs(meshVolume - analyticVolume) < 0.05f * analyticVolume,
           "marching cubes: volume matches analytic sphere within 5%");

    // Verify normals are unit length
    bool allUnit = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float nx = mesh.normals[v * 3 + 0];
        float ny = mesh.normals[v * 3 + 1];
        float nz = mesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.05f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "sphere normals should be unit length");
}

TEST(marching_cubes_all_inside) {
    // Standard SDF: f < iso means inside. An all-negative field is "all inside"
    // and with closeBoundary=true should produce a closed skin around the grid.
    const int N = 4;
    float field[N * N * N];
    for (int i = 0; i < N * N * N; ++i)
        field[i] = -1.0f; // all below iso = all inside

    // closeBoundary=true (default): pads with a super-iso ("outside") sentinel,
    // so an all-inside field produces a closed skin around the entire grid.
    auto closed = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!closed.empty(), "all-inside field with closeBoundary yields a closed skin");

    // closeBoundary=false: no boundary handling, so an all-inside field has no
    // zero crossings inside the grid and produces no triangles.
    auto open = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f, false);
    ASSERT(open.empty(), "all-inside field with closeBoundary=false stays empty");
}

TEST(marching_cubes_all_outside) {
    // Standard SDF: f >= iso means outside. An all-positive field is fully
    // outside; closeBoundary's "outside" sentinel matches, so no surface.
    const int N = 4;
    float field[N * N * N];
    for (int i = 0; i < N * N * N; ++i)
        field[i] = 1.0f; // all above iso = all outside
    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(mesh.empty(), "all-outside field should produce empty mesh");
}

TEST(surface_nets_sphere) {
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 5.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = dist - radius;
            }
        }
    }

    auto snMesh = bromesh::surfaceNets(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!snMesh.empty(), "surface nets sphere should be non-empty");
    ASSERT(snMesh.hasNormals(), "surface nets sphere should have normals");
    ASSERT(snMesh.vertexCount() >= 50, "surface nets sphere should have at least 50 vertices");
    ASSERT(snMesh.triangleCount() >= 50, "surface nets sphere should have at least 50 triangles");

    // Verify all vertices lie within 0.6 cellSize of analytic sphere surface
    float cellSize = 1.0f;
    bool allVertsNearSurface = true;
    for (size_t v = 0; v < snMesh.vertexCount(); ++v) {
        float dx = snMesh.positions[v * 3 + 0] - cx;
        float dy = snMesh.positions[v * 3 + 1] - cy;
        float dz = snMesh.positions[v * 3 + 2] - cz;
        float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (std::fabs(dist - radius) >= 0.6f * cellSize) {
            allVertsNearSurface = false;
            break;
        }
    }
    ASSERT(allVertsNearSurface, "surface nets: all vertices within 0.6 cellSize of sphere surface");

    // Verify mesh volume matches analytic sphere volume within 5%
    float analyticVolume = (4.0f / 3.0f) * 3.14159265f * radius * radius * radius;
    float meshVolume = bromesh::computeVolume(snMesh);
    ASSERT(std::fabs(meshVolume - analyticVolume) < 0.05f * analyticVolume,
           "surface nets: volume matches analytic sphere within 5%");

    // Now that marching cubes welds shared-edge vertices the two algorithms
    // produce comparable counts (surface nets: 1 vert per surface-containing
    // cell; welded MC: 1 vert per intersected cube edge). The original test
    // assumed the unwelded MC produced ~3x more verts. Today the two land
    // within a few percent of each other for a sphere — verify both produce
    // a same-order-of-magnitude mesh instead.
    auto mcMesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(mcMesh.vertexCount() > 0, "marching cubes baseline non-empty");
    double ratio = double(snMesh.vertexCount()) / double(mcMesh.vertexCount());
    ASSERT(ratio > 0.5 && ratio < 2.0,
           "surface nets and welded MC should be within 2x of each other");

    // Verify normals are unit length
    bool allUnit = true;
    for (size_t v = 0; v < snMesh.vertexCount(); ++v) {
        float nx = snMesh.normals[v * 3 + 0];
        float ny = snMesh.normals[v * 3 + 1];
        float nz = snMesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.05f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "surface nets sphere normals should be unit length");
}

TEST(dual_contour_sphere) {
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 5.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = dist - radius;
            }
        }
    }

    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!mesh.empty(), "dual contour sphere should be non-empty");
    ASSERT(mesh.hasNormals(), "dual contour sphere should have normals");
    ASSERT(mesh.vertexCount() >= 50, "dual contour sphere should have at least 50 vertices");
    ASSERT(mesh.triangleCount() >= 50, "dual contour sphere should have at least 50 triangles");

    // Verify all vertices lie within 0.6 cellSize of analytic sphere surface
    float cellSize = 1.0f;
    bool allVertsNearSurface = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float dx = mesh.positions[v * 3 + 0] - cx;
        float dy = mesh.positions[v * 3 + 1] - cy;
        float dz = mesh.positions[v * 3 + 2] - cz;
        float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (std::fabs(dist - radius) >= 0.6f * cellSize) {
            allVertsNearSurface = false;
            break;
        }
    }
    ASSERT(allVertsNearSurface, "dual contour: all vertices within 0.6 cellSize of sphere surface");

    // Verify mesh volume matches analytic sphere volume within 5%
    float analyticVolume = (4.0f / 3.0f) * 3.14159265f * radius * radius * radius;
    float meshVolume = bromesh::computeVolume(mesh);
    ASSERT(std::fabs(meshVolume - analyticVolume) < 0.05f * analyticVolume,
           "dual contour: volume matches analytic sphere within 5%");

    // Verify normals are unit length
    bool allUnit = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float nx = mesh.normals[v * 3 + 0];
        float ny = mesh.normals[v * 3 + 1];
        float nz = mesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.05f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "dual contour sphere normals should be unit length");
}

TEST(dual_contour_box_field) {
    // Box SDF: max of 6 plane distances
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float halfExtent = 4.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = std::fabs(x - cx) - halfExtent;
                float dy = std::fabs(y - cy) - halfExtent;
                float dz = std::fabs(z - cz) - halfExtent;
                // Standard SDF: max of 3 plane distances; negative inside, positive outside.
                float maxDist = dx;
                if (dy > maxDist) maxDist = dy;
                if (dz > maxDist) maxDist = dz;
                field[z * N * N + y * N + x] = maxDist;
            }
        }
    }

    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!mesh.empty(), "dual contour box should be non-empty");
    ASSERT(mesh.hasNormals(), "dual contour box should have normals");
    ASSERT(mesh.vertexCount() >= 8, "dual contour box should have at least 8 vertices");
    ASSERT(mesh.triangleCount() >= 12, "dual contour box should have at least 12 triangles");

    // Verify mesh volume matches analytic box volume (2*4)^3 = 512 within 5%
    float analyticBoxVolume = (2.0f * halfExtent) * (2.0f * halfExtent) * (2.0f * halfExtent);
    float meshVolume = bromesh::computeVolume(mesh);
    ASSERT(std::fabs(meshVolume - analyticBoxVolume) < 0.05f * analyticBoxVolume,
           "dual contour box: volume matches analytic box (512.0) within 5%");

    // Assert that vertices lie near the bounding planes |x-cx| ~= 4.0, |y-cy| ~= 4.0, |z-cz| ~= 4.0 (within 0.5 cell size)
    float cellSize = 1.0f;
    bool allVertsNearBoxFaces = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float dx = std::fabs(mesh.positions[v * 3 + 0] - cx);
        float dy = std::fabs(mesh.positions[v * 3 + 1] - cy);
        float dz = std::fabs(mesh.positions[v * 3 + 2] - cz);
        float maxD = std::max(dx, std::max(dy, dz));
        if (std::fabs(maxD - halfExtent) > 0.5f * cellSize ||
            dx > halfExtent + 0.5f * cellSize ||
            dy > halfExtent + 0.5f * cellSize ||
            dz > halfExtent + 0.5f * cellSize) {
            allVertsNearBoxFaces = false;
            break;
        }
    }
    ASSERT(allVertsNearBoxFaces, "dual contour box: vertices lie near bounding planes within 0.5 cellSize");
}

TEST(transvoxel_uniform_lod) {
    // 17x17x17 sphere field, lod=0, all neighborLods=-1 (no neighbors).
    // Should produce same result as marching cubes.
    const int N = 17;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 6.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = dist - radius;
            }
        }
    }

    int neighborLods[6] = { -1, -1, -1, -1, -1, -1 };
    auto tvMesh = bromesh::transvoxel(field, N, 0, neighborLods, 0.0f, 1.0f);
    // Transvoxel doesn't yet implement boundary closing, so compare against
    // the open-boundary form of marching cubes for an apples-to-apples count.
    auto mcMesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f, false);

    ASSERT(!tvMesh.empty(), "transvoxel uniform lod should produce non-empty mesh");
    ASSERT(tvMesh.hasNormals(), "transvoxel uniform lod should have normals");
    // With lod=0 and no neighbors, transvoxel should produce the same geometry as marching cubes
    ASSERT(tvMesh.vertexCount() == mcMesh.vertexCount(),
           "transvoxel lod=0 no-neighbors should match marching cubes vertex count");
    ASSERT(tvMesh.triangleCount() == mcMesh.triangleCount(),
           "transvoxel lod=0 no-neighbors should match marching cubes triangle count");

    // Verify all vertices lie within 0.6 cellSize of analytic sphere surface
    float cellSize = 1.0f;
    bool allVertsNearSurface = true;
    for (size_t v = 0; v < tvMesh.vertexCount(); ++v) {
        float dx = tvMesh.positions[v * 3 + 0] - cx;
        float dy = tvMesh.positions[v * 3 + 1] - cy;
        float dz = tvMesh.positions[v * 3 + 2] - cz;
        float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (std::fabs(dist - radius) >= 0.6f * cellSize) {
            allVertsNearSurface = false;
            break;
        }
    }
    ASSERT(allVertsNearSurface, "transvoxel: all vertices within 0.6 cellSize of sphere surface");

    // Verify mesh volume matches analytic sphere volume V = 4/3 * pi * R^3 within 5%
    float analyticVolume = (4.0f / 3.0f) * 3.14159265f * radius * radius * radius;
    float meshVolume = bromesh::computeVolume(tvMesh);
    ASSERT(std::fabs(meshVolume - analyticVolume) < 0.05f * analyticVolume,
           "transvoxel: volume matches analytic sphere within 5%");
}

TEST(transvoxel_with_transition) {
    // 17x17x17 sphere field, lod=0, one neighbor at lod=1: an LOD transition
    // across +X. No transition cells are emitted; any vertex on that face is
    // snapped to the neighbour's grid. (This sphere does not reach x = 16, so
    // transvoxel_snaps_to_coarser_neighbor below is the one that sees
    // boundary vertices.)
    const int N = 17;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 6.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = dist - radius;
            }
        }
    }

    // +X neighbor has coarser LOD (lod=1), rest are -1
    int neighborLods[6] = { 1, -1, -1, -1, -1, -1 };
    auto tvMesh = bromesh::transvoxel(field, N, 0, neighborLods, 0.0f, 1.0f);

    ASSERT(!tvMesh.empty(), "transvoxel with transition should produce non-empty mesh");
    ASSERT(tvMesh.hasNormals(), "transvoxel with transition should have normals");
    ASSERT(tvMesh.vertexCount() >= 100, "transvoxel with transition should have reasonable vertex count");
    ASSERT(tvMesh.triangleCount() >= 30, "transvoxel with transition should have reasonable triangle count");

    // Verify that vertices on the +X boundary (x near maxCoord=16) are snapped to
    // the neighbor's grid (stride=2, so y and z should be multiples of 2).
    float maxCoord = (float)(N - 1);
    float neighborStep = 2.0f; // 1 << neighborLod(1) * cellSize(1.0)
    float tolerance = 0.01f;
    bool allSnapped = true;
    int boundaryVertCount = 0;
    for (size_t v = 0; v < tvMesh.vertexCount(); ++v) {
        float vx = tvMesh.positions[v * 3 + 0];
        if (std::fabs(vx - maxCoord) < tolerance) {
            boundaryVertCount++;
            float vy = tvMesh.positions[v * 3 + 1];
            float vz = tvMesh.positions[v * 3 + 2];
            // Check that y and z are snapped to multiples of neighborStep
            float snapY = std::round(vy / neighborStep) * neighborStep;
            float snapZ = std::round(vz / neighborStep) * neighborStep;
            if (std::fabs(vy - snapY) > tolerance || std::fabs(vz - snapZ) > tolerance) {
                allSnapped = false;
                break;
            }
        }
    }
    // Any boundary vertices there are must be snapped (see the comment above:
    // this field has none on x = 16).
    if (boundaryVertCount > 0) {
        ASSERT(allSnapped, "transvoxel +X boundary vertices should be snapped to neighbor grid");
    }
}

namespace {
// A tilted plane that crosses every face of a 17^3 chunk off the lod-1 grid.
void fillTiltedPlane(float* field, int n) {
    for (int z = 0; z < n; ++z)
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
                field[z * n * n + y * n + x] = y + 0.31f * z + 0.17f * x - 7.3f;
}

// Vertices on x = maxCoord, and how many of them sit off the `step` grid in y/z.
void countPlusXBoundary(const bromesh::MeshData& m, float maxCoord, float step, int& onFace, int& offGrid) {
    onFace = 0; offGrid = 0;
    for (size_t v = 0; v < m.vertexCount(); ++v) {
        if (std::fabs(m.positions[v * 3] - maxCoord) > 1e-3f) continue;
        ++onFace;
        const float y = m.positions[v * 3 + 1], z = m.positions[v * 3 + 2];
        if (std::fabs(y - std::round(y / step) * step) > 1e-3f ||
            std::fabs(z - std::round(z / step) * step) > 1e-3f) {
            ++offGrid;
        }
    }
}
} // namespace

TEST(transvoxel_snaps_to_coarser_neighbor) {
    const int N = 17;
    std::vector<float> field(N * N * N);
    fillTiltedPlane(field.data(), N);

    // A same-LOD neighbour leaves the face alone: the plane's crossings are
    // off the lod-1 grid there.
    int sameLods[6] = {0, -1, -1, -1, -1, -1};
    auto same = bromesh::transvoxel(field.data(), N, 0, sameLods, 0.0f, 1.0f);
    int onFace = 0, offGrid = 0;
    countPlusXBoundary(same, 16.0f, 2.0f, onFace, offGrid);
    ASSERT(onFace > 0, "transvoxel: the plane crosses the +X face");
    ASSERT(offGrid > 0, "transvoxel: a same-LOD neighbour does not snap");

    // A coarser neighbour snaps every +X vertex onto its grid.
    int coarseLods[6] = {1, -1, -1, -1, -1, -1};
    auto snapped = bromesh::transvoxel(field.data(), N, 0, coarseLods, 0.0f, 1.0f);
    countPlusXBoundary(snapped, 16.0f, 2.0f, onFace, offGrid);
    ASSERT(onFace > 0, "transvoxel: +X face vertices survive the snap");
    ASSERT(offGrid == 0, "transvoxel: every +X vertex is on the lod-1 grid");
}

TEST(transvoxel_rejects_out_of_range_lod) {
    const int N = 17;
    std::vector<float> field(N * N * N);
    fillTiltedPlane(field.data(), N);
    int none[6] = {-1, -1, -1, -1, -1, -1};
    // 1 << lod is undefined past the int's bits; these used to reach it.
    for (int lod : {-1, -100, 31, 32, 1000}) {
        ASSERT(bromesh::transvoxel(field.data(), N, lod, none).empty(), "transvoxel: out-of-range lod is empty");
    }
    // 2^5 = 32 > N - 1: no cell fits.
    ASSERT(bromesh::transvoxel(field.data(), N, 5, none).empty(), "transvoxel: a lod coarser than the chunk is empty");
    ASSERT(!bromesh::transvoxel(field.data(), N, 4, none).empty(), "transvoxel: lod 4 fits a 17^3 chunk");
    // A neighbour level past the int's bits snaps as the maximum instead.
    int huge[6] = {1000, -1, -1, -1, -1, -1};
    ASSERT(!bromesh::transvoxel(field.data(), N, 0, huge).empty(), "transvoxel: a huge neighbour lod is survivable");
}

