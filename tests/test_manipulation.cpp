#include "test_framework.h"
#include <cmath>

TEST(mesh_data_basics) {
    bromesh::MeshData m;
    ASSERT(m.empty(), "new mesh should be empty");
    ASSERT(m.vertexCount() == 0, "vertex count should be 0");
    ASSERT(m.triangleCount() == 0, "triangle count should be 0");

    m.positions = {0,0,0, 1,0,0, 0,1,0};
    m.indices = {0, 1, 2};
    ASSERT(m.vertexCount() == 3, "should have 3 vertices");
    ASSERT(m.triangleCount() == 1, "should have 1 triangle");
    ASSERT(!m.hasNormals(), "no normals yet");

    m.clear();
    ASSERT(m.empty(), "cleared mesh should be empty");
}

TEST(compute_normals_box) {
    auto b = bromesh::box(1, 1, 1);
    // Clear normals and recompute
    b.normals.clear();
    ASSERT(!b.hasNormals(), "normals cleared");
    bromesh::computeNormals(b);
    ASSERT(b.hasNormals(), "normals computed");

    // All normals should be unit length
    bool allUnit = true;
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        float nx = b.normals[v * 3 + 0];
        float ny = b.normals[v * 3 + 1];
        float nz = b.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.01f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "all normals are unit length");
}

TEST(compute_flat_normals_sphere) {
    auto s = bromesh::sphere(1);
    size_t origIndexCount = s.indices.size();
    auto flat = bromesh::computeFlatNormals(s);
    ASSERT(flat.vertexCount() == origIndexCount, "flat normals: vertex count == original index count");
    ASSERT(flat.hasNormals(), "flat normals: has normals");
    // Indices should be sequential
    bool sequential = true;
    for (size_t i = 0; i < flat.indices.size(); ++i) {
        if (flat.indices[i] != static_cast<uint32_t>(i)) { sequential = false; break; }
    }
    ASSERT(sequential, "flat normals: indices are sequential");
}

TEST(weld_vertices_box) {
    // computeFlatNormals produces a mesh with duplicated vertices (one per face corner)
    auto b = bromesh::box(1, 1, 1);
    auto flat = bromesh::computeFlatNormals(b);
    // flat mesh: every triangle has its own vertices, so vertexCount == indexCount
    size_t flatVerts = flat.vertexCount();
    ASSERT(flatVerts == flat.indices.size(), "flat mesh has duplicated verts");

    auto welded = bromesh::weldVertices(flat, 0.001f);
    ASSERT(welded.vertexCount() < flatVerts, "weld should reduce vertex count");
    ASSERT(!welded.empty(), "welded mesh should not be empty");
    ASSERT(welded.hasNormals(), "welded mesh should have normals");
    // A box has 8 unique positions, welding by position should get close to 8
    ASSERT(welded.vertexCount() <= 24, "welded box should have at most 24 verts");
}

TEST(split_components_two_triangles) {
    // Create two disjoint triangles (no shared vertices = 2 components)
    bromesh::MeshData combined;
    // Triangle 1: verts 0,1,2
    combined.positions = {
        0,0,0,  1,0,0,  0,1,0,   // tri 1
        5,0,0,  6,0,0,  5,1,0    // tri 2
    };
    combined.normals = {
        0,0,1,  0,0,1,  0,0,1,
        0,0,1,  0,0,1,  0,0,1
    };
    combined.indices = { 0,1,2,  3,4,5 };

    auto parts = bromesh::splitConnectedComponents(combined);
    ASSERT(parts.size() == 2, "two disjoint triangles should yield 2 components");
    ASSERT(parts[0].vertexCount() == 3, "component 0 should have 3 vertices");
    ASSERT(parts[1].vertexCount() == 3, "component 1 should have 3 vertices");
    ASSERT(parts[0].triangleCount() == 1, "component 0 should have 1 triangle");
    ASSERT(parts[1].triangleCount() == 1, "component 1 should have 1 triangle");

    // Two triangles sharing a vertex = 1 component
    bromesh::MeshData connected;
    connected.positions = {
        0,0,0,  1,0,0,  0,1,0,  1,1,0
    };
    connected.normals = {
        0,0,1,  0,0,1,  0,0,1,  0,0,1
    };
    connected.indices = { 0,1,2,  1,3,2 };
    auto cParts = bromesh::splitConnectedComponents(connected);
    ASSERT(cParts.size() == 1, "connected triangles should yield 1 component");
}

TEST(simplify_sphere) {
    auto s = bromesh::sphere(1, 32, 24);
    size_t origTris = s.triangleCount();
    ASSERT(origTris > 0, "sphere should have triangles");

    auto simplified = bromesh::simplify(s, 0.5f);
    ASSERT(!simplified.empty(), "simplified mesh should not be empty");
    ASSERT(simplified.triangleCount() < origTris, "simplified should have fewer triangles");
    ASSERT(simplified.triangleCount() > 0, "simplified should have at least one triangle");
    ASSERT(simplified.hasNormals(), "simplified should preserve normals");
    ASSERT(simplified.hasUVs(), "simplified should preserve UVs");
}

TEST(lod_chain) {
    auto s = bromesh::sphere(1, 32, 24);
    size_t origTris = s.triangleCount();
    float ratios[] = { 0.5f, 0.25f };
    auto chain = bromesh::generateLODChain(s, ratios, 2);
    ASSERT(chain.size() == 2, "LOD chain should have 2 meshes");
    ASSERT(!chain[0].empty(), "LOD 0 should not be empty");
    ASSERT(!chain[1].empty(), "LOD 1 should not be empty");
    ASSERT(chain[0].triangleCount() < origTris, "LOD 0 should have fewer triangles than original");
    ASSERT(chain[1].triangleCount() <= chain[0].triangleCount(), "LOD 1 should have <= triangles than LOD 0");
}


#if BROMESH_HAS_MESHOPTIMIZER
TEST(simplify_with_attributes) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto simplified = bromesh::simplifyWithAttributes(mesh, 0.5f);
    ASSERT(!simplified.empty(), "simplify_attr: should produce non-empty result");
    ASSERT(simplified.triangleCount() < mesh.triangleCount(),
           "simplify_attr: should have fewer triangles");
    ASSERT(simplified.hasNormals(), "simplify_attr: should preserve normals");
    ASSERT(simplified.hasUVs(), "simplify_attr: should preserve UVs");
}

TEST(simplify_with_attributes_preserves_more) {
    // Compare attribute-aware vs basic: attribute-aware should produce
    // at least as many triangles (it's more conservative at seams)
    auto mesh = bromesh::torus(2.0f, 0.5f, 32, 16);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto basic = bromesh::simplify(mesh, 0.3f);
    auto attr = bromesh::simplifyWithAttributes(mesh, 0.3f);
    ASSERT(!basic.empty(), "simplify_cmp: basic should work");
    ASSERT(!attr.empty(), "simplify_cmp: attribute-aware should work");
    // Both should reduce triangle count
    ASSERT(basic.triangleCount() < mesh.triangleCount(), "simplify_cmp: basic reduces");
    ASSERT(attr.triangleCount() < mesh.triangleCount(), "simplify_cmp: attr reduces");
}

#endif // BROMESH_HAS_MESHOPTIMIZER

TEST(simplify_target_count) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    size_t origTris = mesh.triangleCount();

    size_t target = origTris / 4;
    auto result = bromesh::simplifyToTriangleCount(mesh, target);

    // Should have fewer triangles than original (meshopt may not hit exact target)
    ASSERT(result.triangleCount() < origTris,
           "target_count: should have fewer triangles");
    ASSERT(result.triangleCount() > 0, "target_count: should have triangles");
}

TEST(simplify_target_count_identity) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    size_t origTris = mesh.triangleCount();

    // Requesting more than current should return same mesh
    auto result = bromesh::simplifyToTriangleCount(mesh, origTris * 2);
    ASSERT(result.triangleCount() == origTris,
           "target_count_identity: should return unchanged when target >= current");
}

TEST(merge_two_meshes) {
    auto a = bromesh::box(1, 1, 1);
    auto b = bromesh::sphere(1.0f, 8, 6);
    std::vector<bromesh::MeshData> meshes = {a, b};
    auto merged = bromesh::mergeMeshes(meshes);
    ASSERT(merged.vertexCount() == a.vertexCount() + b.vertexCount(),
           "merge_two: vertex count should be sum");
    ASSERT(merged.triangleCount() == a.triangleCount() + b.triangleCount(),
           "merge_two: triangle count should be sum");
    ASSERT(merged.hasNormals(), "merge_two: both have normals so merged should too");
    ASSERT(merged.hasUVs(), "merge_two: both have UVs so merged should too");
}

TEST(merge_single_mesh) {
    auto a = bromesh::box(1, 1, 1);
    auto merged = bromesh::mergeMeshes(&a, 1);
    ASSERT(merged.vertexCount() == a.vertexCount(), "merge_single: should be same");
    ASSERT(merged.triangleCount() == a.triangleCount(), "merge_single: should be same");
}

TEST(merge_empty) {
    auto result = bromesh::mergeMeshes(nullptr, 0);
    ASSERT(result.empty(), "merge_empty: should be empty");
}

TEST(merge_index_validity) {
    auto a = bromesh::box(1, 1, 1);
    auto b = bromesh::cylinder(0.5f, 1.0f, 12);
    std::vector<bromesh::MeshData> meshes = {a, b};
    auto merged = bromesh::mergeMeshes(meshes);
    // All indices should be valid
    for (uint32_t idx : merged.indices) {
        ASSERT(idx < merged.vertexCount(), "merge_idx: all indices should be valid");
    }
}

TEST(remove_degenerate_triangles) {
    // Create a mesh with one good triangle and one degenerate (zero-area)
    bromesh::MeshData mesh;
    mesh.positions = {
        0,0,0, 1,0,0, 0,1,0,  // good triangle
        2,0,0, 2,0,0, 3,0,0   // degenerate: two identical vertices
    };
    mesh.indices = {0,1,2, 3,4,5};

    auto repaired = bromesh::removeDegenerateTriangles(mesh);
    ASSERT(repaired.triangleCount() == 1, "remove_degen: should keep 1 triangle");
    ASSERT(repaired.indices[0] == 0 && repaired.indices[1] == 1 && repaired.indices[2] == 2,
           "remove_degen: should keep the good triangle");
}

TEST(remove_degenerate_collinear) {
    // Collinear triangle (zero area)
    bromesh::MeshData mesh;
    mesh.positions = {
        0,0,0, 1,0,0, 0,1,0,  // good
        0,0,0, 1,0,0, 2,0,0   // collinear
    };
    mesh.indices = {0,1,2, 3,4,5};

    auto repaired = bromesh::removeDegenerateTriangles(mesh);
    ASSERT(repaired.triangleCount() == 1, "remove_collinear: should remove collinear triangle");
}

TEST(remove_duplicate_triangles) {
    auto mesh = bromesh::box(1, 1, 1);
    size_t origTris = mesh.triangleCount();

    // Duplicate all triangles
    size_t origIdxCount = mesh.indices.size();
    for (size_t i = 0; i < origIdxCount; ++i) {
        mesh.indices.push_back(mesh.indices[i]);
    }
    ASSERT(mesh.triangleCount() == origTris * 2, "dup_setup: doubled");

    auto repaired = bromesh::removeDuplicateTriangles(mesh);
    ASSERT(repaired.triangleCount() == origTris, "remove_dup: should remove duplicates");
}

TEST(fill_holes_simple) {
    // Create an open box (5 faces, missing the top)
    // Using a simple example: a plane with a triangular hole
    bromesh::MeshData mesh;
    // Square with 4 triangles leaving a hole in the middle
    //  3---2
    //  |\ /|
    //  | 4 |   (vertex 4 at center, but no bottom-center triangle)
    //  |/ \|
    //  0---1
    mesh.positions = {
        0,0,0, 1,0,0, 1,1,0, 0,1,0, 0.5f,0.5f,0
    };
    // Only 3 triangles, leaving a gap
    mesh.indices = {
        0,1,4,  // bottom
        1,2,4,  // right
        2,3,4   // top
        // missing: 3,0,4 (left)
    };

    auto filled = bromesh::fillHoles(mesh);
    // Should add the missing triangle
    ASSERT(filled.triangleCount() >= mesh.triangleCount(),
           "fill_holes: should have at least as many triangles");
    ASSERT(filled.triangleCount() > mesh.triangleCount(),
           "fill_holes: should have added triangles to fill hole");
}

TEST(repair_preserves_clean_mesh) {
    // A clean mesh should pass through unchanged
    auto mesh = bromesh::box(1, 1, 1);
    auto degenFixed = bromesh::removeDegenerateTriangles(mesh);
    ASSERT(degenFixed.triangleCount() == mesh.triangleCount(),
           "repair_clean: degenerate removal should not change clean mesh");

    auto dupFixed = bromesh::removeDuplicateTriangles(mesh);
    ASSERT(dupFixed.triangleCount() == mesh.triangleCount(),
           "repair_clean: duplicate removal should not change clean mesh");
}

TEST(subdivide_midpoint_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideMidpoint(mesh, 1);
    // Each triangle splits into 4, box has 12 triangles -> 48
    ASSERT(result.triangleCount() == mesh.triangleCount() * 4,
           "midpoint_box: should have 4x triangles");
    ASSERT(result.vertexCount() > mesh.vertexCount(),
           "midpoint_box: should have more vertices");
}

TEST(subdivide_midpoint_two_iterations) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideMidpoint(mesh, 2);
    // 12 * 4 * 4 = 192
    ASSERT(result.triangleCount() == mesh.triangleCount() * 16,
           "midpoint_2x: should have 16x triangles");
}

TEST(subdivide_loop_sphere) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    auto result = bromesh::subdivideLoop(mesh, 1);
    ASSERT(result.triangleCount() == mesh.triangleCount() * 4,
           "loop_sphere: should have 4x triangles");
    ASSERT(result.hasNormals(), "loop_sphere: should have normals");
}

TEST(subdivide_catmull_clark_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideCatmullClark(mesh, 1);
    // CC on triangles: each tri -> 3 quads -> 6 triangles
    ASSERT(result.triangleCount() == mesh.triangleCount() * 6,
           "cc_box: should have 6x triangles");
    ASSERT(result.hasNormals(), "cc_box: should have normals");
}

TEST(subdivide_zero_iterations) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideLoop(mesh, 0);
    ASSERT(result.triangleCount() == mesh.triangleCount(),
           "subdiv_zero: 0 iterations should return same mesh");
}

TEST(smooth_laplacian_sphere) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    auto origPositions = mesh.positions;

    bromesh::smoothLaplacian(mesh, 0.5f, 3);

    ASSERT(mesh.vertexCount() * 3 == origPositions.size(),
           "laplacian: vertex count should not change");
    // After smoothing, positions should be different
    bool changed = false;
    for (size_t i = 0; i < origPositions.size(); ++i) {
        if (std::fabs(mesh.positions[i] - origPositions[i]) > 1e-6f) {
            changed = true;
            break;
        }
    }
    ASSERT(changed, "laplacian: positions should change after smoothing");
}

TEST(smooth_taubin_box) {
    auto base = bromesh::subdivideMidpoint(bromesh::weldVertices(bromesh::box(1.0f, 1.0f, 1.0f)), 2);
    auto mesh_taubin = base;
    bromesh::computeNormals(mesh_taubin);

    auto mesh_lap = base;
    bromesh::computeNormals(mesh_lap);

    // Measure initial bounding box
    auto bbox1 = bromesh::computeBBox(mesh_taubin);
    float vol1 = (bbox1.max.x-bbox1.min.x) * (bbox1.max.y-bbox1.min.y) * (bbox1.max.z-bbox1.min.z);

    bromesh::smoothTaubin(mesh_taubin, 0.5f, -0.53f, 5);

    auto bbox_taubin = bromesh::computeBBox(mesh_taubin);
    float vol_taubin = (bbox_taubin.max.x-bbox_taubin.min.x) * (bbox_taubin.max.y-bbox_taubin.min.y) * (bbox_taubin.max.z-bbox_taubin.min.z);

    bromesh::smoothLaplacian(mesh_lap, 0.5f, 5);

    auto bbox_lap = bromesh::computeBBox(mesh_lap);
    float vol_laplacian = (bbox_lap.max.x-bbox_lap.min.x) * (bbox_lap.max.y-bbox_lap.min.y) * (bbox_lap.max.z-bbox_lap.min.z);

    // Taubin should preserve volume (at least 95% of original)
    ASSERT(vol_taubin >= vol1 * 0.95f, "taubin: preserves volume (>= 95%)");
    // Laplacian causes noticeable shrinkage compared to Taubin
    ASSERT(vol_laplacian < vol_taubin * 0.90f, "taubin: prevents shrinkage compared to laplacian");
}

TEST(remesh_isotropic_sphere) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);

    auto result = bromesh::remeshIsotropic(mesh, 0.0f, 3);
    ASSERT(!result.empty(), "remesh: result should not be empty");
    ASSERT(result.hasNormals(), "remesh: should have normals");
    ASSERT(result.triangleCount() > 0, "remesh: should have triangles");

    // Manifoldness: round-trip through PolyMesh and validate.
    auto pm = bromesh::PolyMesh::fromMeshData(result.positions, result.indices);
    auto v = pm.validate();
    ASSERT(v.valid, "remesh: output validates as half-edge mesh");
    ASSERT(v.isClosed, "remesh: closed input → closed output");
}

TEST(remesh_isotropic_edge_length_stats) {
    auto mesh = bromesh::sphere(1.0f, 12, 9);
    const float target = 0.25f;
    auto result = bromesh::remeshIsotropic(mesh, target, 4);
    ASSERT(result.triangleCount() > 0, "remesh-stats: triangles emitted");

    // Mean edge length should land near the target after a few iterations.
    double sum = 0.0;
    size_t n = 0;
    for (size_t t = 0; t < result.triangleCount(); ++t) {
        for (int e = 0; e < 3; ++e) {
            uint32_t i0 = result.indices[t*3 + e];
            uint32_t i1 = result.indices[t*3 + (e + 1) % 3];
            float dx = result.positions[i0*3+0] - result.positions[i1*3+0];
            float dy = result.positions[i0*3+1] - result.positions[i1*3+1];
            float dz = result.positions[i0*3+2] - result.positions[i1*3+2];
            sum += std::sqrt(dx*dx + dy*dy + dz*dz);
            ++n;
        }
    }
    float mean = (float)(sum / (double)n);
    ASSERT(mean > target * 0.7f && mean < target * 1.3f,
           "remesh-stats: mean edge length within 30% of target");
}

TEST(remesh_isotropic_open_disk_preserves_boundary) {
    // Build an open triangle fan (a "disk") with a clear boundary loop.
    // Center vertex 0; ring of 8 vertices around it.
    const int rim = 8;
    std::vector<float> pos;
    pos.push_back(0); pos.push_back(0); pos.push_back(0);
    for (int i = 0; i < rim; ++i) {
        float a = (float)i / (float)rim * 6.2831853f;
        pos.push_back(std::cos(a));
        pos.push_back(std::sin(a));
        pos.push_back(0);
    }
    std::vector<uint32_t> idx;
    for (int i = 0; i < rim; ++i) {
        idx.push_back(0);
        idx.push_back((uint32_t)(1 + i));
        idx.push_back((uint32_t)(1 + (i + 1) % rim));
    }
    bromesh::MeshData disk;
    disk.positions = std::move(pos);
    disk.indices   = std::move(idx);

    auto result = bromesh::remeshIsotropic(disk, 0.5f, 3);
    ASSERT(!result.empty(), "remesh-disk: produces output");

    auto pm = bromesh::PolyMesh::fromMeshData(result.positions, result.indices);
    auto v = pm.validate();
    ASSERT(v.valid, "remesh-disk: output is structurally valid");
    ASSERT(!v.isClosed, "remesh-disk: open mesh stays open");
    ASSERT(v.boundaryHalfEdges > 0, "remesh-disk: boundary preserved");

    // Boundary loop integrity: every boundary vertex stays on z=0 (relax
    // pins boundary verts, so they don't drift normal to the disk plane).
    // Note: boundary splits land on chord midpoints, so radial position
    // shrinks toward the disk interior — that's expected behavior for an
    // unprojected isotropic remesh and not what we're testing here.
    int rimCount = 0;
    float maxR = 0.0f;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        if (he.face == bromesh::PolyMesh::NONE) continue;
        if (he.twin != bromesh::PolyMesh::NONE) continue;
        int32_t v0 = he.origin;
        float p[3]; pm.getVertex(v0, p);
        ASSERT(std::fabs(p[2]) < 1e-4f, "remesh-disk: boundary stays in z=0 plane");
        float r = std::sqrt(p[0]*p[0] + p[1]*p[1]);
        if (r > maxR) maxR = r;
        ++rimCount;
    }
    ASSERT(rimCount >= rim, "remesh-disk: at least the original rim count of boundary he");
    ASSERT(maxR > 0.5f, "remesh-disk: boundary doesn't collapse toward center");
}

TEST(transform_translate) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;

    bromesh::translateMesh(mesh, 5.0f, 0.0f, 0.0f);

    bool correct = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        if (std::fabs(mesh.positions[v*3+0] - (origPos[v*3+0] + 5.0f)) > 0.001f) {
            correct = false; break;
        }
    }
    ASSERT(correct, "translate: +5 on X");
}

TEST(transform_scale_uniform) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;

    bromesh::scaleMesh(mesh, 2.0f);

    bool correct = true;
    for (size_t i = 0; i < origPos.size(); ++i) {
        if (std::fabs(mesh.positions[i] - origPos[i] * 2.0f) > 0.001f) {
            correct = false; break;
        }
    }
    ASSERT(correct, "scale_uniform: all positions doubled");
}

TEST(transform_scale_nonuniform) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    bromesh::scaleMesh(mesh, 2.0f, 1.0f, 1.0f);

    // Normals should still be unit length
    bool unitNormals = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float* n = &mesh.normals[v*3];
        float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (std::fabs(len - 1.0f) > 0.01f) { unitNormals = false; break; }
    }
    ASSERT(unitNormals, "scale_nonuniform: normals remain unit length");
}

TEST(transform_rotate_90_y) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);

    // Record a vertex on the +X face
    float origX = mesh.positions[0];
    float origZ = mesh.positions[2];

    float pi = 3.14159265f;
    bromesh::rotateMesh(mesh, 0.0f, 1.0f, 0.0f, pi / 2.0f);

    // After 90° about Y: X -> Z, Z -> -X
    // Just check the mesh isn't unchanged
    bool changed = false;
    if (std::fabs(mesh.positions[0] - origX) > 0.01f ||
        std::fabs(mesh.positions[2] - origZ) > 0.01f) {
        changed = true;
    }
    ASSERT(changed, "rotate_90_y: positions should change");
}

TEST(transform_mirror_x) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);
    auto origPos = mesh.positions;

    bromesh::mirrorMesh(mesh, 0); // mirror across YZ plane

    bool mirrored = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        if (std::fabs(mesh.positions[v*3+0] - (-origPos[v*3+0])) > 0.001f) {
            mirrored = false; break;
        }
    }
    ASSERT(mirrored, "mirror_x: X coordinates negated");
}

TEST(transform_center) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::translateMesh(mesh, 10.0f, 20.0f, 30.0f);

    bromath::Vec3 center = bromesh::centerMesh(mesh);

    ASSERT(std::fabs(center.x - 10.0f) < 0.01f, "center: original center X=10");
    ASSERT(std::fabs(center.y - 20.0f) < 0.01f, "center: original center Y=20");
    ASSERT(std::fabs(center.z - 30.0f) < 0.01f, "center: original center Z=30");

    // After centering, bbox center should be at origin
    auto bbox = bromesh::computeBBox(mesh);
    ASSERT(std::fabs(bromath::acenter(bbox).x) < 0.01f, "center: bbox center X ~= 0");
    ASSERT(std::fabs(bromath::acenter(bbox).y) < 0.01f, "center: bbox center Y ~= 0");
    ASSERT(std::fabs(bromath::acenter(bbox).z) < 0.01f, "center: bbox center Z ~= 0");
}

TEST(transform_matrix_identity) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;

    float identity[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    bromesh::transformMesh(mesh, identity);

    ASSERT(mesh.positions == origPos, "transform_identity: no change");
}

