#include "test_framework.h"
#include <cmath>

TEST(polymesh_compact_removes_dead_face) {
    // Build a closed cube-equivalent (use sphere with low resolution → fully
    // closed manifold) and tombstone one face via deleteFace. After compact
    // the indices should be dense and validate() should still succeed,
    // though the mesh now has a triangular hole (boundary).
    auto src = bromesh::sphere(1.0f, 8, 6);
    auto pm = bromesh::PolyMesh::fromMeshData(src.positions, src.indices);

    size_t faces0 = pm.faceCount();
    size_t he0    = pm.halfEdgeCount();

    auto v0 = pm.validate();
    ASSERT(v0.valid && v0.isClosed, "polymesh_compact: source mesh validates closed");

    pm.deleteFace(0);
    ASSERT(!pm.isLiveFace(0), "polymesh_compact: face 0 tombstoned");

    // After deleteFace, validate() must skip tombstones and stay valid.
    auto v1 = pm.validate();
    ASSERT(v1.valid, "polymesh_compact: validate skips tombstones");
    ASSERT(!v1.isClosed, "polymesh_compact: deleting a face opens the mesh");
    ASSERT(v1.boundaryHalfEdges == 3, "polymesh_compact: 3 boundary he after one face removed");

    pm.compact();
    ASSERT(pm.faceCount()    == faces0 - 1, "polymesh_compact: face count -1");
    ASSERT(pm.halfEdgeCount() == he0 - 3,    "polymesh_compact: 3 he removed");

    auto v2 = pm.validate();
    ASSERT(v2.valid, "polymesh_compact: post-compact validates");
    ASSERT(v2.boundaryHalfEdges == 3, "polymesh_compact: 3 boundary he preserved across compact");

    // No tombstones remain: every he has face != NONE, every face has he != NONE.
    bool anyDeadHe = false;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i)
        if (!pm.isLiveHalfEdge((int32_t)i)) { anyDeadHe = true; break; }
    ASSERT(!anyDeadHe, "polymesh_compact: no dead halfEdges after compact");

    bool anyDeadFace = false;
    for (size_t i = 0; i < pm.faceCount(); ++i)
        if (!pm.isLiveFace((int32_t)i)) { anyDeadFace = true; break; }
    ASSERT(!anyDeadFace, "polymesh_compact: no dead faces after compact");
}

TEST(polymesh_compact_idempotent) {
    auto src = bromesh::sphere(1.0f, 8, 6);
    auto pm = bromesh::PolyMesh::fromMeshData(src.positions, src.indices);
    size_t v0 = pm.vertexCount(), f0 = pm.faceCount(), h0 = pm.halfEdgeCount();
    pm.compact();
    ASSERT(pm.vertexCount()   == v0, "polymesh_compact_idempotent: verts unchanged on no-op");
    ASSERT(pm.faceCount()     == f0, "polymesh_compact_idempotent: faces unchanged on no-op");
    ASSERT(pm.halfEdgeCount() == h0, "polymesh_compact_idempotent: he unchanged on no-op");
    auto v = pm.validate();
    ASSERT(v.valid && v.isClosed, "polymesh_compact_idempotent: still validates");
}

TEST(polymesh_split_edge_interior) {
    // Two triangles sharing an edge: (a,b,c) + (b,a,d) with a-b shared.
    //   a = (0,0,0), b = (1,0,0), c = (0.5, 1, 0), d = (0.5, -1, 0)
    std::vector<float> pos = {
        0,0,0,    1,0,0,    0.5f,1,0,    0.5f,-1,0
    };
    std::vector<uint32_t> idx = { 0,1,2,  1,0,3 };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    auto v0 = pm.validate();
    ASSERT(v0.valid, "split_interior: setup validates");

    // Find a half-edge on the shared a-b edge (the one with a twin).
    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        if (pm.halfEdges()[i].twin != bromesh::PolyMesh::NONE) { hi = (int32_t)i; break; }
    }
    ASSERT(hi >= 0, "split_interior: found shared edge");

    size_t v_before = pm.vertexCount();
    size_t f_before = pm.faceCount();
    size_t h_before = pm.halfEdgeCount();

    int32_t M = pm.splitEdge(hi);
    ASSERT(M >= 0, "split_interior: returns new vertex");
    ASSERT(pm.vertexCount()   == v_before + 1, "split_interior: +1 vertex");
    ASSERT(pm.faceCount()     == f_before + 2, "split_interior: +2 faces");
    ASSERT(pm.halfEdgeCount() == h_before + 6, "split_interior: +6 halfEdges");

    auto v1 = pm.validate();
    ASSERT(v1.valid, "split_interior: post-split validates");
    // Mesh has 4 boundary half-edges originally (the outer perimeter); after
    // splitting the *interior* shared edge, the same 4 boundary he's plus 0
    // new boundaries. The new midpoint vertex is interior.
    ASSERT(v1.boundaryHalfEdges == 4, "split_interior: boundary count unchanged");

    // Midpoint position check.
    float p[3]; pm.getVertex(M, p);
    ASSERT(std::fabs(p[0] - 0.5f) < 1e-6f &&
           std::fabs(p[1] - 0.0f) < 1e-6f &&
           std::fabs(p[2] - 0.0f) < 1e-6f, "split_interior: midpoint position");
}

TEST(polymesh_split_edge_boundary) {
    // Single triangle (a,b,c). Split a-b (a boundary edge).
    std::vector<float> pos = { 0,0,0,  1,0,0,  0.5f,1,0 };
    std::vector<uint32_t> idx = { 0,1,2 };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    auto v0 = pm.validate();
    ASSERT(v0.valid, "split_boundary: setup validates");
    ASSERT(v0.boundaryHalfEdges == 3, "split_boundary: 3 boundary he initially");

    // Pick the half-edge from vertex 0 (origin a). All 3 he are boundary.
    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        if (pm.halfEdges()[i].origin == 0) { hi = (int32_t)i; break; }
    }
    ASSERT(hi >= 0, "split_boundary: found edge from vertex 0");

    size_t v_before = pm.vertexCount();
    size_t f_before = pm.faceCount();
    size_t h_before = pm.halfEdgeCount();

    int32_t M = pm.splitEdge(hi);
    ASSERT(M >= 0, "split_boundary: returns new vertex");
    ASSERT(pm.vertexCount()   == v_before + 1, "split_boundary: +1 vertex");
    ASSERT(pm.faceCount()     == f_before + 1, "split_boundary: +1 face");
    ASSERT(pm.halfEdgeCount() == h_before + 3, "split_boundary: +3 halfEdges");

    auto v1 = pm.validate();
    ASSERT(v1.valid, "split_boundary: validates post-split");
    // Was 3 boundary half-edges (perimeter); after splitting one boundary
    // segment in two, we have 4.
    ASSERT(v1.boundaryHalfEdges == 4, "split_boundary: 4 boundary he after");
}

TEST(polymesh_flip_edge) {
    // Two triangles sharing edge a-b. Flip; verify the new diagonal is c-d
    // and topology remains valid. Flip again returns to original-ish.
    std::vector<float> pos = {
        0,0,0,    1,0,0,    0.5f,1,0,    0.5f,-1,0
    };
    std::vector<uint32_t> idx = { 0,1,2,  1,0,3 };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        if (pm.halfEdges()[i].twin != bromesh::PolyMesh::NONE) { hi = (int32_t)i; break; }
    }
    ASSERT(hi >= 0, "flip: shared edge found");

    size_t v_before = pm.vertexCount();
    size_t f_before = pm.faceCount();
    size_t h_before = pm.halfEdgeCount();

    bool ok = pm.flipEdge(hi);
    ASSERT(ok, "flip: succeeded");
    ASSERT(pm.vertexCount()   == v_before, "flip: vertex count unchanged");
    ASSERT(pm.faceCount()     == f_before, "flip: face count unchanged");
    ASSERT(pm.halfEdgeCount() == h_before, "flip: he count unchanged");

    auto v1 = pm.validate();
    ASSERT(v1.valid, "flip: validates");
    ASSERT(v1.boundaryHalfEdges == 4, "flip: boundary count unchanged");

    // The new diagonal must connect 2 (=c) and 3 (=d).
    bool foundCD = false;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        int32_t dest = pm.halfEdges()[he.next].origin;
        if ((he.origin == 2 && dest == 3) || (he.origin == 3 && dest == 2)) { foundCD = true; break; }
    }
    ASSERT(foundCD, "flip: new diagonal c-d exists");

    // Original a-b edge should no longer be present.
    bool foundAB = false;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        int32_t dest = pm.halfEdges()[he.next].origin;
        if ((he.origin == 0 && dest == 1) || (he.origin == 1 && dest == 0)) { foundAB = true; break; }
    }
    ASSERT(!foundAB, "flip: old diagonal a-b removed");
}

TEST(polymesh_flip_rejects_boundary) {
    // Single triangle — every edge is boundary.
    std::vector<float> pos = { 0,0,0,  1,0,0,  0.5f,1,0 };
    std::vector<uint32_t> idx = { 0,1,2 };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    bool ok = pm.flipEdge(0);
    ASSERT(!ok, "flip: rejects boundary edge");
    ASSERT(pm.faceCount() == 1, "flip: no mutation on rejection");
}

TEST(polymesh_flip_rejects_existing_diagonal) {
    // Tetrahedron: every edge already has both diagonals' worth of
    // neighbors. Specifically, in a tetrahedron, for edge a-b, the
    // opposite vertices c and d ARE already connected. So any flip
    // request must be refused.
    std::vector<float> pos = {
        0,0,0,  1,0,0,  0,1,0,  0,0,1
    };
    std::vector<uint32_t> idx = {
        0,1,2,  0,2,3,  0,3,1,  1,3,2
    };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    auto v0 = pm.validate();
    ASSERT(v0.valid && v0.isClosed, "flip-tet: setup valid + closed");

    int hadFlipReject = 0, attempts = 0;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        if (pm.halfEdges()[i].twin == bromesh::PolyMesh::NONE) continue;
        ++attempts;
        if (!pm.flipEdge((int32_t)i)) ++hadFlipReject;
    }
    ASSERT(attempts > 0, "flip-tet: tried at least one flip");
    ASSERT(hadFlipReject == attempts, "flip-tet: every flip rejected on tetrahedron");
}

TEST(polymesh_collapse_interior) {
    // Octahedron: 6 verts, 8 faces, 12 edges. Collapsing one interior edge
    // (between two equator vertices) gives 5 verts, 6 faces.
    std::vector<float> pos = {
         1, 0, 0,
        -1, 0, 0,
         0, 1, 0,
         0,-1, 0,
         0, 0, 1,
         0, 0,-1
    };
    std::vector<uint32_t> idx = {
        0,2,4,  2,1,4,  1,3,4,  3,0,4,
        2,0,5,  1,2,5,  3,1,5,  0,3,5
    };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);
    auto v0 = pm.validate();
    ASSERT(v0.valid && v0.isClosed, "collapse: octahedron valid + closed");

    // Find an edge between equator vertices (0..3). Avoid pole-incident
    // edges to keep both endpoints interior with the simple link condition.
    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        int32_t dest = pm.halfEdges()[he.next].origin;
        if (he.origin <= 3 && dest <= 3 && he.origin != dest &&
            he.twin != bromesh::PolyMesh::NONE) {
            hi = (int32_t)i;
            break;
        }
    }
    ASSERT(hi >= 0, "collapse: equator edge found");

    bool ok = pm.collapseEdge(hi);
    ASSERT(ok, "collapse: succeeded");

    auto v1 = pm.validate();
    ASSERT(v1.valid, "collapse: validates with tombstones");

    pm.compact();
    auto v2 = pm.validate();
    ASSERT(v2.valid && v2.isClosed, "collapse: post-compact valid + closed");
    ASSERT(pm.vertexCount() == 5, "collapse: 6 - 1 = 5 vertices");
    ASSERT(pm.faceCount()   == 6, "collapse: 8 - 2 = 6 faces");
}

TEST(polymesh_collapse_link_condition) {
    // Build the bowtie config: vertices a,b share two non-face neighbors x,y.
    // Concretely: triangles (a,b,x), (b,a,y), and *also* (a,b,z), (b,a,w).
    // No — that would just be 4 faces around an edge (non-manifold).
    //
    // Standard link-condition violator: a and b share neighbor x via two
    // separate faces, where x is NOT a tip of the two faces sharing a-b.
    //
    // Setup: two disconnected triangle pairs sharing only the edge a-b in
    // both, plus an extra edge a-x and b-x where x is a fifth vertex not
    // adjacent to either face's tip. We embed this in a closed mesh by
    // building it as: square pyramid base where the apex is shared.
    //
    // Easier: an interior edge of an octahedron whose two endpoints share
    // both their tip-vertices already covers the legal case. To trigger
    // a link-failure, we use a flat "fan" of 4 triangles:
    //   verts: 0=a (center), 1=b, 2=x, 3=y, 4=z
    //   tris: (0,1,2),(0,2,3),(0,3,4),(0,4,1)
    // Now collapse 0-1: a's neighbors {1,2,3,4}; b=1's neighbors {0,2,4}.
    // Common = {2,4}. Tips of faces sharing a-b are (0,1,2) tip=2 and
    // (0,4,1) tip=4. So both common neighbors are tips. Link OK.
    //
    // To force failure: add an extra edge between b and 3 by inserting
    // another tri (1,2,3). That makes 1 and 0 both adjacent to 3, but 3
    // is not a tip of either face containing edge 0-1.
    std::vector<float> pos = {
         0, 0, 0,    // 0 = a (center)
         1, 0, 0,    // 1 = b
         0, 1, 0,    // 2 = x
        -1, 0, 0,    // 3 = y
         0,-1, 0     // 4 = z
    };
    std::vector<uint32_t> idx = {
        0,1,2,  0,2,3,  0,3,4,  0,4,1,
        1,2,3   // extra triangle making 1 ~ 3 (not a tip of an a-b face)
    };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        int32_t dest = pm.halfEdges()[he.next].origin;
        if (he.origin == 0 && dest == 1 && he.twin != bromesh::PolyMesh::NONE) {
            hi = (int32_t)i;
            break;
        }
    }
    ASSERT(hi >= 0, "collapse_link: a-b he found");

    bool ok = pm.collapseEdge(hi);
    ASSERT(!ok, "collapse_link: refused (vertex 3 is non-tip common neighbor)");
}

TEST(polymesh_collapse_boundary) {
    // Two triangles sharing edge a-b; collapse one of the boundary edges.
    //   tris (0,1,2) and (1,0,3). Edge 0-2 is boundary (only in face 0).
    std::vector<float> pos = {
        0,0,0,  1,0,0,  0.5f,1,0,  0.5f,-1,0
    };
    std::vector<uint32_t> idx = { 0,1,2,  1,0,3 };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        if (he.twin != bromesh::PolyMesh::NONE) continue;
        if (he.origin == 0) { hi = (int32_t)i; break; }
    }
    ASSERT(hi >= 0, "collapse_boundary: boundary edge from vertex 0 found");

    bool ok = pm.collapseEdge(hi);
    ASSERT(ok, "collapse_boundary: succeeded");

    pm.compact();
    auto v = pm.validate();
    ASSERT(v.valid, "collapse_boundary: post-compact validates");
    ASSERT(pm.vertexCount() == 3, "collapse_boundary: -1 vertex");
    ASSERT(pm.faceCount()   == 1, "collapse_boundary: -1 face");
}

TEST(polymesh_collapse_refuses_interior_with_boundary_vertex) {
    std::vector<float> pos = {
         0, 0, 0,    // 0 = center
         1, 0, 0,    // 1
         0, 1, 0,    // 2
        -1, 0, 0,    // 3
         0,-1, 0     // 4
    };
    std::vector<uint32_t> idx = {
        0,1,2,  0,2,3,  0,3,4,  0,4,1
    };
    auto pm = bromesh::PolyMesh::fromMeshData(pos, idx);

    int32_t hi = -1;
    for (size_t i = 0; i < pm.halfEdgeCount(); ++i) {
        const auto& he = pm.halfEdges()[i];
        int32_t dest = pm.halfEdges()[he.next].origin;
        if (he.origin == 0 && dest == 1 && he.twin != bromesh::PolyMesh::NONE) {
            hi = (int32_t)i;
            break;
        }
    }
    ASSERT(hi >= 0, "collapse_refuse_boundary: interior 0-1 he found");

    ASSERT(pm.isBoundaryVertex(1), "vertex 1 is on boundary");
    ASSERT(!pm.isBoundaryVertex(0), "vertex 0 is interior");

    bool ok = pm.collapseEdge(hi);
    ASSERT(!ok, "collapse_refuse_boundary: refused");
}

TEST(polymesh_inset_face) {
    // Cube with 8 vertices, 6 quad faces
    std::vector<float> pos = {
        -1.0f, -1.0f,  1.0f,   1.0f, -1.0f,  1.0f,   1.0f,  1.0f,  1.0f,  -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f, -1.0f,   1.0f, -1.0f, -1.0f,   1.0f,  1.0f, -1.0f,  -1.0f,  1.0f, -1.0f
    };
    std::vector<uint32_t> polyVerts = {
        0, 1, 2, 3, // front (+Z)
        5, 4, 7, 6, // back (-Z)
        4, 0, 3, 7, // left (-X)
        1, 5, 6, 2, // right (+X)
        3, 2, 6, 7, // top (+Y)
        4, 5, 1, 0  // bottom (-Y)
    };
    std::vector<uint32_t> polyOffsets = { 0, 4, 8, 12, 16, 20, 24 };
    auto pm = bromesh::PolyMesh::fromPolygons(pos, polyVerts, polyOffsets);

    auto v0 = pm.validate();
    ASSERT(v0.valid && v0.isClosed, "inset_face: cube validates closed");
    ASSERT(pm.faceCount() == 6, "inset_face: 6 initial quad faces");

    auto res = pm.insetFace(0, 0.2f, true);
    ASSERT(res.innerFace >= 0, "inset_face: inner face created");
    ASSERT(res.innerVerts.size() == 4, "inset_face: 4 inner vertices");
    ASSERT(res.bridgeFaces.size() == 4, "inset_face: 4 bridge quad faces created");

    auto v1 = pm.validate();
    ASSERT(v1.valid && v1.isClosed, "inset_face: post-inset validates closed");
    ASSERT(pm.faceCount() == 10, "inset_face: 6 - 1 + 1 + 4 = 10 faces");

    for (int32_t bf : res.bridgeFaces) {
        ASSERT(pm.faceVertexCount(bf) == 4, "inset_face: bridge face is quad");
    }
    ASSERT(pm.faceVertexCount(res.innerFace) == 4, "inset_face: inner face is quad");

    auto tess = pm.tessellate();
    ASSERT(!tess.indices.empty() && tess.indices.size() % 3 == 0, "inset_face: tessellate validates");
}

TEST(polymesh_inset_and_extrude) {
    std::vector<float> pos = {
        -1.0f, -1.0f,  1.0f,   1.0f, -1.0f,  1.0f,   1.0f,  1.0f,  1.0f,  -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f, -1.0f,   1.0f, -1.0f, -1.0f,   1.0f,  1.0f, -1.0f,  -1.0f,  1.0f, -1.0f
    };
    std::vector<uint32_t> polyVerts = {
        0, 1, 2, 3, // front (+Z)
        5, 4, 7, 6, // back (-Z)
        4, 0, 3, 7, // left (-X)
        1, 5, 6, 2, // right (+X)
        3, 2, 6, 7, // top (+Y)
        4, 5, 1, 0  // bottom (-Y)
    };
    std::vector<uint32_t> polyOffsets = { 0, 4, 8, 12, 16, 20, 24 };
    auto pm = bromesh::PolyMesh::fromPolygons(pos, polyVerts, polyOffsets);

    // Inset face 0 (front face)
    auto insetRes = pm.insetFace(0, 0.25f, true);
    ASSERT(insetRes.innerFace >= 0, "inset_extrude: inset created inner face");

    // Extrude inner face along +Z
    float offset[3] = { 0.0f, 0.0f, 0.5f };
    auto extRes = pm.extrudeFace(insetRes.innerFace, offset, false);
    ASSERT(extRes.bridgeFaces.size() == 4, "inset_extrude: 4 bridge faces extruded");

    auto v = pm.validate();
    ASSERT(v.valid && v.isClosed, "inset_extrude: topology validates closed");

    auto tess = pm.tessellate();
    ASSERT(!tess.indices.empty(), "inset_extrude: tessellate non-empty");

    bromesh::MeshData md;
    md.positions = tess.positions;
    md.indices = tess.indices;
    ASSERT(bromesh::isManifold(md), "inset_extrude: manifold verified");
    ASSERT(bromesh::computeVolume(md) > 8.0f, "inset_extrude: extruded volume exceeds base cube (8.0)");
}


