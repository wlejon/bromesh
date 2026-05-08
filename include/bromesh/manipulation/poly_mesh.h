#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <string>
#include <vector>

namespace bromesh {

// ---------------------------------------------------------------------------
// PolyMesh — half-edge adjacency over N-gon faces.
//
// The "interaction" / "edit" mesh sitting alongside the rendered triangle
// mesh. Faces are first-class N-gons (3..N vertices). Tessellation produces
// the render-ready triangle buffer on demand; the edit topology never sees
// arbitrary triangulation choices, so push/pull, bevel, inset, etc. operate
// on the polygon graph and downstream triangulation stays clean.
//
// All references between Vertex/HalfEdge/Face are integer indices into the
// owning vectors — never raw pointers — so vector growth (extrude, bridge)
// doesn't invalidate them, and JS bindings can use the same indices as
// stable handles.
//
// Topology invariants (post-construction or post-rematch):
//   - For every he H: he[H].next.next.next == H if and only if H's face is
//     a triangle; in general, walking .next around face[F].halfEdge returns
//     to itself in faceVertexCount(F) steps.
//   - For every he H with twin T != -1: he[T].twin == H.
//   - Boundary half-edges have twin == -1.
//
// Tombstones (transient state during surgery):
//   - A half-edge with face == NONE is dead. Its origin/next are also set
//     to NONE for clarity. All references (twin, vertex.halfEdge,
//     face.halfEdge) into dead half-edges must be severed by the surgery
//     operation that creates the tombstone.
//   - A face with halfEdge == NONE is dead.
//   - Tombstoned entries are skipped by validate() and other walkers; they
//     are removed (and indices renumbered) by compact().
// ---------------------------------------------------------------------------

class PolyMesh {
public:
    static constexpr int32_t NONE = -1;

    struct Vertex {
        float x = 0, y = 0, z = 0;
        int32_t halfEdge = NONE;   // any outgoing he; NONE if isolated
    };

    struct HalfEdge {
        int32_t origin = NONE;     // vertex index this he points away from
        int32_t twin   = NONE;     // opposing he across the edge; NONE = boundary
        int32_t next   = NONE;     // next he around the same face (CCW)
        int32_t face   = NONE;     // owning face index
    };

    struct Face {
        int32_t halfEdge = NONE;   // any bounding he
        int32_t group    = -1;     // user-assigned group tag (e.g. coplanar grouping)
    };

    // ── Construction ──────────────────────────────────────────────────────

    /// Build from a triangle MeshData. Each tri becomes a triangular face.
    /// `triToGroup` (optional, must be size == triangleCount if non-empty)
    /// tags each tri's face with a group id; later operations can merge or
    /// query faces by group.
    ///
    /// Twin matching runs in two passes: first by vertex *index* equality,
    /// then by *quantized position* equality. The latter closes hard-edge
    /// seams where two faces use distinct vertex indices at the same world
    /// point (e.g. a cube with per-face split corners, or cylinder caps that
    /// duplicate the rim vertex for UV continuity).
    static PolyMesh fromMeshData(const std::vector<float>& positions,
                                  const std::vector<uint32_t>& indices,
                                  const std::vector<int32_t>& triToGroup = {});

    /// Build a single N-gon face from a planar polygon. Vertices are inserted
    /// in the order given; the polygon must be simple (non-self-intersecting)
    /// and CCW as viewed from +normal. The resulting mesh has one face, N
    /// vertices, N half-edges (all boundary).
    static PolyMesh fromPolygon(const std::vector<float>& outerXYZ,
                                 const float normal[3],
                                 int32_t group = 0);

    /// Build from N-gon face soup. `polyOffsets` is length F+1 and gives the
    /// start index of each face's vertex list inside `polyVerts`. Vertices
    /// are global; positions[i] supplies the coords for vertex i. Twins are
    /// matched by vertex identity then by position. Optional per-face group
    /// tags.
    static PolyMesh fromPolygons(const std::vector<float>& positions,
                                  const std::vector<uint32_t>& polyVerts,
                                  const std::vector<uint32_t>& polyOffsets,
                                  const std::vector<int32_t>& faceGroups = {});

    // ── Inspection ────────────────────────────────────────────────────────

    size_t vertexCount()   const { return vertices_.size(); }
    size_t halfEdgeCount() const { return halfEdges_.size(); }
    size_t faceCount()     const { return faces_.size(); }

    /// Liveness checks. After surgery (e.g. collapseEdge, deleteFace) some
    /// entries are tombstoned and walkers should skip them; compact() drops
    /// them and renumbers the survivors.
    bool isLiveHalfEdge(int32_t hi) const {
        return hi >= 0 && hi < (int32_t)halfEdges_.size() &&
               halfEdges_[hi].face != NONE;
    }
    bool isLiveFace(int32_t fi) const {
        return fi >= 0 && fi < (int32_t)faces_.size() &&
               faces_[fi].halfEdge != NONE;
    }

    /// A boundary half-edge has no twin (twin == NONE). True for both the
    /// half-edge representing the outer rim of an open mesh and for any
    /// internal half-edge whose adjacent face has been removed.
    bool isBoundaryHalfEdge(int32_t hi) const {
        return isLiveHalfEdge(hi) && halfEdges_[hi].twin == NONE;
    }

    /// True if any live half-edge incident to `vi` (incoming or outgoing)
    /// is a boundary half-edge. O(halfEdgeCount); used to gate boundary
    /// behavior in surgery operations.
    bool isBoundaryVertex(int32_t vi) const;

    const std::vector<Vertex>&   vertices()  const { return vertices_; }
    const std::vector<HalfEdge>& halfEdges() const { return halfEdges_; }
    const std::vector<Face>&     faces()     const { return faces_; }

    /// Number of vertices bounding a face (3 for a triangle, 4 for a quad…).
    int faceVertexCount(int faceIdx) const;

    /// Ordered list of vertex indices around a face (CCW from face normal).
    std::vector<int32_t> faceVertices(int faceIdx) const;

    /// Ordered list of half-edge indices around a face.
    std::vector<int32_t> faceHalfEdges(int faceIdx) const;

    void getVertex(int vi, float out[3]) const;

    /// Best-fit / cross-product face normal computed from current vertex
    /// positions. Unit-length unless the face is fully degenerate (then zero).
    void computeFaceNormal(int faceIdx, float out[3]) const;

    /// Faces in a given group. Linear scan; O(F).
    std::vector<int32_t> facesInGroup(int groupId) const;

    /// Boundary-of-face-or-group as ordered half-edge loops.
    /// `gIdx == NONE` (default behavior of the int variant) returns the
    /// boundary of a single face by index (faceIdx). Otherwise returns the
    /// outer + hole loops of all faces sharing `gIdx`.
    std::vector<std::vector<int32_t>> findFaceBoundary(int faceIdx) const;
    std::vector<std::vector<int32_t>> findGroupBoundary(int groupId) const;

    // ── Tessellation ──────────────────────────────────────────────────────

    struct Tessellation {
        std::vector<float>    positions;   // xyz, stride 3
        std::vector<float>    normals;     // xyz, stride 3 (per-face flat)
        std::vector<uint32_t> indices;     // stride 3
        /// triToFace[t] = source face index that produced triangle t.
        std::vector<int32_t>  triToFace;
        /// triToGroup[t] = group of triToFace[t] (convenience).
        std::vector<int32_t>  triToGroup;
    };

    /// Tessellate every face into triangles. Each face is projected to its
    /// best-fit plane and triangulated with a robust ear-clipper that
    /// handles concave polygons. Face vertex positions are emitted directly
    /// (no welding across faces — caller can run weld() if desired).
    /// Per-face flat normals are emitted to `normals`.
    Tessellation tessellate() const;

    // ── Validation ────────────────────────────────────────────────────────

    struct Validation {
        bool valid              = true;   // structural ok (faces close, no nulls)
        bool isClosed           = false;  // every he has a twin
        int  boundaryHalfEdges  = 0;
        std::vector<std::string> errors;
    };
    Validation validate() const;

    // ── Surgery ───────────────────────────────────────────────────────────

    /// Append a fresh vertex; returns its index.
    int32_t addVertex(float x, float y, float z);

    /// Append a face from an ordered vertex-index loop (>= 3 verts). Twin
    /// pointers are NOT auto-rematched — call rematchTwins() afterwards.
    /// Returns the new face index.
    int32_t addFace(const std::vector<int32_t>& verts, int32_t group = -1);

    /// Translate a single vertex.
    void translateVertex(int32_t vi, const float offset[3]);

    /// Move every vertex of a face by `offset`. Topology preserved; for an
    /// open boundary face this trivially extrudes one corner. To keep
    /// adjacent geometry attached use translateFaceWithRing for the
    /// vertex-substitution path.
    void translateFace(int32_t faceIdx, const float offset[3]);

    /// Vertex-substitution push/pull for a CLOSED solid: every vertex on
    /// face `faceIdx` (and any seam-duplicate vertices at the same position
    /// in adjacent faces, found via twin.next ring walk) translates by
    /// `offset`. Adjacent faces deform in place via shared vertex
    /// references. Manifold preserved.
    void translateFaceWithRing(int32_t faceIdx, const float offset[3]);

    struct ExtrudeResult {
        std::vector<int32_t> dupVerts;       // boundary-vertex duplicates (parallel to original loop verts)
        std::vector<int32_t> bridgeFaces;    // one new quad face per boundary edge
        int32_t              backFace = -1;  // back-face copy of the moved face (closes the slab)
        /// For each new bridge face, the original adjacent group across the
        /// boundary edge (-1 for mesh boundary). Useful for caller-side
        /// group merging.
        std::vector<int32_t> bridgeAdjGroup;
    };

    /// SketchUp-style extrusion of a single face: duplicate the face's
    /// boundary verts at `offset`, rewire the face to use the duplicates,
    /// emit one quad bridge face per boundary edge, and (if `withBackFace`)
    /// emit a reversed back face at the original positions to close the
    /// slab. New bridge faces get group = `bridgeGroup` (default: assign
    /// a fresh per-bridge group id starting at `nextGroup`); back face gets
    /// `backGroup`. Twin pointers are rematched at the end.
    ExtrudeResult extrudeFace(int32_t faceIdx,
                               const float offset[3],
                               bool withBackFace      = true,
                               int32_t bridgeGroup    = -1,
                               int32_t backGroup      = -1);

    /// Split the edge identified by half-edge `hi`. Inserts a new vertex
    /// at `posOptional` (if non-null) or at the midpoint of the edge's
    /// endpoints. Both adjacent faces (or the single face on a boundary
    /// edge) must be triangles — otherwise returns NONE without modifying
    /// the mesh.
    ///
    /// Net change for an interior edge: +1 vertex, +2 faces, +6 half-edges.
    /// Net change for a boundary edge:  +1 vertex, +1 face,  +3 half-edges.
    /// Twin links across the original edge are preserved (rewired to the
    /// new midpoint vertex on each side); all unrelated twins are unchanged.
    ///
    /// Returns the new vertex index, or NONE on precondition failure.
    int32_t splitEdge(int32_t hi, const float* posOptional = nullptr);

    /// Flip the diagonal of two adjacent triangles. The edge identified by
    /// `hi` separates triangles (a,b,c) and (b,a,d); after flipping, the
    /// diagonal becomes c-d and the triangles become (a,d,c) and (b,c,d).
    ///
    /// Refused (returns false, no mutation) when:
    ///   - hi is dead, on a boundary, or the adjacent faces aren't triangles
    ///   - the new diagonal vertices c,d are already connected by another
    ///     edge in the mesh (flipping would produce a non-manifold edge)
    ///
    /// No vertices, faces, or half-edges are added or deleted — only the
    /// origin/next/face pointers of the six half-edges of the two faces
    /// are rewired.
    bool flipEdge(int32_t hi);

    /// Collapse the edge identified by half-edge `hi`, merging its
    /// destination vertex into its origin. The surviving vertex moves to
    /// `posOptional` (if non-null) or to the midpoint of the original
    /// edge.
    ///
    /// Returns false (no mutation) when:
    ///   - hi is dead, or its adjacent face(s) aren't triangles
    ///   - link condition fails (a and b share a common neighbor that
    ///     isn't an opposite-face tip — collapsing would create a
    ///     duplicated edge or non-manifold)
    ///   - hi is interior (has a twin) but either endpoint is a boundary
    ///     vertex (collapsing would degenerate a boundary loop)
    ///
    /// Net change for an interior edge: -1 vertex, -2 faces, -6 half-edges.
    /// Net change for a boundary edge:  -1 vertex, -1 face,  -3 half-edges.
    /// Tombstones the removed faces/half-edges; call compact() to drop them.
    bool collapseEdge(int32_t hi, const float* posOptional = nullptr);

    /// Low-level: tombstone a face and all its half-edges. Any outside twin
    /// of those half-edges is severed (its `.twin` becomes NONE), so the
    /// surviving side becomes a boundary half-edge. Vertex outgoing-HE
    /// pointers that referenced any of the now-dead half-edges are reseeded
    /// to a surviving outgoing HE if one exists, else set to NONE.
    ///
    /// This is a primitive used by collapseEdge after twin-stitching.
    /// Direct callers must understand it leaves no twin link across the
    /// removed face — typically you stitch the relevant wings yourself
    /// (rewriting twin pointers on the surviving sides) before calling.
    void deleteFace(int32_t faceIdx);

    /// Walk every unpaired half-edge and pair with any other unpaired he
    /// whose endpoints match by vertex identity, then by quantized position.
    void rematchTwins();

    /// Collapse every set of triangular faces sharing a `group` tag into a
    /// single N-gon face per connected component of that group. Boundary
    /// of each component is walked as an ordered loop; the N-gon face uses
    /// that loop as its vertex order. Preserves the group id on the new
    /// face. After merging, internal half-edges of the merged tris are
    /// removed; external twins are re-linked.
    ///
    /// Use this immediately after fromMeshData(positions, indices, triToGroup)
    /// to lift a triangle-soup mesh into proper N-gon face semantics.
    void mergeFacesByGroup();

    /// Drop dead half-edges, dead faces, and unreferenced vertices, and
    /// renumber survivors. After compact(), every index in the mesh is in
    /// [0, count) and there are no tombstones. Surviving relative ordering
    /// is preserved within each array.
    void compact();

private:
    std::vector<Vertex>   vertices_;
    std::vector<HalfEdge> halfEdges_;
    std::vector<Face>     faces_;
};

} // namespace bromesh
