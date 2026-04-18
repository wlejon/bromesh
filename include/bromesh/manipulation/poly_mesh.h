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

    /// Drop unreferenced vertices and renumber. Preserves face/he ordering.
    void compact();

private:
    std::vector<Vertex>   vertices_;
    std::vector<HalfEdge> halfEdges_;
    std::vector<Face>     faces_;
};

} // namespace bromesh
