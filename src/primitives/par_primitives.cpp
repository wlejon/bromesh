#include "bromesh/primitives/par_primitives.h"
#include "bromesh/manipulation/normals.h"

#ifdef BROMESH_HAS_PAR_SHAPES
#include <par_shapes.h>
#include <cmath>
#endif

namespace bromesh {

#ifdef BROMESH_HAS_PAR_SHAPES

// Convert a par_shapes_mesh to MeshData, optionally scaling positions.
static MeshData fromParShapes(par_shapes_mesh* pm, float scale) {
    if (!pm) return {};

    MeshData m;
    size_t vc = pm->npoints;
    size_t tc = pm->ntriangles;

    m.positions.resize(vc * 3);
    for (size_t i = 0; i < vc * 3; ++i) {
        m.positions[i] = pm->points[i] * scale;
    }

    if (pm->normals) {
        m.normals.assign(pm->normals, pm->normals + vc * 3);
    }

    if (pm->tcoords) {
        m.uvs.assign(pm->tcoords, pm->tcoords + vc * 2);
    }

    m.indices.resize(tc * 3);
    for (size_t i = 0; i < tc * 3; ++i) {
        m.indices[i] = static_cast<uint32_t>(pm->triangles[i]);
    }

    par_shapes_free_mesh(pm);
    return m;
}

#endif // BROMESH_HAS_PAR_SHAPES

MeshData geodesicSphere(float radius, int nsubdivisions) {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_subdivided_sphere(nsubdivisions), radius);
#else
    (void)radius; (void)nsubdivisions;
    return {};
#endif
}

MeshData icosahedron() {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_icosahedron(), 1.0f);
#else
    return {};
#endif
}

MeshData dodecahedron() {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_dodecahedron(), 1.0f);
#else
    return {};
#endif
}

MeshData octahedron() {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_octahedron(), 1.0f);
#else
    return {};
#endif
}

MeshData tetrahedron() {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_tetrahedron(), 1.0f);
#else
    return {};
#endif
}

MeshData cone(float radius, float height, int slices, int stacks, bool capBase) {
#ifdef BROMESH_HAS_PAR_SHAPES
    auto* pm = par_shapes_create_cone(slices, stacks);
    if (!pm) return {};

    // par_shapes cone: X,Y are radial in [-1,1], Z is height in [0,1].
    // Remap to bro convention: base disc at Y=0 with `radius`, apex at
    // Y=`height`. par_shapes' radial pair (X,Y) maps to (X,Z); par_shapes'
    // Z (height) maps to Y.
    for (int i = 0; i < pm->npoints; ++i) {
        float pX = pm->points[i * 3 + 0]; // radial
        float pY = pm->points[i * 3 + 1]; // radial
        float pZ = pm->points[i * 3 + 2]; // height [0,1]
        pm->points[i * 3 + 0] = pX * radius;
        pm->points[i * 3 + 1] = pZ * height;
        pm->points[i * 3 + 2] = pY * radius;
    }

    // Recompute normals after transform
    par_shapes_compute_normals(pm);

    MeshData m = fromParShapes(pm, 1.0f);
    if (!capBase || m.empty()) return m;

    // Append a fan cap at Y=0: one center vertex + `slices` ring vertices,
    // wound so triangle normals point -Y. Cap UVs trace the unit disc
    // (cos*0.5+0.5, sin*0.5+0.5). The seam at the lateral/base corner is
    // intentional — we don't want a smooth normal across the crease.
    const bool hasN = !m.normals.empty();
    const bool hasU = !m.uvs.empty();
    const uint32_t centerIdx = (uint32_t)m.vertexCount();

    m.positions.insert(m.positions.end(), { 0.0f, 0.0f, 0.0f });
    if (hasN) m.normals.insert(m.normals.end(), { 0.0f, -1.0f, 0.0f });
    if (hasU) m.uvs.insert(m.uvs.end(), { 0.5f, 0.5f });

    constexpr float kTwoPi = 6.28318530717958647692f;
    for (int i = 0; i < slices; ++i) {
        float a = kTwoPi * (float)i / (float)slices;
        float c = std::cos(a), s = std::sin(a);
        m.positions.insert(m.positions.end(), { c * radius, 0.0f, s * radius });
        if (hasN) m.normals.insert(m.normals.end(), { 0.0f, -1.0f, 0.0f });
        if (hasU) m.uvs.insert(m.uvs.end(), { 0.5f + c * 0.5f, 0.5f + s * 0.5f });
    }
    for (int i = 0; i < slices; ++i) {
        uint32_t a = centerIdx + 1 + (uint32_t)i;
        uint32_t b = centerIdx + 1 + (uint32_t)((i + 1) % slices);
        m.indices.push_back(centerIdx);
        m.indices.push_back(a);
        m.indices.push_back(b);
    }
    return m;
#else
    (void)radius; (void)height; (void)slices; (void)stacks; (void)capBase;
    return {};
#endif
}

MeshData disc(float radius, int slices) {
#ifdef BROMESH_HAS_PAR_SHAPES
    // par_shapes emits the disk in the XY plane with normal +Z. Bro
    // convention is Y-up: leafCard, cone-base disc, etc. all live in the
    // XZ plane with normal +Y. Swap Y/Z so callers don't have to rotate.
    MeshData m = fromParShapes(par_shapes_create_parametric_disk(slices, 1), radius);
    if (m.empty()) return m;
    for (size_t i = 0; i < m.vertexCount(); ++i) {
        std::swap(m.positions[i * 3 + 1], m.positions[i * 3 + 2]);
    }
    if (!m.normals.empty()) {
        for (size_t i = 0; i < m.vertexCount(); ++i) {
            m.normals[i * 3 + 0] = 0.0f;
            m.normals[i * 3 + 1] = 1.0f;
            m.normals[i * 3 + 2] = 0.0f;
        }
    }
    // Swapping Y/Z flips the triangle winding from +Z facing to -Y facing.
    // Reverse each triangle so it faces +Y (matching the new normal).
    for (size_t t = 0; t < m.triangleCount(); ++t) {
        std::swap(m.indices[t * 3 + 1], m.indices[t * 3 + 2]);
    }
    return m;
#else
    (void)radius; (void)slices;
    return {};
#endif
}

MeshData rock(float radius, int seed, int nsubdivisions) {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_rock(seed, nsubdivisions), radius);
#else
    (void)radius; (void)seed; (void)nsubdivisions;
    return {};
#endif
}

MeshData blob(float radius, int seed, int nsubdivisions,
              float scaleX, float scaleY, float scaleZ,
              float centerX, float centerY, float centerZ) {
    MeshData m = rock(radius, seed, nsubdivisions);
    if (m.empty()) return m;
    const bool nonUniform =
        scaleX != 1.0f || scaleY != 1.0f || scaleZ != 1.0f;
    for (size_t i = 0; i < m.vertexCount(); ++i) {
        m.positions[i * 3 + 0] = m.positions[i * 3 + 0] * scaleX + centerX;
        m.positions[i * 3 + 1] = m.positions[i * 3 + 1] * scaleY + centerY;
        m.positions[i * 3 + 2] = m.positions[i * 3 + 2] * scaleZ + centerZ;
    }
    // Non-uniform scale invalidates incoming normals; uniform scale + pure
    // translation preserve them.
    if (nonUniform) computeNormals(m);
    return m;
}

MeshData trefoilKnot(float radius, int slices, int stacks) {
#ifdef BROMESH_HAS_PAR_SHAPES
    auto* pm = par_shapes_create_trefoil_knot(slices, stacks, radius);
    if (!pm) return {};
    par_shapes_compute_normals(pm);
    return fromParShapes(pm, 1.0f);
#else
    (void)radius; (void)slices; (void)stacks;
    return {};
#endif
}

MeshData kleinBottle(int slices, int stacks) {
#ifdef BROMESH_HAS_PAR_SHAPES
    auto* pm = par_shapes_create_klein_bottle(slices, stacks);
    if (!pm) return {};
    par_shapes_compute_normals(pm);
    return fromParShapes(pm, 1.0f);
#else
    (void)slices; (void)stacks;
    return {};
#endif
}

} // namespace bromesh
