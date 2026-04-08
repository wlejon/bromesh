#include "bromesh/primitives/par_primitives.h"

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

MeshData cone(float radius, float height, int slices, int stacks) {
#ifdef BROMESH_HAS_PAR_SHAPES
    auto* pm = par_shapes_create_cone(slices, stacks);
    if (!pm) return {};

    // par_shapes cone is unit-sized: radius 1, height 1 along Z.
    // Scale and orient to our convention (along Y axis).
    for (int i = 0; i < pm->npoints; ++i) {
        float x = pm->points[i * 3 + 0] * radius;
        float z_orig = pm->points[i * 3 + 2] * radius;
        float y = pm->points[i * 3 + 1] * height;
        pm->points[i * 3 + 0] = x;
        pm->points[i * 3 + 1] = y;
        pm->points[i * 3 + 2] = z_orig;
    }

    // Recompute normals after transform
    par_shapes_compute_normals(pm);

    return fromParShapes(pm, 1.0f);
#else
    (void)radius; (void)height; (void)slices; (void)stacks;
    return {};
#endif
}

MeshData disc(float radius, int slices) {
#ifdef BROMESH_HAS_PAR_SHAPES
    return fromParShapes(par_shapes_create_parametric_disk(slices, 1), radius);
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
