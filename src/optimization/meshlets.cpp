#include "bromesh/optimization/meshlets.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#endif

namespace bromesh {

std::vector<Meshlet> buildMeshlets(const MeshData& mesh,
                                   const MeshletParams& params) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.empty() || mesh.indices.empty()) return {};

    size_t indexCount = mesh.indices.size();
    size_t vertexCount = mesh.vertexCount();
    size_t maxVerts = params.maxVertices;
    size_t maxTris = params.maxTriangles;

    size_t maxMeshlets = meshopt_buildMeshletsBound(indexCount, maxVerts, maxTris);

    std::vector<meshopt_Meshlet> rawMeshlets(maxMeshlets);
    std::vector<unsigned int> meshletVertices(maxMeshlets * maxVerts);
    std::vector<unsigned char> meshletTriangles(maxMeshlets * maxTris * 3);

    size_t meshletCount = meshopt_buildMeshlets(
        rawMeshlets.data(),
        meshletVertices.data(),
        meshletTriangles.data(),
        mesh.indices.data(),
        indexCount,
        mesh.positions.data(),
        vertexCount,
        sizeof(float) * 3,
        maxVerts,
        maxTris,
        params.coneWeight
    );

    std::vector<Meshlet> result;
    result.reserve(meshletCount);

    for (size_t i = 0; i < meshletCount; ++i) {
        const auto& raw = rawMeshlets[i];
        Meshlet m;

        // Copy vertex indices
        m.vertices.assign(
            meshletVertices.data() + raw.vertex_offset,
            meshletVertices.data() + raw.vertex_offset + raw.vertex_count
        );

        // Copy micro triangle indices
        m.triangles.assign(
            meshletTriangles.data() + raw.triangle_offset,
            meshletTriangles.data() + raw.triangle_offset + raw.triangle_count * 3
        );

        // Compute bounds for culling
        meshopt_Bounds bounds = meshopt_computeMeshletBounds(
            meshletVertices.data() + raw.vertex_offset,
            meshletTriangles.data() + raw.triangle_offset,
            raw.triangle_count,
            mesh.positions.data(),
            vertexCount,
            sizeof(float) * 3
        );

        m.bounds.center[0] = bounds.center[0];
        m.bounds.center[1] = bounds.center[1];
        m.bounds.center[2] = bounds.center[2];
        m.bounds.radius = bounds.radius;
        m.bounds.coneApex[0] = bounds.cone_apex[0];
        m.bounds.coneApex[1] = bounds.cone_apex[1];
        m.bounds.coneApex[2] = bounds.cone_apex[2];
        m.bounds.coneAxis[0] = bounds.cone_axis[0];
        m.bounds.coneAxis[1] = bounds.cone_axis[1];
        m.bounds.coneAxis[2] = bounds.cone_axis[2];
        m.bounds.coneCutoff = bounds.cone_cutoff;

        result.push_back(std::move(m));
    }

    return result;
#else
    (void)mesh; (void)params;
    return {};
#endif
}

} // namespace bromesh
