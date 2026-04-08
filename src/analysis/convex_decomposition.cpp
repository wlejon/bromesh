#include "bromesh/analysis/convex_decomposition.h"

#if BROMESH_HAS_VHACD
#include "VHACD.h"
#endif

namespace bromesh {

std::vector<MeshData> convexDecomposition(const MeshData& mesh,
                                          const ConvexDecompParams& params) {
#if BROMESH_HAS_VHACD
    if (mesh.empty()) return {};

    VHACD::IVHACD* vhacd = VHACD::CreateVHACD();
    if (!vhacd) return {};

    VHACD::IVHACD::Parameters vparams;
    vparams.m_maxConvexHulls = static_cast<uint32_t>(params.maxHulls);
    vparams.m_maxNumVerticesPerCH = static_cast<uint32_t>(params.maxVerticesPerHull);
    vparams.m_resolution = static_cast<uint32_t>(params.resolution);
    vparams.m_asyncACD = false;

    bool ok = vhacd->Compute(mesh.positions.data(),
                             static_cast<uint32_t>(mesh.vertexCount()),
                             mesh.indices.data(),
                             static_cast<uint32_t>(mesh.triangleCount()),
                             vparams);

    std::vector<MeshData> result;
    if (ok) {
        uint32_t nHulls = vhacd->GetNConvexHulls();
        result.reserve(nHulls);

        for (uint32_t i = 0; i < nHulls; ++i) {
            VHACD::IVHACD::ConvexHull hull;
            if (!vhacd->GetConvexHull(i, hull)) continue;

            MeshData hullMesh;
            hullMesh.positions.resize(hull.m_points.size() * 3);
            for (size_t v = 0; v < hull.m_points.size(); ++v) {
                hullMesh.positions[v * 3 + 0] = static_cast<float>(hull.m_points[v].mX);
                hullMesh.positions[v * 3 + 1] = static_cast<float>(hull.m_points[v].mY);
                hullMesh.positions[v * 3 + 2] = static_cast<float>(hull.m_points[v].mZ);
            }

            hullMesh.indices.resize(hull.m_triangles.size() * 3);
            for (size_t t = 0; t < hull.m_triangles.size(); ++t) {
                hullMesh.indices[t * 3 + 0] = hull.m_triangles[t].mI0;
                hullMesh.indices[t * 3 + 1] = hull.m_triangles[t].mI1;
                hullMesh.indices[t * 3 + 2] = hull.m_triangles[t].mI2;
            }

            result.push_back(std::move(hullMesh));
        }
    }

    vhacd->Release();
    return result;
#else
    (void)mesh; (void)params;
    return {};
#endif
}

MeshData convexHull(const MeshData& mesh) {
    // Single hull = convex decomposition with maxHulls=1
    ConvexDecompParams params;
    params.maxHulls = 1;
    auto hulls = convexDecomposition(mesh, params);
    if (hulls.empty()) return {};
    return std::move(hulls[0]);
}

} // namespace bromesh
