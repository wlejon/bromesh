#include "bromesh/uv/unwrap.h"

#ifdef BROMESH_HAS_XATLAS
#include <xatlas.h>
#include <cstring>
#endif

namespace bromesh {

UnwrapResult unwrapUVs(MeshData& mesh,
                       const UnwrapParams& chartParams,
                       const PackParams& packParams) {
#ifdef BROMESH_HAS_XATLAS
    if (mesh.empty() || mesh.indices.empty()) return {};

    xatlas::Atlas* atlas = xatlas::Create();

    // Set up mesh declaration
    xatlas::MeshDecl meshDecl;
    meshDecl.vertexCount = static_cast<uint32_t>(mesh.vertexCount());
    meshDecl.vertexPositionData = mesh.positions.data();
    meshDecl.vertexPositionStride = sizeof(float) * 3;

    if (mesh.hasNormals()) {
        meshDecl.vertexNormalData = mesh.normals.data();
        meshDecl.vertexNormalStride = sizeof(float) * 3;
    }

    if (mesh.hasUVs()) {
        meshDecl.vertexUvData = mesh.uvs.data();
        meshDecl.vertexUvStride = sizeof(float) * 2;
    }

    meshDecl.indexCount = static_cast<uint32_t>(mesh.indices.size());
    meshDecl.indexData = mesh.indices.data();
    meshDecl.indexFormat = xatlas::IndexFormat::UInt32;

    xatlas::AddMeshError error = xatlas::AddMesh(atlas, meshDecl);
    if (error != xatlas::AddMeshError::Success) {
        xatlas::Destroy(atlas);
        return {};
    }

    // Configure chart options
    xatlas::ChartOptions cOpts;
    cOpts.maxChartArea = 0; // no limit
    cOpts.maxBoundaryLength = 0;
    cOpts.normalDeviationWeight = chartParams.normalDeviationWeight;
    cOpts.roundnessWeight = chartParams.roundnessWeight;
    cOpts.straightnessWeight = chartParams.straightnessWeight;
    cOpts.normalSeamWeight = chartParams.normalSeamWeight;
    cOpts.textureSeamWeight = chartParams.textureSeamWeight;
    cOpts.maxCost = chartParams.maxStretch;
    cOpts.maxIterations = 1;

    // Configure pack options
    xatlas::PackOptions pOpts;
    pOpts.resolution = packParams.resolution;
    pOpts.padding = packParams.padding;
    pOpts.bilinear = packParams.bilinear;
    pOpts.blockAlign = packParams.blockAlign;
    pOpts.bruteForce = packParams.bruteForce;

    xatlas::Generate(atlas, cOpts, pOpts);

    if (atlas->meshCount == 0 || atlas->meshes[0].vertexCount == 0) {
        xatlas::Destroy(atlas);
        return {};
    }

    // Extract result: xatlas may have split vertices at UV seams
    const xatlas::Mesh& outMesh = atlas->meshes[0];
    uint32_t newVertCount = outMesh.vertexCount;
    uint32_t newIdxCount = outMesh.indexCount;

    MeshData result;
    result.positions.resize(newVertCount * 3);
    if (mesh.hasNormals()) result.normals.resize(newVertCount * 3);
    if (mesh.hasColors()) result.colors.resize(newVertCount * 4);
    result.uvs.resize(newVertCount * 2);
    result.indices.resize(newIdxCount);

    float atlasW = static_cast<float>(atlas->width);
    float atlasH = static_cast<float>(atlas->height);
    if (atlasW < 1.0f) atlasW = 1.0f;
    if (atlasH < 1.0f) atlasH = 1.0f;

    for (uint32_t v = 0; v < newVertCount; ++v) {
        const xatlas::Vertex& xv = outMesh.vertexArray[v];
        uint32_t origIdx = xv.xref; // index into original vertex buffer

        result.positions[v * 3 + 0] = mesh.positions[origIdx * 3 + 0];
        result.positions[v * 3 + 1] = mesh.positions[origIdx * 3 + 1];
        result.positions[v * 3 + 2] = mesh.positions[origIdx * 3 + 2];

        if (mesh.hasNormals()) {
            result.normals[v * 3 + 0] = mesh.normals[origIdx * 3 + 0];
            result.normals[v * 3 + 1] = mesh.normals[origIdx * 3 + 1];
            result.normals[v * 3 + 2] = mesh.normals[origIdx * 3 + 2];
        }

        if (mesh.hasColors()) {
            result.colors[v * 4 + 0] = mesh.colors[origIdx * 4 + 0];
            result.colors[v * 4 + 1] = mesh.colors[origIdx * 4 + 1];
            result.colors[v * 4 + 2] = mesh.colors[origIdx * 4 + 2];
            result.colors[v * 4 + 3] = mesh.colors[origIdx * 4 + 3];
        }

        // Normalize UVs to [0,1] range
        result.uvs[v * 2 + 0] = xv.uv[0] / atlasW;
        result.uvs[v * 2 + 1] = xv.uv[1] / atlasH;
    }

    for (uint32_t i = 0; i < newIdxCount; ++i) {
        result.indices[i] = outMesh.indexArray[i];
    }

    UnwrapResult res;
    res.atlasWidth = static_cast<int>(atlas->width);
    res.atlasHeight = static_cast<int>(atlas->height);
    res.chartCount = static_cast<int>(atlas->chartCount);
    res.success = true;

    xatlas::Destroy(atlas);

    mesh = std::move(result);
    return res;
#else
    (void)mesh; (void)chartParams; (void)packParams;
    return {};
#endif
}

} // namespace bromesh
