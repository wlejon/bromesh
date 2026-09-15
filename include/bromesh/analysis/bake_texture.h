#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/io/gltf.h"

#include <string>
#include <vector>

namespace bromesh {

/// A simple image buffer for baked results.
/// Pixels are stored row-major, bottom-to-top (OpenGL convention).
/// Each pixel has `channels` floats (1 for grayscale, 4 for RGBA).
struct TextureBuffer {
    std::vector<float> pixels;
    int width = 0;
    int height = 0;
    int channels = 0;

    /// Get a pointer to the pixel at (x, y), or nullptr if out of bounds.
    float* at(int x, int y) {
        if (x < 0 || x >= width || y < 0 || y >= height) return nullptr;
        return &pixels[(y * width + x) * channels];
    }
    const float* at(int x, int y) const {
        if (x < 0 || x >= width || y < 0 || y >= height) return nullptr;
        return &pixels[(y * width + x) * channels];
    }
};

/// Bake ambient occlusion into a texture using the mesh's UV layout.
/// The mesh must have UVs in [0,1]. Returns a single-channel (grayscale) texture.
/// White = unoccluded, black = fully occluded.
TextureBuffer bakeAmbientOcclusionToTexture(const MeshData& mesh,
                                             int texWidth, int texHeight,
                                             int numRays = 64,
                                             float maxDistance = 0.0f);

/// Bake mean curvature into a texture using the mesh's UV layout.
/// Returns a single-channel texture. 0.5 = flat, >0.5 = convex, <0.5 = concave.
TextureBuffer bakeCurvatureToTexture(const MeshData& mesh,
                                      int texWidth, int texHeight,
                                      float scale = 1.0f);

/// Bake mesh thickness into a texture using the mesh's UV layout.
/// Returns a single-channel texture. Brighter = thicker.
TextureBuffer bakeThicknessToTexture(const MeshData& mesh,
                                      int texWidth, int texHeight,
                                      int numRays = 32,
                                      float maxDistance = 0.0f);

/// Bake world-space normals into a texture using the mesh's UV layout.
/// Returns a 4-channel RGBA texture (xyz mapped to [0,1], alpha=1).
TextureBuffer bakeNormalsToTexture(const MeshData& mesh,
                                    int texWidth, int texHeight);

/// Bake a position map into a texture using the mesh's UV layout.
/// Returns a 4-channel RGBA texture (xyz = world position, alpha=1).
TextureBuffer bakePositionToTexture(const MeshData& mesh,
                                     int texWidth, int texHeight);

/// Save a TextureBuffer as an uncompressed TGA file.
/// Channels: 1 = 8-bit grayscale, 3 = 24-bit BGR, 4 = 32-bit BGRA.
/// Clamps float [0, 1] to uint8 [0, 255]. Handles coordinate orientation correctly
/// (TGA lower-left origin / descriptor bit 5 = 0).
bool saveImageTGA(const TextureBuffer& buf, const std::string& path);

/// Convert TextureBuffer (float row-major, bottom-to-top) to bromesh::Image
/// (RGBA8, top-left origin, stride width * 4), suitable for glTF materials or GPU upload.
/// If buf.channels == 1, replicates to RGB with A=255; if 3, sets A=255; if 4, copies RGBA.
Image textureToImage(const TextureBuffer& buf, const std::string& name = "");

} // namespace bromesh
