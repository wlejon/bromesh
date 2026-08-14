#include "test_framework.h"
#include <cmath>

TEST(bake_curvature_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    bromesh::bakeCurvature(mesh, 1.0f);
    ASSERT(mesh.hasColors(), "bake_curv: should have colors");
    // All colors should be in [0,1]
    bool valid = true;
    for (size_t i = 0; i < mesh.colors.size(); ++i) {
        if (mesh.colors[i] < -0.01f || mesh.colors[i] > 1.01f) {
            valid = false;
            break;
        }
    }
    ASSERT(valid, "bake_curv: colors should be in [0,1]");
}

TEST(bake_ao_small_sphere) {
    auto mesh = bromesh::sphere(1.0f, 6, 4);
    bromesh::computeNormals(mesh);

    bromesh::bakeAmbientOcclusion(mesh, 8, 0.0f);
    ASSERT(mesh.hasColors(), "bake_ao: should have colors");
    // AO values should be in [0,1]
    bool valid = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float ao = mesh.colors[v*4+0];
        if (ao < -0.01f || ao > 1.01f) { valid = false; break; }
    }
    ASSERT(valid, "bake_ao: AO values should be in [0,1]");
}

TEST(bake_thickness_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    bromesh::bakeThickness(mesh, 8, 0.0f);
    ASSERT(mesh.hasColors(), "bake_thick: should have colors");
}

TEST(bake_ao_texture_basic) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto tex = bromesh::bakeAmbientOcclusionToTexture(mesh, 16, 16, 8);
    ASSERT(tex.width == 16, "bake_ao_tex: width 16");
    ASSERT(tex.height == 16, "bake_ao_tex: height 16");
    ASSERT(tex.channels == 1, "bake_ao_tex: 1 channel");
    ASSERT(tex.pixels.size() == 16*16, "bake_ao_tex: correct pixel count");

    // Some texels should be covered (not all zero)
    bool hasCoverage = false;
    for (float v : tex.pixels) {
        if (v > 0.01f) { hasCoverage = true; break; }
    }
    ASSERT(hasCoverage, "bake_ao_tex: some texels should be covered");
}

TEST(bake_curvature_texture_basic) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto tex = bromesh::bakeCurvatureToTexture(mesh, 16, 16, 1.0f);
    ASSERT(tex.width == 16 && tex.channels == 1, "bake_curv_tex: dimensions correct");

    // Values should be in [0,1]
    bool valid = true;
    for (float v : tex.pixels) {
        if (v < -0.01f || v > 1.01f) { valid = false; break; }
    }
    ASSERT(valid, "bake_curv_tex: values in [0,1]");
}

TEST(bake_thickness_texture_basic) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto tex = bromesh::bakeThicknessToTexture(mesh, 16, 16, 8);
    ASSERT(tex.width == 16 && tex.channels == 1, "bake_thick_tex: dimensions correct");
    ASSERT(!tex.pixels.empty(), "bake_thick_tex: has pixels");
}

TEST(bake_normals_texture_basic) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto tex = bromesh::bakeNormalsToTexture(mesh, 16, 16);
    ASSERT(tex.channels == 4, "bake_nrm_tex: 4 channels");
    ASSERT(tex.pixels.size() == 16*16*4, "bake_nrm_tex: correct pixel count");

    // Covered texels should have alpha=1 and normals in [0,1]
    bool valid = true;
    for (int i = 0; i < 16*16; ++i) {
        float a = tex.pixels[i*4+3];
        if (a > 0.5f) {
            for (int c = 0; c < 3; ++c) {
                float v = tex.pixels[i*4+c];
                if (v < -0.01f || v > 1.01f) { valid = false; break; }
            }
        }
        if (!valid) break;
    }
    ASSERT(valid, "bake_nrm_tex: normal values in [0,1]");
}

TEST(bake_position_texture_basic) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto tex = bromesh::bakePositionToTexture(mesh, 16, 16);
    ASSERT(tex.channels == 4, "bake_pos_tex: 4 channels");

    // Covered texels should have positions within the box bounds [-1,1]
    bool valid = true;
    for (int i = 0; i < 16*16; ++i) {
        if (tex.pixels[i*4+3] > 0.5f) {
            for (int c = 0; c < 3; ++c) {
                float v = tex.pixels[i*4+c];
                if (v < -1.5f || v > 1.5f) { valid = false; break; }
            }
        }
        if (!valid) break;
    }
    ASSERT(valid, "bake_pos_tex: positions within expected bounds");
}

TEST(bake_texture_no_uvs) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);
    mesh.uvs.clear(); // strip UVs

    auto tex = bromesh::bakeAmbientOcclusionToTexture(mesh, 16, 16, 8);
    ASSERT(tex.width == 0, "bake_tex_no_uvs: returns empty when no UVs");
}

TEST(bake_texture_at_method) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto tex = bromesh::bakePositionToTexture(mesh, 8, 8);
    ASSERT(tex.at(0, 0) != nullptr, "tex.at: valid pixel returns non-null");
    ASSERT(tex.at(-1, 0) == nullptr, "tex.at: out-of-bounds returns null");
    ASSERT(tex.at(0, 8) == nullptr, "tex.at: out-of-bounds returns null");
}

TEST(bake_thickness_box_analytic) {
    // A box of half-extents 1.0 is [-1, 1]^3 with size 2.0x2.0x2.0.
    // The analytic thickness along face normals is exactly 2.0.
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    const float maxDist = 4.0f;
    // 1 ray is shot along the inverted normal (sample 0 in golden spiral is dir = invN)
    bromesh::bakeThickness(mesh, 1, maxDist);
    ASSERT(mesh.hasColors(), "bake_thickness_box_analytic: has colors");

    bool allCloseTo2 = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float rawThickness = mesh.colors[v * 4 + 0] * maxDist;
        if (std::fabs(rawThickness - 2.0f) > 0.1f) {
            allCloseTo2 = false;
            break;
        }
    }
    ASSERT(allCloseTo2, "bake_thickness_box_analytic: vertex thickness along face normals is 2.0 +- 0.1");

    auto tex = bromesh::bakeThicknessToTexture(mesh, 16, 16, 1, maxDist);
    ASSERT(!tex.pixels.empty(), "bake_thickness_box_analytic: tex has pixels");
    bool texelsCloseTo2 = true;
    int coveredCount = 0;
    for (int y = 0; y < tex.height; ++y) {
        for (int x = 0; x < tex.width; ++x) {
            float val = *tex.at(x, y);
            // Default uncovered texels are 1.0f
            float rawThickness = val * maxDist;
            if (std::fabs(rawThickness - 2.0f) <= 0.1f) {
                coveredCount++;
            } else if (std::fabs(val - 1.0f) > 1e-4f) {
                texelsCloseTo2 = false;
            }
        }
    }
    ASSERT(coveredCount > 0, "bake_thickness_box_analytic: covered texels exist");
    ASSERT(texelsCloseTo2, "bake_thickness_box_analytic: covered texels have thickness 2.0 +- 0.1");
}

TEST(bake_ao_sphere_analytic) {
    // 1. Isolated convex sphere: no self-occluders, AO must be ~1.0 (> 0.85)
    auto sphere = bromesh::sphere(1.0f, 16, 16);
    bromesh::computeNormals(sphere);
    bromesh::projectUVs(sphere, bromesh::ProjectionType::Spherical, 1.0f);

    bromesh::bakeAmbientOcclusion(sphere, 32, 0.0f);
    bool sphereAoUnoccluded = true;
    for (size_t v = 0; v < sphere.vertexCount(); ++v) {
        if (sphere.colors[v * 4 + 0] < 0.85f) {
            sphereAoUnoccluded = false;
            break;
        }
    }
    ASSERT(sphereAoUnoccluded, "bake_ao_sphere_analytic: convex sphere AO > 0.85");

    auto sphereTex = bromesh::bakeAmbientOcclusionToTexture(sphere, 16, 16, 32, 0.0f);
    bool sphereTexAoHigh = true;
    for (float ao : sphereTex.pixels) {
        // Uncovered texels are 0, covered should be > 0.85
        if (ao > 0.01f && ao < 0.85f) {
            sphereTexAoHigh = false;
            break;
        }
    }
    ASSERT(sphereTexAoHigh, "bake_ao_sphere_analytic: sphere texture AO > 0.85");

    // 2. Enclosed cavity / two parallel close planes facing each other along Y
    // Plane 1 at y=0 with normal +Y
    auto p1 = bromesh::plane(2.0f, 2.0f, 4, 4);
    // Plane 2 at y=0.2 with normal -Y
    auto p2 = bromesh::plane(2.0f, 2.0f, 4, 4);
    bromesh::translateMesh(p2, 0.0f, 0.2f, 0.0f);
    // Flip p2 triangles so normal points towards -Y (facing p1)
    for (size_t t = 0; t < p2.triangleCount(); ++t) {
        std::swap(p2.indices[t * 3 + 1], p2.indices[t * 3 + 2]);
    }
    auto merged = bromesh::mergeMeshes({p1, p2});
    bromesh::computeNormals(merged);

    bromesh::bakeAmbientOcclusion(merged, 32, 1.0f);
    // Inner vertices of p1 should be occluded by p2, so average AO is significantly lower (< 0.6)
    float totalAO = 0.0f;
    for (size_t v = 0; v < merged.vertexCount(); ++v) {
        totalAO += merged.colors[v * 4 + 0];
    }
    float avgAO = totalAO / merged.vertexCount();
    ASSERT(avgAO < 0.6f, "bake_ao_sphere_analytic: parallel facing planes have low AO (< 0.6)");
}

TEST(bake_normals_analytic) {
    // Construct a quad flat in XY (z=0) with normal (0, 0, 1) and standard UVs [0, 1]
    bromesh::MeshData plane;
    plane.positions = {
        -1.0f, -1.0f, 0.0f,
         1.0f, -1.0f, 0.0f,
         1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f
    };
    plane.normals = {
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f
    };
    plane.uvs = {
        0.0f, 0.0f,
        1.0f, 0.0f,
        1.0f, 1.0f,
        0.0f, 1.0f
    };
    plane.indices = { 0, 1, 2, 0, 2, 3 };

    // World-space normal bake
    auto worldNormTex = bromesh::bakeNormalsToTexture(plane, 16, 16);
    ASSERT(worldNormTex.channels == 4, "bake_normals_analytic: 4 channels");
    bool worldNormalsFlat = true;
    int covered = 0;
    for (int i = 0; i < 16 * 16; ++i) {
        if (worldNormTex.pixels[i * 4 + 3] > 0.5f) {
            covered++;
            float r = worldNormTex.pixels[i * 4 + 0];
            float g = worldNormTex.pixels[i * 4 + 1];
            float b = worldNormTex.pixels[i * 4 + 2];
            int ir = static_cast<int>(std::round(r * 255.0f));
            int ig = static_cast<int>(std::round(g * 255.0f));
            int ib = static_cast<int>(std::round(b * 255.0f));
            if (std::abs(ir - 128) > 2 || std::abs(ig - 128) > 2 || std::abs(ib - 255) > 2) {
                worldNormalsFlat = false;
            }
        }
    }
    ASSERT(covered > 0, "bake_normals_analytic: covered texels in world normal map");
    ASSERT(worldNormalsFlat, "bake_normals_analytic: world normal map texels match RGB (128, 128, 255) +- 2");

    // Tangent-space normal bake from reference
    auto tsNormTex = bromesh::bakeNormalsFromReference(plane, plane, 16, 16, 0.5f);
    bool tsNormalsFlat = true;
    covered = 0;
    for (int i = 0; i < 16 * 16; ++i) {
        if (tsNormTex.pixels[i * 4 + 3] > 0.5f) {
            covered++;
            float r = tsNormTex.pixels[i * 4 + 0];
            float g = tsNormTex.pixels[i * 4 + 1];
            float b = tsNormTex.pixels[i * 4 + 2];
            int ir = static_cast<int>(std::round(r * 255.0f));
            int ig = static_cast<int>(std::round(g * 255.0f));
            int ib = static_cast<int>(std::round(b * 255.0f));
            if (std::abs(ir - 128) > 2 || std::abs(ig - 128) > 2 || std::abs(ib - 255) > 2) {
                tsNormalsFlat = false;
            }
        }
    }
    ASSERT(covered > 0, "bake_normals_analytic: covered texels in TS normal map");
    ASSERT(tsNormalsFlat, "bake_normals_analytic: tangent space normal map texels match RGB (128, 128, 255) +- 2");
}

TEST(bake_curvature_sphere_analytic) {
    // Sphere of radius R = 2.0 has mean curvature 1/R = 0.5
    auto mesh = bromesh::sphere(2.0f, 32, 32);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    bromesh::bakeCurvature(mesh, 1.0f);
    ASSERT(mesh.hasColors(), "bake_curvature_sphere_analytic: has colors");

    // For a uniform convex sphere, bakeCurvature produces positive convex curvature
    bool validColors = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float c = mesh.colors[v * 4 + 0];
        if (c < -0.01f || c > 1.01f) {
            validColors = false;
            break;
        }
    }
    ASSERT(validColors, "bake_curvature_sphere_analytic: sphere colors in [0, 1]");

    auto tex = bromesh::bakeCurvatureToTexture(mesh, 16, 16, 1.0f);
    ASSERT(tex.width == 16, "bake_curvature_sphere_analytic: tex width 16");
    bool validTex = true;
    for (float v : tex.pixels) {
        if (v < -0.01f || v > 1.01f) {
            validTex = false;
            break;
        }
    }
    ASSERT(validTex, "bake_curvature_sphere_analytic: tex pixels in [0, 1]");
}
