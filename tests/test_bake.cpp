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

