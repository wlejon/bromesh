#include "test_framework.h"
#include <cmath>

TEST(normalize_weights_basic) {
    bromesh::SkinData skin;
    // 1 vertex, 4 weights that don't sum to 1
    skin.boneWeights = {0.5f, 0.3f, 0.1f, 0.05f};
    skin.boneIndices = {0, 1, 2, 3};
    skin.inverseBindMatrices.resize(4 * 16, 0.0f);
    skin.boneCount = 4;

    bromesh::normalizeWeights(skin);

    float sum = skin.boneWeights[0] + skin.boneWeights[1] +
                skin.boneWeights[2] + skin.boneWeights[3];
    ASSERT(std::fabs(sum - 1.0f) < 0.001f, "normalize: weights should sum to 1");
    // Weights should be sorted descending
    ASSERT(skin.boneWeights[0] >= skin.boneWeights[1], "normalize: sorted descending");
}

TEST(normalize_weights_zeros) {
    bromesh::SkinData skin;
    skin.boneWeights = {0.0f, 0.0f, 0.0f, 0.0f};
    skin.boneIndices = {0, 1, 2, 3};
    skin.boneCount = 4;

    bromesh::normalizeWeights(skin);

    // Should default to bone 0 with weight 1
    ASSERT(std::fabs(skin.boneWeights[0] - 1.0f) < 0.001f,
           "normalize_zeros: first weight should be 1");
    ASSERT(skin.boneIndices[0] == 0, "normalize_zeros: first index should be 0");
}

TEST(apply_skinning_identity) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    size_t vCount = mesh.vertexCount();

    bromesh::SkinData skin;
    skin.boneCount = 1;
    // Identity inverse bind matrix
    skin.inverseBindMatrices = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };
    skin.boneWeights.resize(vCount * 4, 0.0f);
    skin.boneIndices.resize(vCount * 4, 0);
    for (size_t v = 0; v < vCount; ++v) {
        skin.boneWeights[v * 4] = 1.0f;
    }

    // Save original positions
    auto origPos = mesh.positions;

    // Identity pose matrix
    float pose[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    bromesh::applySkinning(mesh, skin, pose);

    // Positions should be unchanged
    bool same = true;
    for (size_t i = 0; i < origPos.size(); ++i) {
        if (std::fabs(mesh.positions[i] - origPos[i]) > 0.001f) {
            same = false; break;
        }
    }
    ASSERT(same, "skin_identity: positions unchanged with identity transform");
}

TEST(apply_skinning_translation) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    size_t vCount = mesh.vertexCount();

    bromesh::SkinData skin;
    skin.boneCount = 1;
    skin.inverseBindMatrices = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };
    skin.boneWeights.resize(vCount * 4, 0.0f);
    skin.boneIndices.resize(vCount * 4, 0);
    for (size_t v = 0; v < vCount; ++v)
        skin.boneWeights[v * 4] = 1.0f;

    auto origPos = mesh.positions;

    // Translate by (5, 0, 0) via pose matrix (column-major)
    float pose[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 5,0,0,1};
    bromesh::applySkinning(mesh, skin, pose);

    bool shifted = true;
    for (size_t v = 0; v < vCount; ++v) {
        if (std::fabs(mesh.positions[v*3+0] - (origPos[v*3+0] + 5.0f)) > 0.001f) {
            shifted = false; break;
        }
    }
    ASSERT(shifted, "skin_translate: positions shifted by +5 on X");
}

TEST(apply_morph_target) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;
    size_t vCount = mesh.vertexCount();

    bromesh::MorphTarget morph;
    morph.name = "expand";
    morph.deltaPositions.resize(vCount * 3, 0.0f);
    // Move every vertex +1 on X
    for (size_t v = 0; v < vCount; ++v)
        morph.deltaPositions[v * 3 + 0] = 1.0f;

    bromesh::applyMorphTarget(mesh, morph, 0.5f);

    bool correct = true;
    for (size_t v = 0; v < vCount; ++v) {
        if (std::fabs(mesh.positions[v*3+0] - (origPos[v*3+0] + 0.5f)) > 0.001f) {
            correct = false; break;
        }
    }
    ASSERT(correct, "morph: positions shifted +0.5 on X at weight 0.5");
}

TEST(apply_morph_zero_weight) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;

    bromesh::MorphTarget morph;
    morph.deltaPositions.resize(mesh.vertexCount() * 3, 100.0f);

    bromesh::applyMorphTarget(mesh, morph, 0.0f);

    ASSERT(mesh.positions == origPos, "morph_zero: no change at weight 0");
}

TEST(simplify_preserves_bone_weights_and_indices) {
    auto mesh = bromesh::sphere(1.0f, 16, 16);
    const size_t vCount = mesh.vertexCount();
    mesh.boneWeights.resize(vCount * 4, 0.0f);
    mesh.boneIndices.resize(vCount * 4, 0);

    for (size_t v = 0; v < vCount; ++v) {
        float y = mesh.positions[v * 3 + 1];
        if (y > 0.0f) {
            mesh.boneWeights[v * 4 + 0] = 0.75f;
            mesh.boneWeights[v * 4 + 1] = 0.25f;
            mesh.boneIndices[v * 4 + 0] = 1;
            mesh.boneIndices[v * 4 + 1] = 2;
        } else {
            mesh.boneWeights[v * 4 + 0] = 1.0f;
            mesh.boneIndices[v * 4 + 0] = 0;
        }
    }
    ASSERT(mesh.hasBoneWeights(), "mesh has boneWeights");
    ASSERT(mesh.hasBoneIndices(), "mesh has boneIndices");

    auto simplified = bromesh::simplify(mesh, 0.5f);
    ASSERT(simplified.triangleCount() < mesh.triangleCount(), "simplified drops triangles");
    ASSERT(simplified.hasBoneWeights(), "simplified preserves boneWeights");
    ASSERT(simplified.hasBoneIndices(), "simplified preserves boneIndices");
    ASSERT(simplified.boneWeights.size() == simplified.vertexCount() * 4, "weights properly sized");
    ASSERT(simplified.boneIndices.size() == simplified.vertexCount() * 4, "indices properly sized");

    for (size_t v = 0; v < simplified.vertexCount(); ++v) {
        float y = simplified.positions[v * 3 + 1];
        if (y > 0.0f) {
            ASSERT(std::fabs(simplified.boneWeights[v * 4 + 0] - 0.75f) < 1e-4f, "weight 0 preserved");
            ASSERT(std::fabs(simplified.boneWeights[v * 4 + 1] - 0.25f) < 1e-4f, "weight 1 preserved");
            ASSERT(simplified.boneIndices[v * 4 + 0] == 1, "index 0 preserved");
            ASSERT(simplified.boneIndices[v * 4 + 1] == 2, "index 1 preserved");
        } else {
            ASSERT(std::fabs(simplified.boneWeights[v * 4 + 0] - 1.0f) < 1e-4f, "weight 0 preserved");
            ASSERT(simplified.boneIndices[v * 4 + 0] == 0, "index 0 preserved");
        }
    }
}

TEST(validate_skin_standalone) {
    bromesh::SkinData skin;
    skin.boneWeights = { 0.6f, 0.4f, 0.0f, 0.0f,  1.0f, 0.0f, 0.0f, 0.0f };
    skin.boneIndices = { 0, 1, 0, 0,  2, 0, 0, 0 };
    skin.boneCount = 3;

    auto val = bromesh::validateSkin(skin);
    ASSERT(val.clean(), "valid skin reports clean");
    ASSERT(val.vertexCount == 2, "skin reports 2 vertices");
    ASSERT(val.maxInfluencesObserved == 2, "max influences 2");

    // Orphan vertex test
    skin.boneWeights[4] = 0.0f;
    auto val2 = bromesh::validateSkin(skin);
    ASSERT(!val2.clean(), "unweighted vertex is not clean");
    ASSERT(val2.orphanCount == 1, "orphan detected");
}

