#include "test_framework.h"
#include "bromesh/manipulation/lathe.h"
#include "bromesh/analysis/bbox.h"

#include <cmath>
#include <vector>

TEST(lathe_cylinder) {
    float r = 1.5f;
    float h = 3.0f;
    std::vector<bromath::Vec2> profile = { { r, 0.0f }, { r, h } };

    bromesh::LatheOptions opts;
    opts.segments = 32;
    opts.axis = 1; // Y axis

    auto m = bromesh::lathe(profile, opts);
    ASSERT(!m.indices.empty(), "lathe_cylinder: mesh is non-empty");
    // 32 segments * 2 triangles per quad = 64 triangles
    ASSERT(m.triangleCount() == 64, "lathe_cylinder: triangle count is 64");

    auto bbox = bromesh::computeBBox(m);
    ASSERT(std::fabs((bbox.max.y - bbox.min.y) - h) < 1e-4f, "lathe_cylinder: height matches 3.0");

    for (size_t i = 0; i < m.vertexCount(); ++i) {
        float x = m.positions[i * 3 + 0];
        float z = m.positions[i * 3 + 2];
        float dist = std::sqrt(x * x + z * z);
        ASSERT(std::fabs(dist - r) < 1e-4f, "lathe_cylinder: vertex distance matches radius 1.5");
    }
}

TEST(lathe_sphere) {
    float radius = 2.0f;
    int profileSubdiv = 32;
    std::vector<bromath::Vec2> profile;
    profile.reserve(profileSubdiv + 1);
    constexpr float kPi = 3.14159265358979323846f;

    for (int i = 0; i <= profileSubdiv; ++i) {
        float phi = -kPi * 0.5f + kPi * (static_cast<float>(i) / static_cast<float>(profileSubdiv));
        float r = radius * std::cos(phi);
        float h = radius * std::sin(phi);
        profile.push_back({ r, h });
    }

    bromesh::LatheOptions opts;
    opts.segments = 48;
    opts.axis = 1;

    auto m = bromesh::lathe(profile, opts);
    ASSERT(bromesh::isManifold(m), "lathe_sphere: sphere is manifold");

    float vol = bromesh::computeVolume(m);
    float analyticVol = (4.0f / 3.0f) * kPi * radius * radius * radius;
    float err = std::fabs(vol - analyticVol) / analyticVol;
    ASSERT(err < 0.02f, "lathe_sphere: volume matches analytic 4/3*pi*r^3 within 2%");
}

TEST(lathe_partial) {
    constexpr float kPi = 3.14159265358979323846f;

    // Semicircle revolved 180 degrees with caps
    float radius = 1.5f;
    int profileSubdiv = 24;
    std::vector<bromath::Vec2> profile;
    profile.reserve(profileSubdiv + 1);
    for (int i = 0; i <= profileSubdiv; ++i) {
        float phi = -kPi * 0.5f + kPi * (static_cast<float>(i) / static_cast<float>(profileSubdiv));
        float r = radius * std::cos(phi);
        float h = radius * std::sin(phi);
        profile.push_back({ r, h });
    }

    bromesh::LatheOptions opts;
    opts.segments = 32;
    opts.startAngle = 0.0f;
    opts.endAngle = kPi;
    opts.capStart = true;
    opts.capEnd = true;
    opts.axis = 1;

    auto m = bromesh::lathe(profile, opts);
    ASSERT(bromesh::isManifold(m), "lathe_partial: hemisphere with caps is manifold");

    float vol = bromesh::computeVolume(m);
    float expectedVol = 0.5f * (4.0f / 3.0f) * kPi * radius * radius * radius;
    float err = std::fabs(vol - expectedVol) / expectedVol;
    ASSERT(err < 0.02f, "lathe_partial: hemisphere volume matches half sphere within 2%");

    // Solid box profile revolved 180 degrees with caps
    std::vector<bromath::Vec2> boxProf = {
        { 0.0f, 0.0f }, { 1.0f, 0.0f }, { 1.0f, 2.0f }, { 0.0f, 2.0f }
    };
    auto m2 = bromesh::lathe(boxProf, opts);
    ASSERT(bromesh::isManifold(m2), "lathe_partial: half-cylinder with caps is manifold");

    float vol2 = bromesh::computeVolume(m2);
    float expectedVol2 = 0.5f * kPi * 1.0f * 1.0f * 2.0f;
    float err2 = std::fabs(vol2 - expectedVol2) / expectedVol2;
    ASSERT(err2 < 0.02f, "lathe_partial: half-cylinder volume matches analytic within 2%");
}
