#pragma once

#include <bromath/vec.h>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace bromesh {

enum class SdfOp : uint8_t {
    Sphere,
    Box,
    RoundedBox,
    Cylinder,
    Capsule,
    Torus,
    Plane,
    Union,
    Intersection,
    Subtraction,
    SmoothUnion,
    SmoothIntersection,
    SmoothSubtraction,
    Translate,
    RotateY,
    ScaleUniform,
    Noise3D,
    Displace
};

struct SdfNode {
    SdfOp op = SdfOp::Sphere;
    int left = -1;
    int right = -1;
    bromath::Vec3 v0{0.0f, 0.0f, 0.0f};
    bromath::Vec3 v1{0.0f, 0.0f, 0.0f};
    float f0 = 0.0f;
    float f1 = 0.0f;
};

class SdfGraph {
public:
    SdfGraph() = default;

    // Primitives
    int sphere(float radius);
    int box(bromath::Vec3 halfExtents);
    int roundedBox(bromath::Vec3 halfExtents, float radius);
    int cylinder(float radius, float halfHeight);
    int capsule(bromath::Vec3 a, bromath::Vec3 b, float radius);
    int torus(float majorR, float minorR);
    int plane(bromath::Vec3 normal, float d);

    // CSG operations
    int opUnion(int a, int b);
    int opIntersection(int a, int b);
    int opSubtraction(int a, int b);
    int opSmoothUnion(int a, int b, float k);
    int opSmoothIntersection(int a, int b, float k);
    int opSmoothSubtraction(int a, int b, float k);

    // Transforms & Modifiers
    int translate(int child, bromath::Vec3 offset);
    int rotateY(int child, float angleRad);
    int scaleUniform(int child, float s);
    int noise3D(bromath::Vec3 freq, float amplitude);
    int displace(int child, int noiseNode);

    // Root node management
    int root() const;
    int rootNode() const { return root(); }
    void setRoot(int root);

    // Node inspection
    const std::vector<SdfNode>& nodes() const { return m_nodes; }
    const SdfNode& node(int id) const { return m_nodes.at(static_cast<size_t>(id)); }
    size_t size() const { return m_nodes.size(); }
    bool empty() const { return m_nodes.empty(); }
    void clear();

    // Deterministic hash computation
    uint64_t computeHash() const;

private:
    std::vector<SdfNode> m_nodes;
    int m_root = -1;
};

} // namespace bromesh
