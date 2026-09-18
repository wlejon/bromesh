#include "bromesh/isosurface/jit/sdf_node.h"

#include <cstdint>
#include <cstring>

namespace bromesh {

int SdfGraph::sphere(float radius) {
    SdfNode node;
    node.op = SdfOp::Sphere;
    node.f0 = radius;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::box(bromath::Vec3 halfExtents) {
    SdfNode node;
    node.op = SdfOp::Box;
    node.v0 = halfExtents;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::roundedBox(bromath::Vec3 halfExtents, float radius) {
    SdfNode node;
    node.op = SdfOp::RoundedBox;
    node.v0 = halfExtents;
    node.f0 = radius;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::cylinder(float radius, float halfHeight) {
    SdfNode node;
    node.op = SdfOp::Cylinder;
    node.f0 = radius;
    node.f1 = halfHeight;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::capsule(bromath::Vec3 a, bromath::Vec3 b, float radius) {
    SdfNode node;
    node.op = SdfOp::Capsule;
    node.v0 = a;
    node.v1 = b;
    node.f0 = radius;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::torus(float majorR, float minorR) {
    SdfNode node;
    node.op = SdfOp::Torus;
    node.f0 = majorR;
    node.f1 = minorR;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::plane(bromath::Vec3 normal, float d) {
    SdfNode node;
    node.op = SdfOp::Plane;
    node.v0 = normal;
    node.f0 = d;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::opUnion(int a, int b) {
    SdfNode node;
    node.op = SdfOp::Union;
    node.left = a;
    node.right = b;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::opIntersection(int a, int b) {
    SdfNode node;
    node.op = SdfOp::Intersection;
    node.left = a;
    node.right = b;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::opSubtraction(int a, int b) {
    SdfNode node;
    node.op = SdfOp::Subtraction;
    node.left = a;
    node.right = b;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::opSmoothUnion(int a, int b, float k) {
    SdfNode node;
    node.op = SdfOp::SmoothUnion;
    node.left = a;
    node.right = b;
    node.f0 = k;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::opSmoothIntersection(int a, int b, float k) {
    SdfNode node;
    node.op = SdfOp::SmoothIntersection;
    node.left = a;
    node.right = b;
    node.f0 = k;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::opSmoothSubtraction(int a, int b, float k) {
    SdfNode node;
    node.op = SdfOp::SmoothSubtraction;
    node.left = a;
    node.right = b;
    node.f0 = k;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::translate(int child, bromath::Vec3 offset) {
    SdfNode node;
    node.op = SdfOp::Translate;
    node.left = child;
    node.v0 = offset;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::rotateY(int child, float angleRad) {
    SdfNode node;
    node.op = SdfOp::RotateY;
    node.left = child;
    node.f0 = angleRad;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::scaleUniform(int child, float s) {
    SdfNode node;
    node.op = SdfOp::ScaleUniform;
    node.left = child;
    node.f0 = s;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::noise3D(bromath::Vec3 freq, float amplitude) {
    SdfNode node;
    node.op = SdfOp::Noise3D;
    node.v0 = freq;
    node.f0 = amplitude;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::displace(int child, int noiseNode) {
    SdfNode node;
    node.op = SdfOp::Displace;
    node.left = child;
    node.right = noiseNode;
    m_nodes.push_back(node);
    m_root = static_cast<int>(m_nodes.size() - 1);
    return m_root;
}

int SdfGraph::root() const {
    if (m_root >= 0 && m_root < static_cast<int>(m_nodes.size())) {
        return m_root;
    }
    return m_nodes.empty() ? -1 : static_cast<int>(m_nodes.size() - 1);
}

void SdfGraph::setRoot(int root) {
    m_root = root;
}

void SdfGraph::clear() {
    m_nodes.clear();
    m_root = -1;
}

uint64_t SdfGraph::computeHash() const {
    uint64_t h = 14695981039346656037ULL;
    auto mix = [&](const void* data, size_t len) {
        const uint8_t* p = static_cast<const uint8_t*>(data);
        for (size_t i = 0; i < len; ++i) {
            h ^= static_cast<uint64_t>(p[i]);
            h *= 1099511628211ULL;
        }
    };

    int r = root();
    mix(&r, sizeof(r));
    uint64_t count = static_cast<uint64_t>(m_nodes.size());
    mix(&count, sizeof(count));

    for (const auto& n : m_nodes) {
        uint8_t opVal = static_cast<uint8_t>(n.op);
        mix(&opVal, sizeof(opVal));
        mix(&n.left, sizeof(n.left));
        mix(&n.right, sizeof(n.right));
        mix(&n.v0.x, sizeof(float));
        mix(&n.v0.y, sizeof(float));
        mix(&n.v0.z, sizeof(float));
        mix(&n.v1.x, sizeof(float));
        mix(&n.v1.y, sizeof(float));
        mix(&n.v1.z, sizeof(float));
        mix(&n.f0, sizeof(float));
        mix(&n.f1, sizeof(float));
    }

    return h;
}

} // namespace bromesh
