#include "bromesh/optimization/spatial_hash.h"

#include <cmath>
#include <limits>

namespace bromesh {

SpatialHash3D::SpatialHash3D(float cellSize) {
    reset(cellSize);
}

void SpatialHash3D::reset(float cellSize) {
    cellSize_ = (cellSize > 0.0f) ? cellSize : 1.0f;
    invCell_ = 1.0f / cellSize_;
    cells_.clear();
    count_ = 0;
}

void SpatialHash3D::clear() {
    cells_.clear();
    count_ = 0;
}

SpatialHash3D::CellKey SpatialHash3D::makeKey(int ix, int iy, int iz) {
    // Pack 3 x 21-bit signed ints into one 64-bit key (range ~+-1M cells).
    constexpr int64_t mask = (1LL << 21) - 1;
    int64_t a = static_cast<int64_t>(ix) & mask;
    int64_t b = static_cast<int64_t>(iy) & mask;
    int64_t c = static_cast<int64_t>(iz) & mask;
    return a | (b << 21) | (c << 42);
}

void SpatialHash3D::cellOf(Vec3 p, int& ix, int& iy, int& iz) const {
    ix = static_cast<int>(std::floor(p.x * invCell_));
    iy = static_cast<int>(std::floor(p.y * invCell_));
    iz = static_cast<int>(std::floor(p.z * invCell_));
}

void SpatialHash3D::insert(Vec3 p, int32_t id) {
    int ix, iy, iz;
    cellOf(p, ix, iy, iz);
    cells_[makeKey(ix, iy, iz)].push_back({p, id});
    ++count_;
}

void SpatialHash3D::remove(int32_t id) {
    for (auto it = cells_.begin(); it != cells_.end(); ) {
        auto& vec = it->second;
        for (size_t i = vec.size(); i > 0; --i) {
            if (vec[i - 1].id == id) {
                vec[i - 1] = vec.back();
                vec.pop_back();
                --count_;
            }
        }
        if (vec.empty()) it = cells_.erase(it);
        else ++it;
    }
}

void SpatialHash3D::radiusQuery(Vec3 center, float radius, std::vector<int32_t>& out) const {
    if (radius <= 0.0f) return;
    const float r2 = radius * radius;
    const int extent = static_cast<int>(std::ceil(radius * invCell_));
    int cx, cy, cz;
    cellOf(center, cx, cy, cz);
    for (int dz = -extent; dz <= extent; ++dz) {
        for (int dy = -extent; dy <= extent; ++dy) {
            for (int dx = -extent; dx <= extent; ++dx) {
                auto it = cells_.find(makeKey(cx + dx, cy + dy, cz + dz));
                if (it == cells_.end()) continue;
                for (const Entry& e : it->second) {
                    if (vdist2(e.p, center) <= r2) out.push_back(e.id);
                }
            }
        }
    }
}

int32_t SpatialHash3D::nearest(Vec3 center, float maxRadius) const {
    if (maxRadius <= 0.0f) return -1;
    int extent = static_cast<int>(std::ceil(maxRadius * invCell_));
    int cx, cy, cz;
    cellOf(center, cx, cy, cz);
    int32_t bestId = -1;
    float bestD2 = maxRadius * maxRadius;
    for (int dz = -extent; dz <= extent; ++dz) {
        for (int dy = -extent; dy <= extent; ++dy) {
            for (int dx = -extent; dx <= extent; ++dx) {
                auto it = cells_.find(makeKey(cx + dx, cy + dy, cz + dz));
                if (it == cells_.end()) continue;
                for (const Entry& e : it->second) {
                    float d2 = vdist2(e.p, center);
                    if (d2 < bestD2) {
                        bestD2 = d2;
                        bestId = e.id;
                    }
                }
            }
        }
    }
    return bestId;
}

} // namespace bromesh
