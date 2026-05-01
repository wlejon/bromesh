#pragma once

#include "bromesh/procedural/vec_math.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace bromesh {

/// Uniform-grid 3D spatial hash for point queries.
///
/// Points are inserted with an arbitrary `int32_t` payload id; the same id
/// can be used to remove the point later. Queries return matching ids in
/// no particular order. Single-threaded; the project forbids mutexes.
class SpatialHash3D {
public:
    explicit SpatialHash3D(float cellSize = 1.0f);

    /// Reset cell size and clear all entries. Cell size must be > 0.
    void reset(float cellSize);

    /// Remove all entries.
    void clear();

    /// Insert a point with the given payload id. Duplicate ids are allowed
    /// but `remove(id)` will remove all of them.
    void insert(Vec3 p, int32_t id);

    /// Remove every entry whose id matches. O(N) over the touched cells.
    void remove(int32_t id);

    /// Append all ids whose points lie within `radius` of `center`.
    /// Existing contents of `out` are NOT cleared.
    void radiusQuery(Vec3 center, float radius, std::vector<int32_t>& out) const;

    /// Return the id of the nearest point within `maxRadius`, or -1 if none.
    int32_t nearest(Vec3 center, float maxRadius) const;

    /// Number of stored points.
    size_t size() const { return count_; }

    float cellSize() const { return cellSize_; }

private:
    struct Entry {
        Vec3 p;
        int32_t id;
    };
    using CellKey = int64_t;

    static CellKey makeKey(int ix, int iy, int iz);
    void cellOf(Vec3 p, int& ix, int& iy, int& iz) const;

    float cellSize_ = 1.0f;
    float invCell_ = 1.0f;
    size_t count_ = 0;
    std::unordered_map<CellKey, std::vector<Entry>> cells_;
};

} // namespace bromesh
