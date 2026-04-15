#include "bromesh/rigging/weight_smooth.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <unordered_set>
#include <vector>

namespace bromesh {

namespace {

constexpr int K = 4; // glTF-standard influences per vertex

} // namespace

void postProcessWeights(const MeshData& mesh,
                        SkinData& skin,
                        const WeightPostProcessOptions& opts) {
    const size_t nVerts = mesh.vertexCount();
    const size_t B = skin.boneCount;
    if (nVerts == 0 || B == 0) return;
    if (skin.boneWeights.size() != nVerts * K) return;
    if (skin.boneIndices.size() != nVerts * K) return;

    // Mesh-edge adjacency. Dedup via per-vertex set to keep the Laplacian
    // weights uniform (each neighbor counted once even if it shares multiple
    // triangles with v).
    std::vector<std::unordered_set<uint32_t>> adjSet(nVerts);
    for (size_t t = 0; t + 2 < mesh.indices.size(); t += 3) {
        uint32_t a = mesh.indices[t];
        uint32_t b = mesh.indices[t + 1];
        uint32_t c = mesh.indices[t + 2];
        if (a >= nVerts || b >= nVerts || c >= nVerts) continue;
        adjSet[a].insert(b); adjSet[a].insert(c);
        adjSet[b].insert(a); adjSet[b].insert(c);
        adjSet[c].insert(a); adjSet[c].insert(b);
    }
    std::vector<std::vector<uint32_t>> adj(nVerts);
    for (size_t i = 0; i < nVerts; ++i) {
        adj[i].assign(adjSet[i].begin(), adjSet[i].end());
        std::sort(adj[i].begin(), adj[i].end()); // deterministic iteration
    }

    // Expand top-K to dense per-bone. Duplicate indices are summed (harmless).
    std::vector<float> dense(nVerts * B, 0.0f);
    for (size_t v = 0; v < nVerts; ++v) {
        for (int k = 0; k < K; ++k) {
            float w = skin.boneWeights[v * K + k];
            uint32_t bi = skin.boneIndices[v * K + k];
            if (w > 0.0f && bi < B) dense[v * B + bi] += w;
        }
    }

    // Laplacian smoothing. alpha = blend toward neighborhood mean.
    if (opts.iterations > 0 && opts.alpha > 0.0f) {
        std::vector<float> next(dense.size());
        const float a = std::clamp(opts.alpha, 0.0f, 1.0f);
        for (int it = 0; it < opts.iterations; ++it) {
            for (size_t v = 0; v < nVerts; ++v) {
                if (adj[v].empty()) {
                    std::memcpy(&next[v * B], &dense[v * B], B * sizeof(float));
                    continue;
                }
                float inv = 1.0f / (float)adj[v].size();
                for (size_t bi = 0; bi < B; ++bi) {
                    float s = 0.0f;
                    for (uint32_t n : adj[v]) s += dense[(size_t)n * B + bi];
                    float mean = s * inv;
                    next[v * B + bi] = (1.0f - a) * dense[v * B + bi] + a * mean;
                }
            }
            dense.swap(next);
        }
    }

    // Re-select top-K, apply min-weight prune, re-normalize. Bones absent
    // from the smoothed distribution naturally drop out — this is the
    // outlier-rejection step: a grid-induced stray influence with no
    // neighborhood support shrinks below minWeight and gets pruned.
    std::vector<float>   topW(K, 0.0f);
    std::vector<int32_t> topB(K, -1);
    for (size_t v = 0; v < nVerts; ++v) {
        std::fill(topW.begin(), topW.end(), 0.0f);
        std::fill(topB.begin(), topB.end(), -1);
        for (size_t bi = 0; bi < B; ++bi) {
            float w = dense[v * B + bi];
            if (w <= topW[K - 1]) continue;
            topW[K - 1] = w;
            topB[K - 1] = (int32_t)bi;
            for (int i = K - 1; i > 0 && topW[i] > topW[i - 1]; --i) {
                std::swap(topW[i], topW[i - 1]);
                std::swap(topB[i], topB[i - 1]);
            }
        }
        float sum = 0.0f;
        for (int k = 0; k < K; ++k) {
            if (topW[k] < opts.minWeight) { topW[k] = 0.0f; topB[k] = -1; }
            sum += topW[k];
        }
        if (sum > 0.0f) {
            float inv = 1.0f / sum;
            for (int k = 0; k < K; ++k) {
                skin.boneWeights[v * K + k] = topW[k] * inv;
                skin.boneIndices[v * K + k] =
                    (topB[k] < 0) ? 0u : (uint32_t)topB[k];
            }
        }
        // If sum is zero (vertex had no valid weights even pre-smooth) leave
        // the existing per-vertex entry as-is so we don't silently orphan it.
    }
}

} // namespace bromesh
