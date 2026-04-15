#pragma once

#include "bromesh/mesh_data.h"

#include <vector>

namespace bromesh {

/// Sparse matrix in Compressed Sparse Row format. Indexing convention:
///   row i occupies entries [rowStart[i], rowStart[i+1]).
///   colIndex[k] is the column of the k-th entry; values[k] is the value.
/// Rows are not required to be sorted by column, but assembleCotangentLaplacian
/// produces sorted rows.
struct SparseCsr {
    int rows = 0;
    int cols = 0;
    std::vector<int> rowStart;   // length rows+1
    std::vector<int> colIndex;   // length nnz
    std::vector<double> values;  // length nnz

    int nnz() const { return (int)values.size(); }
};

/// Build the cotangent Laplacian L and lumped Voronoi mass M for a triangle
/// mesh. Convention: L is the "integrated" Laplacian — for a scalar field f,
/// (L f)[i] = 0.5 * sum_j (cot α + cot β) (f[j] - f[i]) over edges (i,j).
///
/// Cotangent weights are clamped at zero on obtuse triangles (intrinsic
/// Laplacian trick): preserves the M-matrix property and keeps downstream
/// solves stable on meshes with occasional bad triangles.
///
/// Mass is lumped (diagonal) with "mixed Voronoi" areas: Voronoi for acute
/// triangles, falls back to triangle-area/3 for obtuse ones. Guaranteed
/// positive.
///
/// Mesh must have at least one triangle; output dimensions == vertexCount.
void assembleCotangentLaplacian(const MeshData& mesh,
                                SparseCsr& L,
                                std::vector<double>& massDiag);

/// Conjugate-gradient solve of A x = b, where A is symmetric positive
/// definite in CSR form. Jacobi preconditioner (inverse diagonal). Warm-
/// start with x already populated; pass zero-initialized x for a cold start.
///
/// Returns number of iterations taken. Bails after maxIter or once
/// ||r|| / ||b|| < tol.
int sparseCg(const SparseCsr& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double tol = 1e-8,
             int maxIter = 2000);

/// Sparse matrix-vector product y = A x. Vectors must be sized A.cols / A.rows
/// respectively. Does not zero y first.
void sparseSpmv(const SparseCsr& A,
                const std::vector<double>& x,
                std::vector<double>& y);

} // namespace bromesh
