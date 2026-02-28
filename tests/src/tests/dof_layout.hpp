/// @file
/// Utility helpers that make finite-difference tests work regardless of
/// whether ipc::VERTEX_DERIVATIVE_LAYOUT is RowMajor or ColMajor.
///
/// fd::flatten / fd::unflatten always use **RowMajor** ordering:
///   flatten   :  MatrixXd(n_verts, dim)  →  [x0 y0 z0 x1 y1 z1 …]
///   unflatten :  [x0 y0 z0 x1 y1 z1 …]  →  MatrixXd(n_verts, dim)
///
/// When VERTEX_DERIVATIVE_LAYOUT == ColMajor the library returns gradient /
/// hessian values laid out as
///   [x0 x1 …  y0 y1 …  z0 z1 …]   (columns-of-vertices first)
///
/// The helpers below bridge the two worlds so that every test can be written as
///
///   fd::finite_gradient(dof_flatten(V), [&](auto& x){ … dof_unflatten(x, V) …
///   }, fgrad);
///   fgrad = dof_reorder_grad(fgrad, n_verts, dim);   // RowMajor → layout
///
/// and the comparison with the library result is valid.

#pragma once

#include <ipc/config.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <cassert>

namespace ipc::tests {

// ─── Flatten / Unflatten in VERTEX_DERIVATIVE_LAYOUT order ──────────────────

/// Flatten a vertex matrix (n_verts × dim) to a DOF vector in the order
/// dictated by VERTEX_DERIVATIVE_LAYOUT.
///
///   RowMajor → [x0 y0 z0  x1 y1 z1  …]   (same as fd::flatten)
///   ColMajor → [x0 x1 …   y0 y1 …   z0 z1 …]
inline auto flatten(const Eigen::MatrixXd& V)
{
    return V.reshaped<VERTEX_DERIVATIVE_LAYOUT>();
}

/// Unflatten a DOF vector back to a vertex matrix (n_verts × dim).
/// The inverse of dof_flatten.
inline auto unflatten(const Eigen::VectorXd& x, int dim)
{
    return x.reshaped<VERTEX_DERIVATIVE_LAYOUT>(x.size() / dim, dim);
}

// ─── Reorder FD results from RowMajor to VERTEX_DERIVATIVE_LAYOUT ──────────
//
// These are used when the FD callback is forced to stay in RowMajor order
// (e.g. because it calls fd::flatten/unflatten internally) but the library
// output is in VERTEX_DERIVATIVE_LAYOUT.

/// Reorder a gradient vector computed with RowMajor flattening so that it
/// matches VERTEX_DERIVATIVE_LAYOUT.
///
/// When the layout is already RowMajor this is a no-op (compiler will
/// optimize the copy away or it can be elided).
inline auto
reorder_gradient(const Eigen::VectorXd& g_rowmajor, int n_verts, int dim)
{
    if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
        return g_rowmajor; // Already in the correct order, no copy needed.
    } else {
        return g_rowmajor.reshaped<Eigen::RowMajor>(n_verts, dim)
            .reshaped<VERTEX_DERIVATIVE_LAYOUT>()
            .eval();
    }
}

/// Reorder a sparse Hessian from RowMajor DOF ordering to
/// VERTEX_DERIVATIVE_LAYOUT DOF ordering.
inline auto reorder_hessian(
    const Eigen::SparseMatrix<double>& H_rowmajor, int n_verts, int dim)
{
    if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
        return H_rowmajor; // Already in the correct order, no copy needed.
    } else {
        const int n = n_verts * dim;
        assert(H_rowmajor.rows() == n && H_rowmajor.cols() == n);
        std::vector<Eigen::Triplet<double>> trips;
        trips.reserve(H_rowmajor.nonZeros());
        for (int k = 0; k < H_rowmajor.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(H_rowmajor, k);
                 it; ++it) {
                // row_major → col_major
                trips.emplace_back(
                    (it.row() % dim) * n_verts + it.row() / dim,
                    (it.col() % dim) * n_verts + it.col() / dim, it.value());
            }
        }
        Eigen::SparseMatrix<double> H(n, n);
        H.setFromTriplets(trips.begin(), trips.end());
        return H;
    }
}

} // namespace ipc::tests