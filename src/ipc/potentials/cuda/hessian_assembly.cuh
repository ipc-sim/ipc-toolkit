// GPU assembly of a symmetric sparse hessian from per-collision dense blocks,
// using Thrust (sort + segmented reduce) to merge duplicate entries directly
// into CSC arrays -- no host-side setFromTriplets. CUDA-only.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/config.hpp> // index_t

#include <Eigen/Sparse>

namespace ipc::cuda {

/// @brief Assemble a symmetric sparse hessian from per-collision blocks.
///
/// Each collision contributes an ndof×ndof block (ndof = 3 × #stencil vertices,
/// determined per collision by counting its non-negative vertex ids), stored
/// column-major in the top-left of a zero-padded 12×12 slot at
/// @p d_blocks + i·144. Entry (3a+k, 3b+l) of collision i scatters to global
/// (dof(vertex_a, k), dof(vertex_b, l)) per VERTEX_DERIVATIVE_LAYOUT; duplicate
/// (row, col) contributions are summed.
///
/// @param d_blocks Device buffer of n × 144 doubles (padded 12×12 blocks).
/// @param d_vertex_id_0..3 Device arrays (n each) of stencil vertex ids (-1 for unused slots).
/// @param n Number of collisions.
/// @param num_vertices Total number of mesh vertices (n_dofs = 3·num_vertices).
/// @return The assembled sparse hessian (n_dofs × n_dofs, column-major).
Eigen::SparseMatrix<double> assemble_sparse_hessian(
    const double* d_blocks,
    const index_t* d_vertex_id_0,
    const index_t* d_vertex_id_1,
    const index_t* d_vertex_id_2,
    const index_t* d_vertex_id_3,
    const size_t n,
    const index_t num_vertices);

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
