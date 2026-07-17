// Internal helper: CPU-side assembly of local hessian blocks computed on the
// GPU. Kept in a separate host-compiled translation unit so the CUDA file
// does not need TBB or the sparse-assembly headers.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <array>
#include <vector>

namespace ipc::cuda::internal {

/// @brief Size of one local hessian block slot (12 × 12 doubles).
constexpr int HESSIAN_BLOCK_SLOT_SIZE = 144;

/// @brief Append a slab of GPU-computed local hessian blocks to triplets.
///
/// Each block slot holds an ndof × ndof column-major matrix (ndof = 3 ×
/// number of stencil vertices, determined by counting non-negative vertex
/// ids). Blocks are PSD-projected (per @p project_hessian_to_psd) before
/// being scattered, matching the CPU ordering (weights are already
/// multiplied in on the GPU).
///
/// @param[in] blocks Host buffer of @p slab_count block slots.
/// @param[in] slab_begin Flat index of the first collision in the slab.
/// @param[in] slab_count Number of collisions in the slab.
/// @param[in] vertex_ids Vertex ids of every collision (flattened order).
/// @param[in] project_hessian_to_psd Method of projecting the local hessians
/// to the positive semi-definite cone.
/// @param[in] n_total_vertices Total number of vertices in the mesh.
/// @param[out] triplets Vector to which the block triplets are appended.
void append_local_hessian_blocks(
    const double* blocks,
    const size_t slab_begin,
    const size_t slab_count,
    const std::vector<std::array<index_t, 4>>& vertex_ids,
    const PSDProjectionMethod project_hessian_to_psd,
    const index_t n_total_vertices,
    std::vector<Eigen::Triplet<double>>& triplets);

} // namespace ipc::cuda::internal

#endif // IPC_TOOLKIT_WITH_CUDA
