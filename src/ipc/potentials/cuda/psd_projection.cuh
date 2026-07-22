// Batched positive-semi-definite projection of small dense symmetric matrices
// on the GPU. Each block is fully eigendecomposed by a single CUDA thread with
// Eigen's SelfAdjointEigenSolver (QR iteration), then its eigenvalues are
// clamped (CLAMP) or made positive (ABS) and the block is reconstructed in
// place. CUDA-only; include from .cu files.
//
// This one-thread-per-matrix model (following libuipc's make_spd /
// muda::eigen::evd) is ~10x faster than cuSOLVER's batched Jacobi for a large
// batch of tiny 12x12 blocks, whose per-matrix API overhead dominates, and it
// uses the same eigensolver as the CPU path so results match to reduction-order
// noise.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/eigen_ext.hpp> // PSDProjectionMethod

namespace ipc::cuda {

/// @brief Projects a batch of symmetric matrices onto the PSD cone on device.
///
/// Every matrix is stored column-major as a fixed DIM×DIM block (DIM = 12),
/// contiguously in the batch buffer; blocks smaller than DIM must be
/// zero-padded (the padding contributes zero eigenvalues that project to zero).
/// Projection is done in place, one CUDA thread per block. The class is
/// stateless (no retained device scratch); it is kept as a small RAII object so
/// callers can hold one across repeated calls (e.g. Newton iterations) without
/// changing the call site.
class BatchedPSDProjection {
public:
    /// @brief Fixed block dimension of the eigensolver.
    static constexpr int DIM = 12;
    /// @brief Number of doubles per (column-major) block slot.
    static constexpr int BLOCK_SIZE = DIM * DIM;

    BatchedPSDProjection() = default;
    ~BatchedPSDProjection() = default;

    BatchedPSDProjection(const BatchedPSDProjection&) = delete;
    BatchedPSDProjection& operator=(const BatchedPSDProjection&) = delete;

    /// @brief Project @p num_blocks symmetric DIM×DIM blocks in place.
    /// @param[in,out] d_blocks Device pointer to num_blocks × BLOCK_SIZE doubles
    /// (column-major blocks). On return each block holds its PSD projection.
    /// @param num_blocks Number of blocks in the batch.
    /// @param method Projection method (CLAMP or ABS; NONE is a no-op).
    void project(
        double* d_blocks,
        const int num_blocks,
        const PSDProjectionMethod method);
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
