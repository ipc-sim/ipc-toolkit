// Batched positive-semi-definite projection of small dense symmetric matrices
// on the GPU, backed by cuSOLVER's batched Jacobi eigensolver
// (cusolverDnDsyevjBatched). CUDA-only; include from .cu files.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/eigen_ext.hpp> // PSDProjectionMethod

#include <cusolverDn.h>

namespace ipc::cuda {

/// @brief Projects a batch of symmetric matrices onto the PSD cone on device.
///
/// Every matrix is stored column-major as a fixed DIM×DIM block (DIM = 12),
/// contiguously in the batch buffer; blocks smaller than DIM must be
/// zero-padded (the padding contributes zero eigenvalues that project to zero).
/// The cuSOLVER handle and eigenvalue/workspace buffers are retained so the
/// setup cost is amortized across repeated calls (e.g. Newton iterations).
class BatchedPSDProjection {
public:
    /// @brief Fixed block dimension of the batched eigensolver.
    static constexpr int DIM = 12;
    /// @brief Number of doubles per (column-major) block slot.
    static constexpr int BLOCK_SIZE = DIM * DIM;

    BatchedPSDProjection();
    ~BatchedPSDProjection();

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

private:
    void ensure_capacity(const int num_blocks);

    cusolverDnHandle_t m_handle = nullptr;
    syevjInfo_t m_syevj_params = nullptr;

    // Retained device scratch (grown on demand).
    double* m_d_eigenvalues = nullptr;  // num_blocks × DIM
    double* m_d_eigenvectors = nullptr; // num_blocks × BLOCK_SIZE (syevj input)
    double* m_d_work = nullptr;         // cuSOLVER workspace
    int* m_d_info = nullptr;            // num_blocks
    int m_capacity_blocks = 0;          // blocks the scratch is sized for
    int m_lwork = 0;                    // doubles in m_d_work
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
