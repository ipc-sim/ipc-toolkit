// Device-side utilities shared by the ipc::cuda implementation files.
// This header is CUDA-only and must be included from .cu files exclusively.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

// atomicAdd(double*, double) requires compute capability 6.0+.
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
#error "ipc::cuda requires compute capability 6.0+ (atomicAdd on double)."
#endif

/// @brief Throw a std::runtime_error if a CUDA runtime call fails.
#define IPC_TOOLKIT_CUDA_CHECK(expr)                                           \
    do {                                                                       \
        const cudaError_t ipc_cuda_check_err = (expr);                         \
        if (ipc_cuda_check_err != cudaSuccess) {                               \
            throw std::runtime_error(                                          \
                std::string("CUDA error at " __FILE__ ":")                     \
                + std::to_string(__LINE__) + ": "                              \
                + cudaGetErrorString(ipc_cuda_check_err));                     \
        }                                                                      \
    } while (false)

namespace ipc::cuda {

/// @brief Number of threads per block used by the ipc::cuda kernels.
constexpr int KERNEL_BLOCK_SIZE = 256;

/// @brief Compute the launch grid size for @p n threads.
inline int kernel_grid_size(const size_t n)
{
    return static_cast<int>((n + KERNEL_BLOCK_SIZE - 1) / KERNEL_BLOCK_SIZE);
}

/// @brief Global DOF index of component @p d of vertex @p vertex_id.
/// Mirrors the index math of local_gradient_to_global_gradient()
/// (see src/ipc/utils/local_to_global.hpp) for dim=3.
__device__ inline index_t global_dof_index(
    const index_t vertex_id, const int d, const index_t n_total_vertices)
{
    if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
        return 3 * vertex_id + d;
    } else {
        return n_total_vertices * d + vertex_id;
    }
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
