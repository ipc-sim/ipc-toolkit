// Device-side utilities shared by the ipc::cuda implementation files.
// This header is CUDA-only and must be included from .cu files exclusively.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/logger.hpp>

#include <cuda_runtime.h>

/// @brief Log and throw (ipc::log_and_throw_error) if a CUDA runtime call
/// fails.
#define IPC_TOOLKIT_CUDA_CHECK(expr)                                           \
    do {                                                                       \
        const cudaError_t ipc_cuda_check_err = (expr);                         \
        if (ipc_cuda_check_err != cudaSuccess) {                               \
            ::ipc::log_and_throw_error(                                        \
                "CUDA error at {}:{}: {}", __FILE__, __LINE__,                 \
                cudaGetErrorString(ipc_cuda_check_err));                       \
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

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
