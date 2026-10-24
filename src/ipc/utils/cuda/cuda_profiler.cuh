#pragma once

#include <ipc/utils/profiler.hpp>

#ifdef IPC_TOOLKIT_WITH_PROFILER

#include <ipc/utils/cuda/device_utils.cuh>

#include <cuda_runtime.h>

namespace ipc::cuda {

/// @brief A ProfilePoint timer that measures device time with CUDA events.
///
/// Kernel launches are asynchronous, so a host clock around one measures the
/// launch and not the work. Recording events on the stream instead measures
/// when the device ran the enclosed stage.
///
/// Reading the elapsed time waits on the stop event, so a profiling build
/// serializes stages that would otherwise overlap and its totals are larger
/// than the same build without the profiler. Stage *shares* are what this is
/// for; take absolute build and detect times from an uninstrumented build.
class CudaEventTimer {
public:
    CudaEventTimer()
    {
        IPC_TOOLKIT_CUDA_CHECK(cudaEventCreate(&m_start));
        IPC_TOOLKIT_CUDA_CHECK(cudaEventCreate(&m_stop));
    }

    ~CudaEventTimer()
    {
        cudaEventDestroy(m_start);
        cudaEventDestroy(m_stop);
    }

    CudaEventTimer(const CudaEventTimer&) = delete;
    CudaEventTimer& operator=(const CudaEventTimer&) = delete;

    void start() { IPC_TOOLKIT_CUDA_CHECK(cudaEventRecord(m_start)); }

    void stop() { IPC_TOOLKIT_CUDA_CHECK(cudaEventRecord(m_stop)); }

    // NOLINTNEXTLINE(readability-identifier-naming)
    double getElapsedTimeInMilliSec() const
    {
        IPC_TOOLKIT_CUDA_CHECK(cudaEventSynchronize(m_stop));
        float ms = 0;
        IPC_TOOLKIT_CUDA_CHECK(cudaEventElapsedTime(&ms, m_start, m_stop));
        return static_cast<double>(ms);
    }

private:
    cudaEvent_t m_start {}, m_stop {};
};

} // namespace ipc::cuda

/// @brief Profile a device-side stage, timed with CUDA events.
#define IPC_TOOLKIT_PROFILE_BLOCK_CUDA(...)                                    \
    ipc::ProfilePoint<ipc::cuda::CudaEventTimer>                               \
        IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(                                      \
            __ipc_cuda_profile_point_, __COUNTER__)(__VA_ARGS__);              \
    ZoneScopedN(__VA_ARGS__)

#else

#define IPC_TOOLKIT_PROFILE_BLOCK_CUDA(...) ZoneScopedN(__VA_ARGS__)

#endif
