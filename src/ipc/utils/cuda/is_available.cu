#include "is_available.hpp"

#include <cuda_runtime.h>

namespace ipc::cuda {

bool is_available()
{
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
        // Without a driver or a device the query fails (e.g.,
        // cudaErrorInsufficientDriver or cudaErrorNoDevice). We clear the
        // error so a later IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError()) does not
        // report it as its own.
        static_cast<void>(cudaGetLastError());
        return false;
    }
    return device_count > 0;
}

} // namespace ipc::cuda
