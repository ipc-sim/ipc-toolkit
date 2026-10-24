// Exercises the __host__ __device__ barrier free functions from a CUDA
// kernel.
//
// Compilation proves the framework works: nvcc must instantiate the barrier
// free-function trio as device code. On a machine with a CUDA device the
// kernel also runs and its outputs are checked against the host
// implementations — host/device parity is the contract. Without a device
// (e.g. a GPU-less container) the runtime checks are skipped but the device
// code has still been compiled.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/barrier/barrier.hpp>

#include <tests/gpu_utils.hpp>

#include <cuda_runtime.h>

#include <vector>

using namespace ipc::tests;

namespace {

__global__ void barrier_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const double dhat = in[0];
    // in[1..3]: interior, inactive (d >= dhat), and invalid (d <= 0) samples
    for (int i = 0; i < 3; ++i) {
        const double d = in[1 + i];
        out[3 * i + 0] = ipc::barrier(d, dhat);
        out[3 * i + 1] = ipc::barrier_first_derivative(d, dhat);
        out[3 * i + 2] = ipc::barrier_second_derivative(d, dhat);
    }
}

} // namespace

TEST_CASE("GPU barrier", "[barrier][gpu]")
{
    skip_if_no_cuda_device();

    const double dhat = 1e-2;
    const std::vector<double> in = { dhat, 0.5 * dhat, 2 * dhat, -1.0 };

    std::vector<double> expected;
    for (int i = 0; i < 3; ++i) {
        const double d = in[1 + i];
        expected.push_back(ipc::barrier(d, dhat));
        expected.push_back(ipc::barrier_first_derivative(d, dhat));
        expected.push_back(ipc::barrier_second_derivative(d, dhat));
    }

    check_gpu_matches_host(barrier_kernel, in, expected);
}

#endif // IPC_TOOLKIT_WITH_CUDA
