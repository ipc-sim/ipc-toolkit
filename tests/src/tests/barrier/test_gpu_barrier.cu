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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace {

#define REQUIRE_CUDA(expr) REQUIRE((expr) == cudaSuccess)

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
    int device_count = 0;
    const cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        SKIP("No CUDA device available; kernel compiled but not executed.");
    }

    const double dhat = 1e-2;
    const std::vector<double> in = { dhat, 0.5 * dhat, 2 * dhat, -1.0 };

    std::vector<double> expected;
    for (int i = 0; i < 3; ++i) {
        const double d = in[1 + i];
        expected.push_back(ipc::barrier(d, dhat));
        expected.push_back(ipc::barrier_first_derivative(d, dhat));
        expected.push_back(ipc::barrier_second_derivative(d, dhat));
    }

    double *d_in = nullptr, *d_out = nullptr;
    REQUIRE_CUDA(cudaMalloc(&d_in, in.size() * sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_out, expected.size() * sizeof(double)));
    REQUIRE_CUDA(cudaMemcpy(
        d_in, in.data(), in.size() * sizeof(double), cudaMemcpyHostToDevice));

    barrier_kernel<<<1, 1>>>(d_in, d_out);
    REQUIRE_CUDA(cudaGetLastError());
    REQUIRE_CUDA(cudaDeviceSynchronize());

    std::vector<double> out(expected.size());
    REQUIRE_CUDA(cudaMemcpy(
        out.data(), d_out, out.size() * sizeof(double),
        cudaMemcpyDeviceToHost));
    cudaFree(d_in);
    cudaFree(d_out);

    for (size_t i = 0; i < expected.size(); ++i) {
        CAPTURE(i);
        if (std::isinf(expected[i])) {
            CHECK(std::isinf(out[i]));
            CHECK((out[i] > 0) == (expected[i] > 0));
        } else {
            CHECK(out[i] == Catch::Approx(expected[i]).margin(1e-12));
        }
    }
}

#endif // IPC_TOOLKIT_WITH_CUDA
