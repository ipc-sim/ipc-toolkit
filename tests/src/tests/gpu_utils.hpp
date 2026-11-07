#pragma once

// Shared scaffolding for the GPU parity tests
// (tests/src/tests/**/test_gpu_*.cu).
//
// Each of those files launches one-thread kernels over a flat double buffer and
// checks the results against the host implementations. The launch, transfer,
// and comparison are identical everywhere, so they live here; only the kernels
// and their fixtures belong in the individual files.
//
// This header defines __device__ code, so include it from a .cu, inside the
// file's `#ifdef IPC_TOOLKIT_WITH_CUDA` guard.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

#ifndef __CUDACC__
#error "tests/gpu_utils.hpp defines device code; include it from a .cu file."
#endif

#define REQUIRE_CUDA(expr) REQUIRE((expr) == cudaSuccess)

namespace ipc::tests {

/// @brief Uniform kernel signature: unpack inputs from `in`, write results
/// flattened (row-major) to `out`.
using GpuTestKernel = void (*)(const double* in, double* out);

/// @brief Copy a fixed-size Eigen object into `out` row-major; returns the
/// number of doubles written.
template <typename Derived>
__device__ int
write_out(const Eigen::MatrixBase<Derived>& m, double* out, int offset)
{
    for (int r = 0; r < m.rows(); ++r) {
        for (int c = 0; c < m.cols(); ++c) {
            out[offset++] = m(r, c);
        }
    }
    return offset;
}

/// @brief Append a fixed-size Eigen object to a host vector row-major.
template <typename Derived>
inline void write_expected(
    const Eigen::MatrixBase<Derived>& m, std::vector<double>& expected)
{
    for (int r = 0; r < m.rows(); ++r) {
        for (int c = 0; c < m.cols(); ++c) {
            expected.push_back(m(r, c));
        }
    }
}

/// @brief SKIP the test if no CUDA device is present (kernels have still been
/// compiled as device code, which is the primary verification).
inline void skip_if_no_cuda_device()
{
    int device_count = 0;
    const cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }
}

/// @brief Launch `kernel` with the given inputs; returns the device outputs.
inline std::vector<double> run_gpu_kernel(
    GpuTestKernel kernel, const std::vector<double>& in, const size_t n_out)
{
    double *d_in = nullptr, *d_out = nullptr;
    REQUIRE_CUDA(cudaMalloc(&d_in, in.size() * sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_out, n_out * sizeof(double)));
    REQUIRE_CUDA(cudaMemcpy(
        d_in, in.data(), in.size() * sizeof(double), cudaMemcpyHostToDevice));

    kernel<<<1, 1>>>(d_in, d_out);
    REQUIRE_CUDA(cudaGetLastError());
    REQUIRE_CUDA(cudaDeviceSynchronize());

    std::vector<double> out(n_out);
    REQUIRE_CUDA(cudaMemcpy(
        out.data(), d_out, n_out * sizeof(double), cudaMemcpyDeviceToHost));
    cudaFree(d_in);
    cudaFree(d_out);
    return out;
}

/// @brief Run `kernel` and check every output against the host reference.
///
/// An infinite reference is compared by sign rather than by margin, so a
/// barrier evaluated at d <= 0 can be checked alongside finite values.
inline void check_gpu_matches_host(
    GpuTestKernel kernel,
    const std::vector<double>& in,
    const std::vector<double>& expected)
{
    const std::vector<double> out = run_gpu_kernel(kernel, in, expected.size());
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

} // namespace ipc::tests
