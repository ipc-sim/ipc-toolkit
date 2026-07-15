// Exercises the __host__ __device__ point-point distance framework from an
// actual CUDA kernel.
//
// The value of this test is twofold:
//   1. Compilation proves the framework works: nvcc must be able to instantiate
//      ipc::point_point_distance (and its gradient/hessian) as *device* code,
//      including the Eigen operations they rely on. If the IPC_TOOLKIT_INLINE
//      __host__ __device__ machinery is wrong, this file fails to compile.
//   2. On a machine with a CUDA device, it also runs the kernel and checks the
//      results against the known-good CPU values. Where no device is present
//      (e.g. compiling in a GPU-less container), the run is skipped but the
//      kernel has still been compiled as device code.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/distance/point_point.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cuda_runtime.h>

#include <array>
#include <vector>

namespace {

/// @brief Evaluate the point-point distance, gradient, and hessian on device.
///
/// A single thread constructs the Eigen inputs from raw arrays, calls the
/// toolkit functions, and writes the results back to raw output arrays.
__global__ void point_point_kernel(
    const double* __restrict__ p0_data,
    const double* __restrict__ p1_data,
    const int dim,
    double* __restrict__ dist_out,
    double* __restrict__ grad_out,
    double* __restrict__ hess_out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }

    ipc::VectorMax3d p0(dim), p1(dim);
    for (int i = 0; i < dim; ++i) {
        p0(i) = p0_data[i];
        p1(i) = p1_data[i];
    }

    *dist_out = ipc::point_point_distance(p0, p1);

    const ipc::VectorMax6d grad = ipc::point_point_distance_gradient(p0, p1);
    std::memcpy(grad_out, grad.data(), grad.size() * sizeof(double));

    const ipc::MatrixMax6d hess = ipc::point_point_distance_hessian(p0, p1);
    std::memcpy(hess_out, hess.data(), hess.size() * sizeof(double));
}

/// @brief RAII helper so we bail out cleanly if any CUDA call fails.
#define REQUIRE_CUDA(expr) REQUIRE((expr) == cudaSuccess)

} // namespace

TEST_CASE("GPU point-point distance", "[distance][point_point][gpu]")
{
    // If no device is present the kernel above has still been compiled as
    // device code (the primary thing we are verifying); we just cannot run it.
    int device_count = 0;
    const cudaError_t count_err = cudaGetDeviceCount(&device_count);
    if (count_err != cudaSuccess || device_count == 0) {
        SKIP("No CUDA device available; kernel compiled but not executed.");
    }

    constexpr int dim = 3;
    constexpr int grad_size = 2 * dim;
    constexpr int hess_size = grad_size * grad_size;

    // p0 = (0, 0, 0), p1 = (1, 2, 2)  =>  squared distance = 1 + 4 + 4 = 9
    const std::vector<double> p0 = { 0.0, 0.0, 0.0 };
    const std::vector<double> p1 = { 1.0, 2.0, 2.0 };

    double *d_p0 = nullptr, *d_p1 = nullptr;
    double *d_dist = nullptr, *d_grad = nullptr, *d_hess = nullptr;
    REQUIRE_CUDA(cudaMalloc(&d_p0, dim * sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_p1, dim * sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_dist, sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_grad, grad_size * sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_hess, hess_size * sizeof(double)));

    REQUIRE_CUDA(cudaMemcpy(
        d_p0, p0.data(), dim * sizeof(double), cudaMemcpyHostToDevice));
    REQUIRE_CUDA(cudaMemcpy(
        d_p1, p1.data(), dim * sizeof(double), cudaMemcpyHostToDevice));

    point_point_kernel<<<1, 1>>>(d_p0, d_p1, dim, d_dist, d_grad, d_hess);
    REQUIRE_CUDA(cudaGetLastError());
    REQUIRE_CUDA(cudaDeviceSynchronize());

    double dist = 0.0;
    std::vector<double> grad(grad_size), hess(hess_size);
    REQUIRE_CUDA(
        cudaMemcpy(&dist, d_dist, sizeof(double), cudaMemcpyDeviceToHost));
    REQUIRE_CUDA(cudaMemcpy(
        grad.data(), d_grad, grad_size * sizeof(double),
        cudaMemcpyDeviceToHost));
    REQUIRE_CUDA(cudaMemcpy(
        hess.data(), d_hess, hess_size * sizeof(double),
        cudaMemcpyDeviceToHost));

    cudaFree(d_p0);
    cudaFree(d_p1);
    cudaFree(d_dist);
    cudaFree(d_grad);
    cudaFree(d_hess);

    // Compare against the CPU reference values.
    CHECK(dist == Catch::Approx(9.0));

    // grad.head(3) = 2 * (p0 - p1) = (-2, -4, -4); grad.tail(3) = -head.
    const std::array<double, grad_size> expected_grad = { -2.0, -4.0, -4.0,
                                                          2.0,  4.0,  4.0 };
    for (int i = 0; i < grad_size; ++i) {
        CHECK(grad[i] == Catch::Approx(expected_grad[i]));
    }

    // hess = [[2 I, -2 I], [-2 I, 2 I]] with 3x3 identity blocks.
    for (int r = 0; r < grad_size; ++r) {
        for (int c = 0; c < grad_size; ++c) {
            double expected = 0.0;
            if (r == c) {
                expected = 2.0;
            } else if (r == c + dim || c == r + dim) {
                expected = -2.0;
            }
            CHECK(hess[r * grad_size + c] == Catch::Approx(expected));
        }
    }
}

#endif // IPC_TOOLKIT_WITH_CUDA
