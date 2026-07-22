// Verifies the cuSOLVER-backed BatchedPSDProjection matches Eigen's
// project_to_psd on random symmetric matrices (dims 6, 9, 12 -- the barrier
// hessian stencil sizes -- zero-padded to 12x12).

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/potentials/cuda/psd_projection.cuh>
#include <ipc/utils/eigen_ext.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cuda_runtime.h>

#include <Eigen/Eigenvalues>

#include <vector>

namespace {

bool has_cuda_device()
{
    int n = 0;
    return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

// Deterministic symmetric matrix with a mix of positive and negative
// eigenvalues, embedded (zero-padded) in the top-left ndof block of a 12x12.
Eigen::Matrix<double, 12, 12> random_symmetric(const int ndof, const int seed)
{
    Eigen::Matrix<double, 12, 12> M;
    M.setZero();
    // Simple LCG so the test is reproducible without <random> ordering issues.
    uint64_t s = uint64_t(seed) * 2654435761u + 1013904223u;
    auto next = [&]() {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        return double(int32_t(s >> 33)) / double(1u << 30); // ~[-2, 2]
    };
    Eigen::MatrixXd A(ndof, ndof);
    for (int i = 0; i < ndof; i++) {
        for (int j = 0; j < ndof; j++) {
            A(i, j) = next();
        }
    }
    M.topLeftCorner(ndof, ndof) = A + A.transpose(); // symmetric, indefinite
    return M;
}

} // namespace

TEST_CASE(
    "GPU batched PSD projection", "[potential][psd_projection][cuda][gpu]")
{
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

    const ipc::PSDProjectionMethod method = GENERATE(
        ipc::PSDProjectionMethod::CLAMP, ipc::PSDProjectionMethod::ABS);
    CAPTURE(static_cast<int>(method));

    // A batch mixing the three stencil sizes.
    constexpr int DIM = ipc::cuda::BatchedPSDProjection::DIM;
    constexpr int BLOCK = ipc::cuda::BatchedPSDProjection::BLOCK_SIZE;
    const int ndofs[] = { 6, 9, 12 };
    const int num_blocks = 30;

    std::vector<Eigen::Matrix<double, 12, 12>> mats(num_blocks);
    std::vector<double> host_blocks(size_t(num_blocks) * BLOCK);
    for (int b = 0; b < num_blocks; b++) {
        mats[b] = random_symmetric(ndofs[b % 3], b + 1);
        // Column-major copy into the batch buffer.
        for (int c = 0; c < DIM; c++) {
            for (int r = 0; r < DIM; r++) {
                host_blocks[size_t(b) * BLOCK + c * DIM + r] = mats[b](r, c);
            }
        }
    }

    double* d_blocks = nullptr;
    REQUIRE(
        cudaMalloc(&d_blocks, sizeof(double) * host_blocks.size())
        == cudaSuccess);
    REQUIRE(
        cudaMemcpy(
            d_blocks, host_blocks.data(), sizeof(double) * host_blocks.size(),
            cudaMemcpyHostToDevice)
        == cudaSuccess);

    ipc::cuda::BatchedPSDProjection projector;
    projector.project(d_blocks, num_blocks, method);

    REQUIRE(
        cudaMemcpy(
            host_blocks.data(), d_blocks, sizeof(double) * host_blocks.size(),
            cudaMemcpyDeviceToHost)
        == cudaSuccess);
    cudaFree(d_blocks);

    for (int b = 0; b < num_blocks; b++) {
        const int ndof = ndofs[b % 3];
        CAPTURE(b, ndof);

        // Reference: Eigen's project_to_psd on the real ndof×ndof submatrix.
        const ipc::MatrixMax12d expected = ipc::project_to_psd(
            ipc::MatrixMax12d(mats[b].topLeftCorner(ndof, ndof)), method);

        Eigen::Matrix<double, 12, 12> got;
        for (int c = 0; c < DIM; c++) {
            for (int r = 0; r < DIM; r++) {
                got(r, c) = host_blocks[size_t(b) * BLOCK + c * DIM + r];
            }
        }

        // Real block matches the reference; the padded region stays zero.
        CHECK(
            (got.topLeftCorner(ndof, ndof) - expected).norm()
            <= 1e-9 * (1 + expected.norm()));
        CHECK(got.bottomRows(DIM - ndof).norm() <= 1e-9);
        CHECK(got.rightCols(DIM - ndof).norm() <= 1e-9);
    }
}

#endif // IPC_TOOLKIT_WITH_CUDA
