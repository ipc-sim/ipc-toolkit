#include "psd_projection.cuh"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/cuda/device_utils.cuh>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace ipc::cuda {

namespace {

    constexpr int DIM = BatchedPSDProjection::DIM;
    constexpr int BLOCK_SIZE = BatchedPSDProjection::BLOCK_SIZE;

    /// @brief One thread per block: full eigendecomposition + eigenvalue
    /// clamp/abs + reconstruction, in place. Mirrors libuipc's make_spd.
    __global__ void psd_project_kernel(
        double* __restrict__ blocks,
        const int num_blocks,
        const PSDProjectionMethod method)
    {
        const int b = blockIdx.x * blockDim.x + threadIdx.x;
        if (b >= num_blocks) {
            return;
        }

        // Column-major DIM×DIM block owned by this thread.
        Eigen::Map<Eigen::Matrix<double, DIM, DIM>> H(
            blocks + static_cast<size_t>(b) * BLOCK_SIZE);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, DIM, DIM>> es;
        es.compute(H); // eigenvalues + eigenvectors (QR iteration)

        Eigen::Matrix<double, DIM, 1> w = es.eigenvalues();
#pragma unroll
        for (int i = 0; i < DIM; i++) {
            w[i] = (method == PSDProjectionMethod::CLAMP) ? fmax(w[i], 0.0)
                                                          : fabs(w[i]);
        }

        H = es.eigenvectors() * w.asDiagonal() * es.eigenvectors().transpose();
    }

} // namespace

void BatchedPSDProjection::project(
    double* d_blocks, const int num_blocks, const PSDProjectionMethod method)
{
    if (num_blocks == 0 || method == PSDProjectionMethod::NONE) {
        return;
    }

    psd_project_kernel<<<kernel_grid_size(num_blocks), KERNEL_BLOCK_SIZE>>>(
        d_blocks, num_blocks, method);
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
