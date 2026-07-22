#include "psd_projection.cuh"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/cuda/device_utils.cuh>
#include <ipc/utils/logger.hpp>

#include <string>

namespace ipc::cuda {

namespace {

    /// @brief Throw if a cuSOLVER call fails.
    inline void check_cusolver(cusolverStatus_t status, const int line)
    {
        if (status != CUSOLVER_STATUS_SUCCESS) {
            log_and_throw_error(
                "cuSOLVER error at psd_projection.cu:{} (status {})", line,
                int(status));
        }
    }
#define IPC_CUSOLVER_CHECK(expr) check_cusolver((expr), __LINE__)

    /// @brief Reconstruct each projected block A_proj = V g(W) Vᵀ in place.
    /// One thread per output entry (column-major): blocks[b][c*DIM+r].
    __global__ void reconstruct_kernel(
        const double* __restrict__ eigenvectors, // n × DIM×DIM (columns = V_k)
        const double* __restrict__ eigenvalues,  // n × DIM (ascending)
        const size_t num_blocks,
        const PSDProjectionMethod method,
        double* __restrict__ blocks)
    {
        constexpr int DIM = BatchedPSDProjection::DIM;
        constexpr int BLOCK_SIZE = BatchedPSDProjection::BLOCK_SIZE;

        const size_t tid =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (tid >= num_blocks * size_t(BLOCK_SIZE)) {
            return;
        }

        const size_t b = tid / BLOCK_SIZE;
        const int entry = int(tid % BLOCK_SIZE);
        const int r = entry % DIM;
        const int c = entry / DIM;

        const double* V = eigenvectors + b * BLOCK_SIZE;
        const double* W = eigenvalues + b * DIM;

        double sum = 0.0;
        for (int k = 0; k < DIM; k++) {
            double eig = W[k];
            if (method == PSDProjectionMethod::CLAMP) {
                eig = fmax(eig, 0.0);
            } else { // ABS
                eig = fabs(eig);
            }
            // V column k: V[k*DIM + row]
            sum += V[k * DIM + r] * eig * V[k * DIM + c];
        }
        blocks[b * BLOCK_SIZE + entry] = sum;
    }

} // namespace

BatchedPSDProjection::BatchedPSDProjection()
{
    IPC_CUSOLVER_CHECK(cusolverDnCreate(&m_handle));
    IPC_CUSOLVER_CHECK(cusolverDnCreateSyevjInfo(&m_syevj_params));
}

BatchedPSDProjection::~BatchedPSDProjection()
{
    if (m_d_eigenvalues) {
        cudaFree(m_d_eigenvalues);
    }
    if (m_d_eigenvectors) {
        cudaFree(m_d_eigenvectors);
    }
    if (m_d_work) {
        cudaFree(m_d_work);
    }
    if (m_d_info) {
        cudaFree(m_d_info);
    }
    if (m_syevj_params) {
        cusolverDnDestroySyevjInfo(m_syevj_params);
    }
    if (m_handle) {
        cusolverDnDestroy(m_handle);
    }
}

void BatchedPSDProjection::ensure_capacity(const int num_blocks)
{
    if (num_blocks <= m_capacity_blocks) {
        return;
    }

    if (m_d_eigenvalues) {
        cudaFree(m_d_eigenvalues);
    }
    if (m_d_eigenvectors) {
        cudaFree(m_d_eigenvectors);
    }
    if (m_d_info) {
        cudaFree(m_d_info);
    }

    IPC_TOOLKIT_CUDA_CHECK(cudaMalloc(
        &m_d_eigenvalues, sizeof(double) * size_t(num_blocks) * DIM));
    IPC_TOOLKIT_CUDA_CHECK(cudaMalloc(
        &m_d_eigenvectors, sizeof(double) * size_t(num_blocks) * BLOCK_SIZE));
    IPC_TOOLKIT_CUDA_CHECK(
        cudaMalloc(&m_d_info, sizeof(int) * size_t(num_blocks)));

    // Query and (re)allocate the cuSOLVER workspace for this batch size.
    int lwork = 0;
    IPC_CUSOLVER_CHECK(cusolverDnDsyevjBatched_bufferSize(
        m_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, DIM,
        m_d_eigenvectors, DIM, m_d_eigenvalues, &lwork, m_syevj_params,
        num_blocks));
    if (m_d_work) {
        cudaFree(m_d_work);
    }
    IPC_TOOLKIT_CUDA_CHECK(
        cudaMalloc(&m_d_work, sizeof(double) * size_t(lwork)));

    m_lwork = lwork;
    m_capacity_blocks = num_blocks;
}

void BatchedPSDProjection::project(
    double* d_blocks, const int num_blocks, const PSDProjectionMethod method)
{
    if (num_blocks == 0 || method == PSDProjectionMethod::NONE) {
        return;
    }

    ensure_capacity(num_blocks);

    // syevj overwrites its input with the eigenvectors, so decompose a copy and
    // keep d_blocks as the destination for the reconstructed projection.
    IPC_TOOLKIT_CUDA_CHECK(cudaMemcpy(
        m_d_eigenvectors, d_blocks,
        sizeof(double) * size_t(num_blocks) * BLOCK_SIZE,
        cudaMemcpyDeviceToDevice));

    IPC_CUSOLVER_CHECK(cusolverDnDsyevjBatched(
        m_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, DIM,
        m_d_eigenvectors, DIM, m_d_eigenvalues, m_d_work, m_lwork, m_d_info,
        m_syevj_params, num_blocks));

    const size_t n_entries = size_t(num_blocks) * BLOCK_SIZE;
    reconstruct_kernel<<<kernel_grid_size(n_entries), KERNEL_BLOCK_SIZE>>>(
        m_d_eigenvectors, m_d_eigenvalues, size_t(num_blocks), method,
        d_blocks);
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
