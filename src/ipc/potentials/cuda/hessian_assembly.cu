#include "hessian_assembly.cuh"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/cuda/device_utils.cuh>

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <vector>

namespace ipc::cuda {

namespace {

    constexpr int PADDED_DIM = 12; // block slot leading dimension
    constexpr int BLOCK_SLOT = PADDED_DIM * PADDED_DIM;

    /// Number of real dofs of each collision = 3 × (#non-negative vertex ids).
    __global__ void compute_ndof_kernel(
        const index_t* __restrict__ vid0,
        const index_t* __restrict__ vid1,
        const index_t* __restrict__ vid2,
        const index_t* __restrict__ vid3,
        const size_t n,
        int* __restrict__ ndof,
        int* __restrict__ entry_count)
    {
        const size_t i =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
        int nv = 0;
        nv += vid0[i] >= 0;
        nv += vid1[i] >= 0;
        nv += vid2[i] >= 0;
        nv += vid3[i] >= 0;
        const int d = 3 * nv;
        ndof[i] = d;
        entry_count[i] = d * d;
    }

    /// Emit each block entry as a (CSC key, value) pair. One thread per
    /// collision writes its ndof² entries at entry_offset[i].
    __global__ void coo_gen_kernel(
        const double* __restrict__ blocks,
        const index_t* __restrict__ vid0,
        const index_t* __restrict__ vid1,
        const index_t* __restrict__ vid2,
        const index_t* __restrict__ vid3,
        const int* __restrict__ ndof,
        const long long* __restrict__ entry_offset,
        const size_t n,
        const index_t num_vertices,
        long long* __restrict__ keys,
        double* __restrict__ vals)
    {
        const size_t i =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }

        const index_t vids[4] = { vid0[i], vid1[i], vid2[i], vid3[i] };
        const int d = ndof[i];
        const long long base = entry_offset[i];
        const long long n_dofs = 3ll * num_vertices;
        const double* block = blocks + i * BLOCK_SLOT;

        long long e = base;
        for (int lc = 0; lc < d; lc++) { // local column = 3*b + l
            const index_t gcol =
                global_dof_index(vids[lc / 3], lc % 3, num_vertices);
            for (int lr = 0; lr < d; lr++) { // local row = 3*a + k
                const index_t grow =
                    global_dof_index(vids[lr / 3], lr % 3, num_vertices);
                keys[e] = static_cast<long long>(gcol) * n_dofs + grow;
                vals[e] = block[lc * PADDED_DIM + lr];
                e++;
            }
        }
    }

    struct KeyToRow {
        long long n_dofs;
        __host__ __device__ int operator()(const long long key) const
        {
            return static_cast<int>(key % n_dofs);
        }
    };

    struct KeyToCol {
        long long n_dofs;
        __host__ __device__ int operator()(const long long key) const
        {
            return static_cast<int>(key / n_dofs);
        }
    };

} // namespace

Eigen::SparseMatrix<double> assemble_sparse_hessian(
    const double* d_blocks,
    const index_t* d_vertex_id_0,
    const index_t* d_vertex_id_1,
    const index_t* d_vertex_id_2,
    const index_t* d_vertex_id_3,
    const size_t n,
    const index_t num_vertices)
{
    const int n_dofs = 3 * static_cast<int>(num_vertices);
    Eigen::SparseMatrix<double> H(n_dofs, n_dofs);
    if (n == 0) {
        return H;
    }

    // Per-collision dof count and entry count.
    thrust::device_vector<int> d_ndof(n), d_entry_count(n);
    compute_ndof_kernel<<<kernel_grid_size(n), KERNEL_BLOCK_SIZE>>>(
        d_vertex_id_0, d_vertex_id_1, d_vertex_id_2, d_vertex_id_3, n,
        thrust::raw_pointer_cast(d_ndof.data()),
        thrust::raw_pointer_cast(d_entry_count.data()));
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    // Exclusive-scan the entry counts to get per-collision offsets and E.
    thrust::device_vector<long long> d_entry_offset(n);
    thrust::exclusive_scan(
        d_entry_count.begin(), d_entry_count.end(), d_entry_offset.begin(),
        0ll);
    const long long last_count = d_entry_count.back();
    const long long E = d_entry_offset.back() + last_count;
    if (E == 0) {
        return H;
    }

    // Generate (key, value) COO entries.
    thrust::device_vector<long long> d_keys(E);
    thrust::device_vector<double> d_vals(E);
    coo_gen_kernel<<<kernel_grid_size(n), KERNEL_BLOCK_SIZE>>>(
        d_blocks, d_vertex_id_0, d_vertex_id_1, d_vertex_id_2, d_vertex_id_3,
        thrust::raw_pointer_cast(d_ndof.data()),
        thrust::raw_pointer_cast(d_entry_offset.data()), n, num_vertices,
        thrust::raw_pointer_cast(d_keys.data()),
        thrust::raw_pointer_cast(d_vals.data()));
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    // Sort by (col, row) key and sum duplicates -> CSC-ordered unique entries.
    thrust::sort_by_key(d_keys.begin(), d_keys.end(), d_vals.begin());

    thrust::device_vector<long long> d_ukeys(E);
    thrust::device_vector<double> d_uvals(E);
    const auto new_end = thrust::reduce_by_key(
        d_keys.begin(), d_keys.end(), d_vals.begin(), d_ukeys.begin(),
        d_uvals.begin());
    const size_t nnz = new_end.first - d_ukeys.begin();

    // Split unique keys into CSC inner (row) indices and column indices.
    const long long n_dofs_ll = n_dofs;
    thrust::device_vector<int> d_inner(nnz), d_col(nnz);
    thrust::transform(
        d_ukeys.begin(), d_ukeys.begin() + nnz, d_inner.begin(),
        KeyToRow { n_dofs_ll });
    thrust::transform(
        d_ukeys.begin(), d_ukeys.begin() + nnz, d_col.begin(),
        KeyToCol { n_dofs_ll });

    // CSC outer pointers: outer[j] = #nonzeros with col < j.
    thrust::device_vector<int> d_outer(n_dofs + 1);
    thrust::lower_bound(
        d_col.begin(), d_col.begin() + nnz, thrust::counting_iterator<int>(0),
        thrust::counting_iterator<int>(n_dofs + 1), d_outer.begin());

    // Copy CSC structure + values to host and build the Eigen matrix.
    std::vector<int> h_outer(n_dofs + 1), h_inner(nnz);
    std::vector<double> h_vals(nnz);
    thrust::copy(d_outer.begin(), d_outer.end(), h_outer.begin());
    thrust::copy(d_inner.begin(), d_inner.begin() + nnz, h_inner.begin());
    thrust::copy(d_uvals.begin(), d_uvals.begin() + nnz, h_vals.begin());

    return Eigen::Map<const Eigen::SparseMatrix<double>>(
        n_dofs, n_dofs, static_cast<Eigen::Index>(nnz), h_outer.data(),
        h_inner.data(), h_vals.data());
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
