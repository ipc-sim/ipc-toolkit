#include "local_hessian_assembly.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/local_to_global.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace ipc::cuda::internal {

void append_local_hessian_blocks(
    const double* blocks,
    const size_t slab_begin,
    const size_t slab_count,
    const std::vector<std::array<index_t, 4>>& vertex_ids,
    const PSDProjectionMethod project_hessian_to_psd,
    const index_t n_total_vertices,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(size_t(0), slab_count, [&](const size_t k) {
        std::vector<Eigen::Triplet<double>>& local_triplets = storage.local();
        const std::array<index_t, 4>& ids = vertex_ids[slab_begin + k];

        // Unused vertex id entries are -1.
        int n_verts = 0;
        while (n_verts < 4 && ids[n_verts] >= 0) {
            n_verts++;
        }
        const int ndof = 3 * n_verts;

        MatrixMax12d local_hessian = Eigen::Map<const Eigen::MatrixXd>(
            blocks + k * HESSIAN_BLOCK_SLOT_SIZE, ndof, ndof);

        local_hessian = project_to_psd(local_hessian, project_hessian_to_psd);

        local_hessian_to_global_triplets(
            local_hessian, ids, /*dim=*/3, local_triplets, n_total_vertices);
    });

    for (const std::vector<Eigen::Triplet<double>>& local_triplets : storage) {
        triplets.insert(
            triplets.end(), local_triplets.begin(), local_triplets.end());
    }
}

} // namespace ipc::cuda::internal

#endif // IPC_TOOLKIT_WITH_CUDA
