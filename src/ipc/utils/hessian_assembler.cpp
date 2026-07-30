#include "hessian_assembler.hpp"

#include <ipc/utils/logger.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include <algorithm>
#include <vector>

namespace ipc {

namespace {
    void set_triplets(
        const Eigen::SparseMatrix<double>& M,
        std::vector<Eigen::Triplet<double>>& triplets,
        const size_t start_index)
    {
        using InnerIterator = Eigen::SparseMatrix<double>::InnerIterator;
        assert(start_index + M.nonZeros() <= triplets.size());
        int count = 0;
        for (int k = 0; k < M.outerSize(); ++k) {
            for (InnerIterator it(M, k); it; ++it) {
                assert(count < M.nonZeros());
                triplets[start_index + count++] =
                    Eigen::Triplet<double>(it.row(), it.col(), it.value());
            }
        }
    }
} // namespace

void TripletHessianAssembler::begin(
    const int ndof,
    const int dim,
    const size_t num_stencils,
    const StencilGetter& /*stencil*/)
{
    m_ndof = ndof;
    m_dim = dim;

    constexpr int MAX_TRIPLETS_SIZE = 10'000'000;
    const int buffer_size = std::min(MAX_TRIPLETS_SIZE, ndof);

    m_storage = std::make_unique<
        tbb::enumerable_thread_specific<LocalThreadMatStorage>>(
        LocalThreadMatStorage(buffer_size, ndof, ndof));
}

void TripletHessianAssembler::add_local_hessian(
    const MatrixMax12d& local_hess, const std::array<index_t, 4>& vertex_ids)
{
    assert(m_storage != nullptr);
    IPC_TOOLKIT_PROFILE_BLOCK("Map Local Hessian to Global Triplets");
    local_hessian_to_global_triplets(
        local_hess, vertex_ids, m_dim, *(m_storage->local().cache),
        m_ndof / m_dim);
}

Eigen::SparseMatrix<double> TripletHessianAssembler::get_matrix()
{
    assert(m_storage != nullptr);
    tbb::enumerable_thread_specific<LocalThreadMatStorage>& storage =
        *m_storage;

    Eigen::SparseMatrix<double> hess(m_ndof, m_ndof);
    if (storage.empty()) {
        return hess;
    }

    // Assemble the stiffness matrix by concatenating the triplets in each
    // local storage.

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Prune Local Storages");
        tbb::parallel_for_each(
            storage.begin(), storage.end(),
            [](const auto& local_storage) { local_storage.cache->prune(); });
    }

    // Prepares for parallel concatenation
    std::vector<size_t> offsets(storage.size());

    size_t index = 0;
    size_t triplet_count = 0;
    for (auto& local_storage : storage) {
        offsets[index++] = triplet_count;
        triplet_count += local_storage.cache->triplet_count();
    }

    std::vector<Eigen::Triplet<double>> triplets;

    if (triplet_count >= triplets.max_size()) {
        // Serial fallback version in case the vector of triplets cannot be
        // allocated
        logger().warn(
            "Unable to allocate sufficient memory for triplets. "
            "Falling back to serial assembly, which may impact performance. "
            "Consider reducing the problem size or optimizing memory usage.");
        // Serially merge local storages
        for (LocalThreadMatStorage& local_storage : storage) {
            hess += local_storage.cache->get_matrix(false); // will also prune
        }
        hess.makeCompressed();
        return hess;
    }

    // Allocate triplets
    {
        IPC_TOOLKIT_PROFILE_BLOCK("Allocate Triplets");
        triplets.resize(triplet_count);
    }

    // Parallel copy into triplets
    {
        IPC_TOOLKIT_PROFILE_BLOCK("Parallel Copy into Triplets");
        tbb::parallel_for(size_t(0), storage.size(), [&](size_t i) {
            const SparseMatrixCache& cache =
                dynamic_cast<const SparseMatrixCache&>(
                    *((storage.begin() + i)->cache));
            size_t offset = offsets[i];

            std::copy(
                cache.entries().begin(), cache.entries().end(),
                triplets.begin() + offset);
            offset += cache.entries().size();

            if (cache.mat().nonZeros() > 0) {
                set_triplets(cache.mat(), triplets, offset);
            }
        });
    }

    // Sort and assemble
    {
        IPC_TOOLKIT_PROFILE_BLOCK("Assemble Hessian from Triplets");
        hess.setFromTriplets(triplets.begin(), triplets.end());
    }

    return hess;
}

} // namespace ipc
