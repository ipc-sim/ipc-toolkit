#include "potential.hpp"

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_reduce.h>

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

template <class TCollisions>
double Potential<TCollisions>::operator()(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    assert(X.rows() == mesh.num_vertices());

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial_sum) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Quadrature weight is premultiplied by local potential
                partial_sum += (*this)(
                    collisions[i],
                    collisions[i].dof(X, mesh.edges(), mesh.faces()));
            }
            return partial_sum;
        },
        [](double a, double b) { return a + b; });
}

template <class TCollisions>
Eigen::VectorXd Potential<TCollisions>::gradient(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(X.size());
    }

    const int dim = X.cols();

    tbb::combinable<Eigen::VectorXd> grad(Eigen::VectorXd::Zero(X.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const TCollision& collision = collisions[i];

                const VectorMaxNd local_grad = this->gradient(
                    collision, collision.dof(X, mesh.edges(), mesh.faces()));

                const std::array<index_t, TCollision::STENCIL_SIZE> vids =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(
                    local_grad, vids, dim, grad.local());
            }
        });

    return grad.combine([](const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
        return a + b;
    });
}

template <class TCollisions>
Eigen::SparseMatrix<double> Potential<TCollisions>::hessian(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(X.size(), X.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = X.cols();
    const int ndof = X.size();

    constexpr int MAX_TRIPLETS_SIZE = 10'000'000;
    const int buffer_size = std::min(MAX_TRIPLETS_SIZE, ndof);

    tbb::enumerable_thread_specific<LocalThreadMatStorage> storage(
        LocalThreadMatStorage(buffer_size, ndof, ndof));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {

                const TCollision& collision = collisions[i];

                const MatrixMaxNd local_hess = this->hessian(
                    collisions[i], collisions[i].dof(X, edges, faces),
                    project_hessian_to_psd);

                const std::array<index_t, TCollision::STENCIL_SIZE> vids =
                    collision.vertex_ids(edges, faces);

                local_hessian_to_global_triplets(
                    local_hess, vids, dim, *(hess_triplets.cache));
            }
        });

    if (storage.empty()) {
        return Eigen::SparseMatrix<double>();
    }

    // Assemble the stiffness matrix by concatenating the tuples in each local
    // storage

    tbb::parallel_for_each(
        storage.begin(), storage.end(),
        [](const auto& local_storage) { local_storage.cache->prune(); });

    // Prepares for parallel concatenation
    std::vector<size_t> offsets(storage.size());

    size_t index = 0;
    size_t triplet_count = 0;
    for (auto& local_storage : storage) {
        offsets[index++] = triplet_count;
        triplet_count += local_storage.cache->triplet_count();
    }

    std::vector<Eigen::Triplet<double>> triplets;

    Eigen::SparseMatrix<double> hess(ndof, ndof);
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

    triplets.resize(triplet_count);

    // Parallel copy into triplets
    tbb::parallel_for(size_t(0), storage.size(), [&](size_t i) {
        const SparseMatrixCache& cache = dynamic_cast<const SparseMatrixCache&>(
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

    // Sort and assemble
    hess.setFromTriplets(triplets.begin(), triplets.end());

    return hess;
}

template class Potential<NormalCollisions>;
template class Potential<TangentialCollisions>;

} // namespace ipc