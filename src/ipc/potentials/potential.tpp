#pragma once

#include "potential.hpp"

#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

template <class TCollisions>
double Potential<TCollisions>::operator()(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X) const
{
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return 0;
    }

    tbb::enumerable_thread_specific<double> storage(0);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Quadrature weight is premultiplied by local potential
                local_potential += (*this)(
                    collisions[i],
                    collisions[i].dof(X, mesh.edges(), mesh.faces()));
            }
        });

    return storage.combine([](double a, double b) { return a + b; });
}

template <class TCollisions>
Eigen::VectorXd Potential<TCollisions>::gradient(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X) const
{
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(X.size());
    }

    const int dim = X.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(X.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& global_grad = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const TCollision& collision = collisions[i];

                const VectorMax12d local_grad = this->gradient(
                    collision, collision.dof(X, mesh.edges(), mesh.faces()));

                const std::array<long, 4> vids =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(
                    local_grad, vids, dim, global_grad);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

template <class TCollisions>
Eigen::SparseMatrix<double> Potential<TCollisions>::hessian(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const bool project_hessian_to_psd) const
{
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(X.size(), X.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = X.cols();
    const int ndof = X.size();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const TCollision& collision = collisions[i];

                const MatrixMax12d local_hess = this->hessian(
                    collisions[i], collisions[i].dof(X, edges, faces),
                    project_hessian_to_psd);

                const std::array<long, 4> vids =
                    collision.vertex_ids(edges, faces);

                local_hessian_to_global_triplets(
                    local_hess, vids, dim, hess_triplets);
            }
        });

    Eigen::SparseMatrix<double> hess(ndof, ndof);
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(ndof, ndof);
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

} // namespace ipc