#pragma once

#include "potential.hpp"

#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

template <class Contacts>
double Potential<Contacts>::operator()(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Contacts& contacts) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (contacts.empty()) {
        return 0;
    }

    tbb::enumerable_thread_specific<double> storage(0);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), contacts.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Quadrature weight is premultiplied by compute_potential
                local_potential += (*this)(
                    contacts[i],
                    contacts[i].dof(vertices, mesh.edges(), mesh.faces()));
            }
        });

    return storage.combine([](double a, double b) { return a + b; });
}

template <class Contacts>
Eigen::VectorXd Potential<Contacts>::gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Contacts& contacts) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (contacts.empty()) {
        return Eigen::VectorXd::Zero(vertices.size());
    }

    const int dim = vertices.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(vertices.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), contacts.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& global_grad = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const Contact& contact = contacts[i];

                const VectorMax12d local_grad = this->gradient(
                    contact, contact.dof(vertices, mesh.edges(), mesh.faces()));

                const std::array<long, 4> vids =
                    contact.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(
                    local_grad, vids, dim, global_grad);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

template <class Contacts>
Eigen::SparseMatrix<double> Potential<Contacts>::hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Contacts& contacts,
    const bool project_hessian_to_psd) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (contacts.empty()) {
        return Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = vertices.cols();
    const int ndof = vertices.size();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), contacts.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const Contact& contact = contacts[i];

                const MatrixMax12d local_hess = this->hessian(
                    contacts[i], contacts[i].dof(vertices, edges, faces),
                    project_hessian_to_psd);

                const std::array<long, 4> vids =
                    contact.vertex_ids(edges, faces);

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