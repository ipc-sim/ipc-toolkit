#include "distance_based_potential.hpp"

#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

double DistanceBasedPotential::operator()(
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
                local_potential += potential(
                    contacts[i],
                    contacts[i].dof(vertices, mesh.edges(), mesh.faces()));
            }
        });

    return storage.combine([](double a, double b) { return a + b; });
}

Eigen::VectorXd DistanceBasedPotential::gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Contacts& contacts) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (contacts.empty()) {
        return Eigen::VectorXd::Zero(vertices.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = vertices.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(vertices.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), contacts.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_grad = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                const VectorMax12d x =
                    contacts[i].dof(vertices, mesh.edges(), mesh.faces());
                local_gradient_to_global_gradient(
                    potential_gradient(contacts[i], x),
                    contacts[i].vertex_ids(edges, faces), dim, local_grad);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

Eigen::SparseMatrix<double> DistanceBasedPotential::hessian(
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

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), contacts.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const VectorMax12d x =
                    contacts[i].dof(vertices, mesh.edges(), mesh.faces());
                local_hessian_to_global_triplets(
                    potential_hessian(contacts[i], x, project_hessian_to_psd),
                    contacts[i].vertex_ids(edges, faces), dim,
                    local_hess_triplets);
            }
        });

    Eigen::SparseMatrix<double> hess(vertices.size(), vertices.size());
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(
            vertices.size(), vertices.size());
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

// ----------------------------------------------------------------------------

double DistanceBasedPotential::potential(
    const Contact& contact, const VectorMax12d& x) const
{
    // w * f(d(x))
    return contact.weight
        * distance_based_potential(contact.compute_distance(x), contact.dmin);
}

VectorMax12d DistanceBasedPotential::potential_gradient(
    const Contact& contact, const VectorMax12d& x) const
{
    // ∇f(d(x)) = f'(d(x)) * ∇d(x)
    return contact.weight
        * distance_based_potential_gradient(
               contact.compute_distance(x), contact.dmin)
        * contact.compute_distance_gradient(x);
}

MatrixMax12d DistanceBasedPotential::potential_hessian(
    const Contact& contact,
    const VectorMax12d& x,
    const bool project_hessian_to_psd) const
{
    // ∇²[f(d(x))] = ∇(f'(d(x)) * ∇d(x))
    //             = f"(d(x)) * ∇d(x) * ∇d(x)ᵀ + f'(d(x)) * ∇²d(x)

    const double d = contact.compute_distance(x);
    const VectorMax12d grad_d = contact.compute_distance_gradient(x);
    const MatrixMax12d hess_d = contact.compute_distance_hessian(x);

    const double grad_f = distance_based_potential_gradient(d, contact.dmin);
    const double hess_f = distance_based_potential_hessian(d, contact.dmin);

    // f"(d(x)) ≥ 0 ⟹ f"(d(x)) * ∇d(x) * ∇d(x)ᵀ is PSD
    assert(hess_f >= 0);
    MatrixMax12d term1 = hess_f * grad_d * grad_d.transpose();
    MatrixMax12d term2 = grad_f * hess_d;
    if (project_hessian_to_psd) {
        term2 = project_to_psd(term2);
    }

    return contact.weight * (term1 + term2);
}

} // namespace ipc