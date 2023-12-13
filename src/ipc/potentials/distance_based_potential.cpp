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
    // w * m(x) * f(d(x))
    // NOTE: can save a multiplication by checking if !contact.is_mollified()
    return contact.weight * contact.mollifier(x)
        * distance_based_potential(contact.compute_distance(x), contact.dmin);
}

VectorMax12d DistanceBasedPotential::potential_gradient(
    const Contact& contact, const VectorMax12d& x) const
{
    const double d = contact.compute_distance(x);                     // d(x)
    const VectorMax12d grad_d = contact.compute_distance_gradient(x); // ∇d(x)

    // f(d(x))
    const double f = distance_based_potential(d, contact.dmin);
    // f'(d(x))
    const double grad_f = distance_based_potential_gradient(d, contact.dmin);

    if (!contact.is_mollified()) {
        // ∇[f(d(x))] = f'(d(x)) * ∇d(x)
        return (contact.weight * grad_f) * grad_d;
    }

    const double m = contact.mollifier(x);                     // m(x)
    const VectorMax12d grad_m = contact.mollifier_gradient(x); // ∇m(x)

    // ∇[m(x) * f(d(x))] = f(d(x)) * ∇m(x) + m(x) * ∇ f(d(x))
    return (contact.weight * f) * grad_m
        + (contact.weight * m * grad_f) * grad_d;
}

MatrixMax12d DistanceBasedPotential::potential_hessian(
    const Contact& contact,
    const VectorMax12d& x,
    const bool project_hessian_to_psd) const
{
    const double d = contact.compute_distance(x);                     // d(x)
    const VectorMax12d grad_d = contact.compute_distance_gradient(x); // ∇d(x)
    const MatrixMax12d hess_d = contact.compute_distance_hessian(x); // ∇²d(x)

    // f'(d(x))
    const double grad_f = distance_based_potential_gradient(d, contact.dmin);
    // f"(d(x))
    const double hess_f = distance_based_potential_hessian(d, contact.dmin);

    MatrixMax12d hess;
    if (!contact.is_mollified()) {
        // ∇²[f(d(x))] = ∇(f'(d(x)) * ∇d(x))
        //             = f"(d(x)) * ∇d(x) * ∇d(x)ᵀ + f'(d(x)) * ∇²d(x)
        hess = (contact.weight * hess_f) * grad_d * grad_d.transpose()
            + (contact.weight * grad_f) * hess_d;
    } else {
        const double f = distance_based_potential(d, contact.dmin); // f(d(x))

        const double m = contact.mollifier(x);                     // m(x)
        const VectorMax12d grad_m = contact.mollifier_gradient(x); // ∇ m(x)
        const MatrixMax12d hess_m = contact.mollifier_hessian(x);  // ∇² m(x)

        const double weighted_m = contact.weight * m;

        // ∇f(d(x)) * ∇m(x)ᵀ
        const MatrixMax12d grad_f_grad_m =
            (contact.weight * grad_f) * grad_d * grad_m.transpose();

        // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * ∇b(d(x))]
        //                    = ∇²m(x) * b(d(x)) + ∇b(d(x)) * ∇m(x)ᵀ
        //                      + ∇m(x) * ∇b(d(x))ᵀ + m(x) * ∇²b(d(x))
        hess = (contact.weight * f) * hess_m + grad_f_grad_m
            + grad_f_grad_m.transpose()
            + (weighted_m * hess_f) * grad_d * grad_d.transpose()
            + (weighted_m * grad_f) * hess_d;
    }

    // Need to project entire hessian because w can be negative
    return project_hessian_to_psd ? project_to_psd(hess) : hess;
}

} // namespace ipc