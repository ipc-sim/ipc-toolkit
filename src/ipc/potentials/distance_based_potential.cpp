#include "distance_based_potential.hpp"

namespace ipc {

// -- Cumulative methods ---------------------------------------------------

Eigen::SparseMatrix<double> DistanceBasedPotential::shape_derivative(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Contacts& contacts) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (contacts.empty()) {
        return Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
    }

    const Eigen::MatrixXd& rest_positions = mesh.rest_positions();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int ndof = vertices.size();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), contacts.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                this->shape_derivative(
                    contacts[i], contacts[i].vertex_ids(edges, faces),
                    contacts[i].dof(rest_positions, edges, faces),
                    contacts[i].dof(vertices, edges, faces), local_triplets);
            }
        });

    Eigen::SparseMatrix<double> shape_derivative(ndof, ndof);
    for (const auto& local_triplets : storage) {
        Eigen::SparseMatrix<double> local_shape_derivative(ndof, ndof);
        local_shape_derivative.setFromTriplets(
            local_triplets.begin(), local_triplets.end());
        shape_derivative += local_shape_derivative;
    }
    return shape_derivative;
}

// -- Single contact methods -----------------------------------------------

double DistanceBasedPotential::operator()(
    const Contact& contact, const VectorMax12d& positions) const
{
    // w * m(x) * f(d(x))
    // NOTE: can save a multiplication by checking if !contact.is_mollified()
    const double d = contact.compute_distance(positions);
    return contact.weight * contact.mollifier(positions)
        * distance_based_potential(d, contact.dmin);
}

VectorMax12d DistanceBasedPotential::gradient(
    const Contact& contact, const VectorMax12d& positions) const
{
    // d(x)
    const double d = contact.compute_distance(positions);
    // ∇d(x)
    const VectorMax12d grad_d = contact.compute_distance_gradient(positions);

    // f(d(x))
    const double f = distance_based_potential(d, contact.dmin);
    // f'(d(x))
    const double grad_f = distance_based_potential_gradient(d, contact.dmin);

    if (!contact.is_mollified()) {
        // ∇[f(d(x))] = f'(d(x)) * ∇d(x)
        return (contact.weight * grad_f) * grad_d;
    }

    const double m = contact.mollifier(positions);                     // m(x)
    const VectorMax12d grad_m = contact.mollifier_gradient(positions); // ∇m(x)

    // ∇[m(x) * f(d(x))] = f(d(x)) * ∇m(x) + m(x) * ∇ f(d(x))
    return (contact.weight * f) * grad_m
        + (contact.weight * m * grad_f) * grad_d;
}

MatrixMax12d DistanceBasedPotential::hessian(
    const Contact& contact,
    const VectorMax12d& positions,
    const bool project_hessian_to_psd) const
{
    // d(x)
    const double d = contact.compute_distance(positions);
    // ∇d(x)
    const VectorMax12d grad_d = contact.compute_distance_gradient(positions);
    // ∇²d(x)
    const MatrixMax12d hess_d = contact.compute_distance_hessian(positions);

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

        // m(x)
        const double m = contact.mollifier(positions);
        // ∇ m(x)
        const VectorMax12d grad_m = contact.mollifier_gradient(positions);
        // ∇² m(x)
        const MatrixMax12d hess_m = contact.mollifier_hessian(positions);

        const double weighted_m = contact.weight * m;

        // ∇f(d(x)) * ∇m(x)ᵀ
        const MatrixMax12d grad_f_grad_m =
            (contact.weight * grad_f) * grad_d * grad_m.transpose();

        // ∇²[m(x) * f(d(x))] = ∇[∇m(x) * f(d(x)) + m(x) * ∇f(d(x))]
        //                    = ∇²m(x) * f(d(x)) + ∇f(d(x)) * ∇m(x)ᵀ
        //                      + ∇m(x) * ∇f(d(x))ᵀ + m(x) * ∇²f(d(x))
        hess = (contact.weight * f) * hess_m + grad_f_grad_m
            + grad_f_grad_m.transpose()
            + (weighted_m * hess_f) * grad_d * grad_d.transpose()
            + (weighted_m * grad_f) * hess_d;
    }

    // Need to project entire hessian because w can be negative
    return project_hessian_to_psd ? project_to_psd(hess) : hess;
}

void DistanceBasedPotential::shape_derivative(
    const Contact& contact,
    const std::array<long, 4>& vertex_ids,
    const VectorMax12d& rest_positions, // = x̄
    const VectorMax12d& positions,      // = x̄ + u
    std::vector<Eigen::Triplet<double>>& out) const
{
    assert(rest_positions.size() == positions.size());

    const int dim = positions.size() / contact.num_vertices();
    assert(positions.size() % contact.num_vertices() == 0);

    // Compute:
    // ∇ₓ (w ∇ᵤf(d(x̄+u))) = (∇ₓw)(∇ᵤf(d(x̄+u)))ᵀ + w ∇ₓ∇ᵤf(d(x̄+u))
    //                         (first term)        (second term)

    // First term:
    if (contact.weight_gradient.size() <= 0) {
        throw std::runtime_error(
            "Shape derivative is not computed for contact constraint!");
    }

    if (contact.weight_gradient.nonZeros()) {
        VectorMax12d grad_b = gradient(contact, positions);
        assert(contact.weight != 0);
        grad_b.array() /= contact.weight; // remove weight

        for (int i = 0; i < contact.num_vertices(); i++) {
            for (int d = 0; d < dim; d++) {
                using Itr = Eigen::SparseVector<double>::InnerIterator;
                for (Itr j(contact.weight_gradient); j; ++j) {
                    out.emplace_back(
                        vertex_ids[i] * dim + d, j.index(),
                        grad_b[dim * i + d] * j.value());
                }
            }
        }
    }

    // Second term:
    MatrixMax12d local_hess;
    if (!contact.is_mollified()) {
        // w ∇ₓ∇ᵤf = w ∇ᵤ²f
        local_hess =
            hessian(contact, positions, /*project_hessian_to_psd=*/false);
    } else {
        // d(x̄+u)
        const double d = contact.compute_distance(positions);
        // ∇d(x̄+u)
        const VectorMax12d grad_d =
            contact.compute_distance_gradient(positions);
        // ∇²d(x̄+u)
        const MatrixMax12d hess_d = contact.compute_distance_hessian(positions);

        // f(d(x̄+u))
        const double f =
            contact.weight * distance_based_potential(d, contact.dmin);
        // ∇ᵤ f(d(x̄+u))
        const Vector12d gradu_f =
            (contact.weight
             * distance_based_potential_gradient(d, contact.dmin))
            * grad_d;
        // ∇ᵤ² f(d(x̄+u))
        const Matrix12d hessu_f =
            (contact.weight * distance_based_potential_hessian(d, contact.dmin))
                * grad_d * grad_d.transpose()
            + (contact.weight
               * distance_based_potential_gradient(d, contact.dmin))
                * hess_d;

        // ε(x̄)
        const double eps_x = contact.mollifier_threshold(rest_positions);

        // m(x̄,u)
        const double m = contact.mollifier(positions, eps_x);
        // ∇ᵤ m(x̄,u)
        const Vector12d gradu_m = contact.mollifier_gradient(positions, eps_x);
        // ∇ₓ m(x̄,u)
        const Vector12d gradx_m =
            contact.mollifier_gradient_wrt_x(rest_positions, positions);
        // ∇ₓ∇ᵤ m(x̄,u)
        const Matrix12d jac_m = contact.mollifier_gradient_jacobian_wrt_x(
            rest_positions, positions);

        // Only compute the second term of the shape derivative
        // ∇ₓ (f ∇ᵤm + m ∇ᵤf) = f ∇ₓ∇ᵤm + (∇ₓm)(∇ᵤf)ᵀ + (∇ᵤf)(∇ᵤm)ᵀ + m ∇ᵤ²f
        local_hess = f * jac_m + gradx_m * gradu_f.transpose()
            + gradu_f * gradu_m.transpose() + m * hessu_f;
    }

    local_hessian_to_global_triplets(local_hess, vertex_ids, dim, out);
}

} // namespace ipc