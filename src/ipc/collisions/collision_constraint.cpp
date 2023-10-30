#include "collision_constraint.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/local_to_global.hpp>

namespace ipc {

CollisionConstraint::CollisionConstraint(
    const double _weight, const Eigen::SparseVector<double>& _weight_gradient)
    : weight(_weight)
    , weight_gradient(_weight_gradient)
{
}

double CollisionConstraint::compute_potential(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    // Squared distance
    const double d = compute_distance(vertices, edges, faces);
    return weight * barrier(d - dmin * dmin, 2 * dmin * dhat + dhat * dhat);
}

VectorMax12d CollisionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    // ∇b(d(x)) = b'(d(x)) * ∇d(x)
    const VectorMax12d positions = dof(vertices, edges, faces);
    const double d = compute_distance(positions);
    const VectorMax12d grad_d = compute_distance_gradient(positions);

    const double grad_b =
        barrier_gradient(d - dmin * dmin, 2 * dmin * dhat + dhat * dhat);
    return weight * grad_b * grad_d;
}

MatrixMax12d CollisionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const bool project_hessian_to_psd) const
{
    const double adjusted_dhat = 2 * dmin * dhat + dhat * dhat;
    const double dmin_sq = dmin * dmin;

    // ∇²[b(d(x))] = ∇(b'(d(x)) * ∇d(x))
    //             = b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)

    const VectorMax12d positions = dof(vertices, edges, faces);
    const double d = compute_distance(positions);
    const VectorMax12d grad_d = compute_distance_gradient(positions);
    const MatrixMax12d hess_d = compute_distance_hessian(positions);

    const double grad_b = barrier_gradient(d - dmin_sq, adjusted_dhat);
    const double hess_b = barrier_hessian(d - dmin_sq, adjusted_dhat);

    // b"(x) ≥ 0 ⟹ b"(x) * ∇d(x) * ∇d(x)ᵀ is PSD
    assert(hess_b >= 0);
    MatrixMax12d term1 = hess_b * grad_d * grad_d.transpose();
    MatrixMax12d term2 = grad_b * hess_d;
    if (project_hessian_to_psd) {
        term2 = project_to_psd(term2);
    }

    return weight * (term1 + term2);
}

void CollisionConstraint::compute_shape_derivative(
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    std::vector<Eigen::Triplet<double>>& triplets) const
{
    compute_shape_derivative_first_term(
        rest_positions, vertices, edges, faces, dhat, triplets);

    local_hessian_to_global_triplets(
        compute_shape_derivative_second_term(
            rest_positions, vertices, edges, faces, dhat),
        vertex_ids(edges, faces), vertices.cols(), triplets);
}

void CollisionConstraint::compute_shape_derivative_first_term(
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    std::vector<Eigen::Triplet<double>>& triplets) const
{
    if (weight_gradient.size() != vertices.size()) {
        throw std::runtime_error(
            "Shape derivative is not computed for contact constraint!");
    }

    if (weight_gradient.nonZeros() == 0) {
        return;
    }

    const int dim = vertices.cols();

    VectorMax12d grad_b =
        compute_potential_gradient(vertices, edges, faces, dhat);
    assert(weight != 0);
    grad_b.array() /= weight; // remove weight

    const std::array<long, 4> ids = vertex_ids(edges, faces);
    for (int i = 0; i < num_vertices(); i++) {
        for (int d = 0; d < dim; d++) {
            using Itr = Eigen::SparseVector<double>::InnerIterator;
            for (Itr j(weight_gradient); j; ++j) {
                triplets.emplace_back(
                    ids[i] * dim + d, j.index(),
                    grad_b[dim * i + d] * j.value());
            }
        }
    }
}

MatrixMax12d CollisionConstraint::compute_shape_derivative_second_term(
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    // w ∇ₓ∇ᵤb = w ∇ᵤ²b
    return compute_potential_hessian(vertices, edges, faces, dhat, false);
}

} // namespace ipc
