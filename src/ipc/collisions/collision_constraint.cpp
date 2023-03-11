#include "collision_constraint.hpp"

#include <ipc/barrier/barrier.hpp>

namespace ipc {

double CollisionConstraint::compute_potential(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    const double distance =
        compute_distance(positions, edges, faces); // Squared distance
    return weight
        * barrier(
               distance - minimum_distance * minimum_distance,
               2 * minimum_distance * dhat + dhat * dhat);
}

VectorMax12d CollisionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    // ∇b(d(x)) = b'(d(x)) * ∇d(x)
    const double distance = compute_distance(positions, edges, faces);
    const VectorMax12d distance_grad =
        compute_distance_gradient(positions, edges, faces);

    const double grad_b = barrier_gradient(
        distance - minimum_distance * minimum_distance,
        2 * minimum_distance * dhat + dhat * dhat);
    return weight * grad_b * distance_grad;
}

MatrixMax12d CollisionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const bool project_hessian_to_psd) const
{
    const double dhat_squared = dhat * dhat;
    const double min_dist_squrared = minimum_distance * minimum_distance;

    // ∇²[b(d(x))] = ∇(b'(d(x)) * ∇d(x))
    //             = b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)

    const double distance = compute_distance(positions, edges, faces);
    const VectorMax12d distance_grad =
        compute_distance_gradient(positions, edges, faces);
    const MatrixMax12d distance_hess =
        compute_distance_hessian(positions, edges, faces);

    const double grad_b = barrier_gradient(
        distance - min_dist_squrared,
        2 * minimum_distance * dhat + dhat_squared);
    const double hess_b = barrier_hessian(
        distance - min_dist_squrared,
        2 * minimum_distance * dhat + dhat_squared);

    // b"(x) ≥ 0 ⟹ b"(x) * ∇d(x) * ∇d(x)ᵀ is PSD
    assert(hess_b >= 0);
    MatrixMax12d term1 = hess_b * distance_grad * distance_grad.transpose();
    MatrixMax12d term2 = grad_b * distance_hess;
    if (project_hessian_to_psd) {
        term2 = project_to_psd(term2);
    }

    return weight * (term1 + term2);
}

} // namespace ipc
