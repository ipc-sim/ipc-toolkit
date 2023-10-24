#include "collision_constraint.hpp"

#include <ipc/barrier/barrier.hpp>

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

} // namespace ipc
