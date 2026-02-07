#include "point_line.hpp"

namespace ipc {

Vector6d point_line_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector2d> p,
    Eigen::ConstRef<Eigen::Vector2d> e0,
    Eigen::ConstRef<Eigen::Vector2d> e1)
{
    const Eigen::Vector2d n = point_line_normal(p, e0, e1);
    const Eigen::Matrix<double, 2, 6> jac_n =
        point_line_normal_jacobian(p, e0, e1);

    Vector6d grad = jac_n.transpose() * (p - e0);
    grad.segment<2>(0) += n;
    grad.segment<2>(2) -= n;

    return grad;
}

Matrix6d point_line_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector2d> p,
    Eigen::ConstRef<Eigen::Vector2d> e0,
    Eigen::ConstRef<Eigen::Vector2d> e1)
{
    // Precompute normal's Jacobian and Hessian
    const Eigen::Matrix<double, 2, 6> jac_n =
        point_line_normal_jacobian(p, e0, e1);
    const Eigen::Matrix<double, 12, 6> hess_n =
        point_line_normal_hessian(p, e0, e1);

    // Vector from e0 to p
    const Eigen::Vector2d v = p - e0;

    Matrix6d hess;

    // ---------------------------------------------------------
    // 1. Tensor Contraction (Curvature Term)
    // ---------------------------------------------------------
    // Contract the normal Hessian (2x36) with vector v (2x1).
    // Result is 1x36, mapped to 6x6.
    for (int i = 0; i < 6; ++i) {
        hess.row(i) = hess_n.middleRows<2>(2 * i).transpose() * v;
    }

    // ---------------------------------------------------------
    // 2. Add Jacobian Terms (Product Rule Corrections)
    // ---------------------------------------------------------
    // Formula: H += (Jₙᵀ Jᵥ) + (Jᵥᵀ Jₙ)
    // Jᵥ w.r.t [p, e0, e1] is [I, -I, 0]

    // Extract 2x2 Jacobian blocks for e₀ and e₁
    // Note: ∇ n has columns 0-1 (p), 2-3 (e₀), 4-5 (e₁).
    // Normal usually doesn't depend on p, so block(0,0) is zero.
    assert(bool(jac_n.leftCols<2>().isZero()));
    const auto J_e0 = jac_n.middleCols<2>(2);
    const auto J_e1 = jac_n.rightCols<2>();

    // --- Block Row 0 (p) ---
    // d/dp [Jₙᵀ v] -> Jₙᵀ => Adds Jᵀ to the row
    // (p, e0) += J_e0
    hess.block<2, 2>(0, 2) += J_e0;
    // (p, e1) += J_e1
    hess.block<2, 2>(0, 4) += J_e1;

    // --- Block Row 1 (e0) ---
    // d/de0 [Jₙᵀ * v + n] -> Jₙᵀ * (-I) + (-I)ᵀ * Jₙ

    // (e0, p) += J_e0^T
    hess.block<2, 2>(2, 0) += J_e0.transpose();

    // (e0, e0) -= (J_e0 + J_e0^T)
    hess.block<2, 2>(2, 2) -= (J_e0 + J_e0.transpose());

    // (e0, e1) -= J_e1
    hess.block<2, 2>(2, 4) -= J_e1;

    // --- Block Row 2 (e1) ---
    // d/de1 [Jₙᵀ * v] -> Jₙᵀ (-I)

    // (e1, p) += J_e1^T
    hess.block<2, 2>(4, 0) += J_e1.transpose();

    // (e1, e0) -= J_e1^T
    hess.block<2, 2>(4, 2) -= J_e1.transpose();

    return hess;
}

} // namespace ipc