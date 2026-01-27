#include "line_line.hpp"

namespace ipc {

Vector12d line_line_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d n = line_line_normal(ea0, ea1, eb0, eb1);
    const Eigen::Matrix<double, 3, 12> jac_n =
        line_line_normal_jacobian(ea0, ea1, eb0, eb1);

    Vector12d grad = jac_n.transpose() * (ea0 - eb0);
    grad.segment<3>(0) += n;
    grad.segment<3>(6) -= n;

    return grad;
}

Matrix12d line_line_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    // Precompute normal's Jacobian and Hessian
    const Eigen::Matrix<double, 3, 12> jac_n =
        line_line_normal_jacobian(ea0, ea1, eb0, eb1);
    const Eigen::Matrix<double, 3, 144> hess_n =
        line_line_normal_hessian(ea0, ea1, eb0, eb1);

    // Vector from eb0 to ea0 (Distance vector)
    const Eigen::Vector3d v = ea0 - eb0;

    Matrix12d hess;

    // ---------------------------------------------------------
    // 1. Tensor Contraction (Curvature Term)
    // ---------------------------------------------------------
    // Contract the normal Hessian (3x144) with vector v (3x1).
    // This computes (v ⋅ d²n/dx²).
    // The result is a 1x144 vector, which maps to the 12x12 Hessian matrix.
    hess = (v.transpose() * hess_n).reshaped(12, 12);

    // ---------------------------------------------------------
    // 2. Add Jacobian Terms (Product Rule Corrections)
    // ---------------------------------------------------------
    // The gradient is Jₙᵀ v + Jᵥᵀ n.
    // The Hessian correction involves Jₙᵀ * Jᵥ + Jᵥᵀ * Jₙ.
    // Since v = ea0 - eb0:
    // d(v)/d(ea0) = I, d(v)/d(eb0) = -I, others are 0.

    // Extract 3x3 sub-blocks for cleaner code
    // Indices: 0->ea0, 1->ea1, 2->eb0, 3->eb1
    const auto J_a0 = jac_n.leftCols<3>();
    const auto J_a1 = jac_n.middleCols<3>(3);
    const auto J_b0 = jac_n.middleCols<3>(6);
    const auto J_b1 = jac_n.rightCols<3>();

    // --- Block Row 0 (ea0) ---
    // d/d(ea0) [Jₙᵀ v + n] -> adds J + Jᵀ terms
    hess.block<3, 3>(0, 0) += J_a0 + J_a0.transpose();
    hess.block<3, 3>(0, 3) += J_a1;
    hess.block<3, 3>(0, 6) += J_b0 - J_a0.transpose();
    hess.block<3, 3>(0, 9) += J_b1;

    // --- Block Row 1 (ea1) ---
    // d/d(ea1) [Jₙᵀ v] -> adds Jᵀ terms
    hess.block<3, 3>(3, 0) += J_a1.transpose();
    hess.block<3, 3>(3, 6) -= J_a1.transpose();

    // --- Block Row 2 (eb0) ---
    // d/d(eb0) [Jₙᵀ v - n] -> adds -J - Jᵀ terms
    hess.block<3, 3>(6, 0) += J_b0.transpose() - J_a0;
    hess.block<3, 3>(6, 3) -= J_a1;
    hess.block<3, 3>(6, 6) -= J_b0 + J_b0.transpose();
    hess.block<3, 3>(6, 9) -= J_b1;

    // --- Block Row 3 (eb1) ---
    // d/d(eb1) [Jₙᵀ v] -> adds Jᵀ terms
    hess.block<3, 3>(9, 0) += J_b1.transpose();
    hess.block<3, 3>(9, 6) -= J_b1.transpose();

    return hess;
}

} // namespace ipc