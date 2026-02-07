#include "point_plane.hpp"

namespace ipc {

Vector12d point_plane_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    const Eigen::Vector3d n = triangle_normal(t0, t1, t2);
    const Eigen::Matrix<double, 3, 9> jac_n =
        triangle_normal_jacobian(t0, t1, t2);

    Vector12d grad;
    grad.segment<3>(0) = n;
    grad.segment<3>(3) = jac_n.leftCols<3>().transpose() * (p - t0) - n;
    grad.segment<3>(6) = jac_n.middleCols<3>(3).transpose() * (p - t0);
    grad.segment<3>(9) = jac_n.rightCols<3>().transpose() * (p - t0);

    return grad;
}

Matrix12d point_plane_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    // Precompute normal's Jacobian and Hessian
    const Eigen::Matrix<double, 3, 9> jac_n =
        triangle_normal_jacobian(t0, t1, t2);
    const Eigen::Matrix<double, 27, 9> hess_n =
        triangle_normal_hessian(t0, t1, t2);

    // Vector from t0 to p
    const Eigen::Vector3d v = p - t0;

    Matrix12d hess;

    // ---------------------------------------------------------
    // 0. Fill Point-Point Block (p, p)
    // ---------------------------------------------------------
    // The second derivative w.r.t p is zero since the normal is constant w.r.t
    // p.
    hess.block<3, 3>(0, 0).setZero();

    // ---------------------------------------------------------
    // 1. Fill Mixed Derivatives (p, t) and (t, p)
    // ---------------------------------------------------------
    // The gradient w.r.t p is n.
    // The mixed derivative is the Jacobian of n w.r.t t.
    hess.block<3, 9>(0, 3) = jac_n;
    hess.block<9, 3>(3, 0) = jac_n.transpose();

    // ---------------------------------------------------------
    // 2. Fill Triangle-Triangle Block (t, t)
    // ---------------------------------------------------------
    // Formula: Hₜₜ = (v ⋅ Hₙ) - (δⱼ₀ Jᵢᵀ + δᵢ₀ Jⱼ)

    // A. Contraction of the normal Hessian tensor with vector v
    // hess_n is 3x81. v is 3x1. Result is 1x81, which maps to 9x9.
    for (int i = 0; i < 9; ++i) {
        hess.block<1, 9>(i + 3, 3) =
            hess_n.middleRows<3>(3 * i).transpose() * v;
    }

    // B. Subtract first derivative terms (Product Rule corrections)
    // Extract 3x3 Jacobian blocks for t0, t1, t2
    const auto J0 = jac_n.leftCols<3>();
    const auto J1 = jac_n.middleCols<3>(3);
    const auto J2 = jac_n.rightCols<3>();

    // Apply corrections for terms involving t0 (index 0)

    // Block (t0, t0): i=0, j=0. Subtract J0 + J0^T
    hess.block<3, 3>(3, 3) -= (J0 + J0.transpose());

    // Block (t0, t1): i=0, j=1. Subtract J1
    hess.block<3, 3>(3, 6) -= J1;

    // Block (t0, t2): i=0, j=2. Subtract J2
    hess.block<3, 3>(3, 9) -= J2;

    // Block (t1, t0): i=1, j=0. Subtract J1^T
    hess.block<3, 3>(6, 3) -= J1.transpose();

    // Block (t2, t0): i=2, j=0. Subtract J2^T
    hess.block<3, 3>(9, 3) -= J2.transpose();

    return hess;
}

} // namespace ipc