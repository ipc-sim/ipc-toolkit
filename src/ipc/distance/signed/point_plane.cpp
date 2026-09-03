#include "point_plane.hpp"

#include <ipc/utils/simd.hpp>

namespace ipc::detail {

template <typename T>
Eigen::Matrix<T, 12, 12> point_plane_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    // Precompute normal's Jacobian and Hessian
    const Eigen::Matrix<T, 3, 9> jac_n = triangle_normal_jacobian(t0, t1, t2);
    const Eigen::Matrix<T, 27, 9> hess_n = triangle_normal_hessian(t0, t1, t2);

    // Vector from t0 to p
    const Eigen::Vector3<T> v = p - t0;

    Eigen::Matrix<T, 12, 12> hess;

    // ---------------------------------------------------------
    // 0. Fill Point-Point Block (p, p)
    // ---------------------------------------------------------
    // The second derivative w.r.t p is zero since the normal is constant w.r.t
    // p.
    hess.template block<3, 3>(0, 0).setZero();

    // ---------------------------------------------------------
    // 1. Fill Mixed Derivatives (p, t) and (t, p)
    // ---------------------------------------------------------
    // The gradient w.r.t p is n.
    // The mixed derivative is the Jacobian of n w.r.t t.
    hess.template block<3, 9>(0, 3) = jac_n;
    hess.template block<9, 3>(3, 0) = jac_n.transpose();

    // ---------------------------------------------------------
    // 2. Fill Triangle-Triangle Block (t, t)
    // ---------------------------------------------------------
    // Formula: Hₜₜ = (v ⋅ Hₙ) - (δⱼ₀ Jᵢᵀ + δᵢ₀ Jⱼ)

    // A. Contraction of the normal Hessian tensor with vector v
    // hess_n is 3x81. v is 3x1. Result is 1x81, which maps to 9x9.
    hess.template block<9, 9>(3, 3) =
        (hess_n.reshaped(3, 81).transpose() * v).reshaped(9, 9);

    // B. Subtract first derivative terms (Product Rule corrections)
    // Extract 3x3 Jacobian blocks for t0, t1, t2
    const auto J0 = jac_n.template leftCols<3>();
    const auto J1 = jac_n.template middleCols<3>(3);
    const auto J2 = jac_n.template rightCols<3>();

    // Apply corrections for terms involving t0 (index 0)

    // Block (t0, t0): i=0, j=0. Subtract J0 + J0^T
    hess.template block<3, 3>(3, 3) -= (J0 + J0.transpose());

    // Block (t0, t1): i=0, j=1. Subtract J1
    hess.template block<3, 3>(3, 6) -= J1;

    // Block (t0, t2): i=0, j=2. Subtract J2
    hess.template block<3, 3>(3, 9) -= J2;

    // Block (t1, t0): i=1, j=0. Subtract J1^T
    hess.template block<3, 3>(6, 3) -= J1.transpose();

    // Block (t2, t0): i=2, j=0. Subtract J2^T
    hess.template block<3, 3>(9, 3) -= J2.transpose();

    return hess;
}

#define IPC_INSTANTIATE_POINT_PLANE_SIGNED_DISTANCE_HESSIAN(T)                 \
    template Eigen::Matrix<T, 12, 12> point_plane_signed_distance_hessian<T>(  \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>)

IPC_INSTANTIATE_POINT_PLANE_SIGNED_DISTANCE_HESSIAN(float);
IPC_INSTANTIATE_POINT_PLANE_SIGNED_DISTANCE_HESSIAN(double);
#ifdef IPC_TOOLKIT_WITH_SIMD
IPC_INSTANTIATE_POINT_PLANE_SIGNED_DISTANCE_HESSIAN(SimdBatch<float>);
IPC_INSTANTIATE_POINT_PLANE_SIGNED_DISTANCE_HESSIAN(SimdBatch<double>);
#endif

#undef IPC_INSTANTIATE_POINT_PLANE_SIGNED_DISTANCE_HESSIAN

} // namespace ipc::detail
