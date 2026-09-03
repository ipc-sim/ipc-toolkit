#include "line_line.hpp"

#include <ipc/utils/simd.hpp>

namespace ipc::detail {

template <typename T>
Eigen::Matrix<T, 12, 12> line_line_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    // Precompute normal's Jacobian and Hessian
    const Eigen::Matrix<T, 3, 12> jac_n =
        line_line_normal_jacobian(ea0, ea1, eb0, eb1);
    const Eigen::Matrix<T, 36, 12> hess_n =
        line_line_normal_hessian(ea0, ea1, eb0, eb1);

    // Vector from eb0 to ea0 (Distance vector)
    const Eigen::Vector3<T> v = ea0 - eb0;

    Eigen::Matrix<T, 12, 12> hess;

    // ---------------------------------------------------------
    // 1. Tensor Contraction (Curvature Term)
    // ---------------------------------------------------------
    // Contract the normal Hessian (3x12x12) with vector v (3x1).
    // This computes (v ⋅ d²n/dx²).
    // The result is a 1x12x12 vector, which maps to the 12x12 Hessian matrix.
    hess = (hess_n.reshaped(3, 144).transpose() * v).reshaped(12, 12);

    // ---------------------------------------------------------
    // 2. Add Jacobian Terms (Product Rule Corrections)
    // ---------------------------------------------------------
    // The gradient is Jₙᵀ v + Jᵥᵀ n.
    // The Hessian correction involves Jₙᵀ * Jᵥ + Jᵥᵀ * Jₙ.
    // Since v = ea0 - eb0:
    // d(v)/d(ea0) = I, d(v)/d(eb0) = -I, others are 0.

    // Extract 3x3 sub-blocks for cleaner code
    // Indices: 0->ea0, 1->ea1, 2->eb0, 3->eb1
    const auto J_a0 = jac_n.template leftCols<3>();
    const auto J_a1 = jac_n.template middleCols<3>(3);
    const auto J_b0 = jac_n.template middleCols<3>(6);
    const auto J_b1 = jac_n.template rightCols<3>();

    // --- Block Row 0 (ea0) ---
    // d/d(ea0) [Jₙᵀ v + n] -> adds J + Jᵀ terms
    hess.template block<3, 3>(0, 0) += J_a0 + J_a0.transpose();
    hess.template block<3, 3>(0, 3) += J_a1;
    hess.template block<3, 3>(0, 6) += J_b0 - J_a0.transpose();
    hess.template block<3, 3>(0, 9) += J_b1;

    // --- Block Row 1 (ea1) ---
    // d/d(ea1) [Jₙᵀ v] -> adds Jᵀ terms
    hess.template block<3, 3>(3, 0) += J_a1.transpose();
    hess.template block<3, 3>(3, 6) -= J_a1.transpose();

    // --- Block Row 2 (eb0) ---
    // d/d(eb0) [Jₙᵀ v - n] -> adds -J - Jᵀ terms
    hess.template block<3, 3>(6, 0) += J_b0.transpose() - J_a0;
    hess.template block<3, 3>(6, 3) -= J_a1;
    hess.template block<3, 3>(6, 6) -= J_b0 + J_b0.transpose();
    hess.template block<3, 3>(6, 9) -= J_b1;

    // --- Block Row 3 (eb1) ---
    // d/d(eb1) [Jₙᵀ v] -> adds Jᵀ terms
    hess.template block<3, 3>(9, 0) += J_b1.transpose();
    hess.template block<3, 3>(9, 6) -= J_b1.transpose();

    return hess;
}

#define IPC_INSTANTIATE_LINE_LINE_SIGNED_DISTANCE_HESSIAN(T)                   \
    template Eigen::Matrix<T, 12, 12> line_line_signed_distance_hessian<T>(    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>)

IPC_INSTANTIATE_LINE_LINE_SIGNED_DISTANCE_HESSIAN(float);
IPC_INSTANTIATE_LINE_LINE_SIGNED_DISTANCE_HESSIAN(double);
#ifdef IPC_TOOLKIT_WITH_SIMD
IPC_INSTANTIATE_LINE_LINE_SIGNED_DISTANCE_HESSIAN(SimdBatch<float>);
IPC_INSTANTIATE_LINE_LINE_SIGNED_DISTANCE_HESSIAN(SimdBatch<double>);
#endif

#undef IPC_INSTANTIATE_LINE_LINE_SIGNED_DISTANCE_HESSIAN

} // namespace ipc::detail
