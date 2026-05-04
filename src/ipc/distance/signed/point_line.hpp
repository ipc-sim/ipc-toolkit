#pragma once

#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// Compute the signed distance from a point to a directed line segment.
///
/// The signed distance is computed as d = n · (p - e0),
/// where n = point_line_normal(p, e0, e1) is the unit normal associated with
/// the directed edge (e0 -> e1) chosen consistently for the point p.
/// Positive/negative sign indicates which side of the directed edge the point
/// lies on.
///
/// @param p   The query point (2D).
/// @param e0  The first endpoint of the directed edge (2D).
/// @param e1  The second endpoint of the directed edge (2D).
/// @return    The signed scalar distance from p to the infinite line through e0 and e1.
/// @note      The edge must be non-degenerate (e0 != e1).
template <typename T>
inline T point_line_signed_distance(
    Eigen::ConstRef<Eigen::Vector2<T>> p,
    Eigen::ConstRef<Eigen::Vector2<T>> e0,
    Eigen::ConstRef<Eigen::Vector2<T>> e1)
{
    return point_line_normal(p, e0, e1).dot(p - e0);
}

/// Compute the gradient of the signed point-to-line distance with respect to
/// all input coordinates packed as [p_x, p_y, e0_x, e0_y, e1_x, e1_y]^T.
///
/// The returned Vector6d contains partial derivatives of the scalar signed
/// distance d with respect to each coordinate in the above ordering:
///   grad = [∂d/∂p, ∂d/∂e0, ∂d/∂e1]^T.
///
/// @param p   The query point (2D).
/// @param e0  The first endpoint of the directed edge (2D).
/// @param e1  The second endpoint of the directed edge (2D).
/// @return    A 6-vector containing the gradient of the signed distance.
/// @note      The edge must be non-degenerate (e0 != e1).
template <typename T>
inline Eigen::Vector<T, 6> point_line_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector2<T>> p,
    Eigen::ConstRef<Eigen::Vector2<T>> e0,
    Eigen::ConstRef<Eigen::Vector2<T>> e1)
{
    const Eigen::Vector2<T> n = point_line_normal(p, e0, e1);
    const Eigen::Matrix<T, 2, 6> jac_n = point_line_normal_jacobian(p, e0, e1);
    Eigen::Vector<T, 6> grad = jac_n.transpose() * (p - e0);
    grad.template segment<2>(0) += n;
    grad.template segment<2>(2) -= n;
    return grad;
}

/// Compute the Hessian (second derivatives) of the signed point-to-line
/// distance with respect to all input coordinates packed as [p_x, p_y, e0_x,
/// e0_y, e1_x, e1_y]^T.
///
/// The returned Matrix6d is the symmetric 6x6 matrix of second partial
/// derivatives:
///   H_ij = ∂^2 d / (∂x_i ∂x_j),
/// with the same coordinate ordering as in the gradient.
///
/// @param p   The query point (2D).
/// @param e0  The first endpoint of the directed edge (2D).
/// @param e1  The second endpoint of the directed edge (2D).
/// @return    A 6x6 Hessian matrix of the signed distance.
/// @note      The edge must be non-degenerate (e0 != e1).
template <typename T>
Eigen::Matrix<T, 6, 6> point_line_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector2<T>> p,
    Eigen::ConstRef<Eigen::Vector2<T>> e0,
    Eigen::ConstRef<Eigen::Vector2<T>> e1);

// --- EigenExpression wrappers ---

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_signed_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1) -> typename DerivedP::Scalar
{
    using T = typename DerivedP::Scalar;
    return point_line_signed_distance(
        Eigen::Ref<const Eigen::Vector2<T>>(p),
        Eigen::Ref<const Eigen::Vector2<T>>(e0),
        Eigen::Ref<const Eigen::Vector2<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_signed_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> Eigen::Vector<typename DerivedP::Scalar, 6>
{
    using T = typename DerivedP::Scalar;
    return point_line_signed_distance_gradient(
        Eigen::Ref<const Eigen::Vector2<T>>(p),
        Eigen::Ref<const Eigen::Vector2<T>>(e0),
        Eigen::Ref<const Eigen::Vector2<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_signed_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> Eigen::Matrix<typename DerivedP::Scalar, 6, 6>
{
    using T = typename DerivedP::Scalar;
    return point_line_signed_distance_hessian(
        Eigen::Ref<const Eigen::Vector2<T>>(p),
        Eigen::Ref<const Eigen::Vector2<T>>(e0),
        Eigen::Ref<const Eigen::Vector2<T>>(e1));
}

} // namespace ipc
