#pragma once

#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// Compute the signed distance from a point to the plane of a triangle.
///
/// The signed distance is computed as:
///   d = triangle_normal(t0, t1, t2) · (p - t0)
/// The sign corresponds to the orientation of the triangle normal.
///
/// @param p   The query point (3D).
/// @param t0  First vertex of the triangle (3D), used as a point on the plane.
/// @param t1  Second vertex of the triangle (3D).
/// @param t2  Third vertex of the triangle (3D).
/// @return    The signed distance from p to the plane of the triangle.
template <typename T>
inline T point_plane_signed_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    return triangle_normal(t0, t1, t2).dot(p - t0);
}

/// Compute the gradient of the signed point-to-plane distance.
///
/// The returned vector contains derivatives in the order:
///   [ ∂d/∂p (3), ∂d/∂t0 (3), ∂d/∂t1 (3), ∂d/∂t2 (3) ]
///
/// @param p   The query point (3D).
/// @param t0  First vertex of the triangle (3D).
/// @param t1  Second vertex of the triangle (3D).
/// @param t2  Third vertex of the triangle (3D).
/// @return    A 12-vector containing the gradient of the signed distance.
template <typename T>
inline Eigen::Vector<T, 12> point_plane_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    const Eigen::Vector3<T> n = triangle_normal(t0, t1, t2);
    const Eigen::Matrix<T, 3, 9> jac_n = triangle_normal_jacobian(t0, t1, t2);
    Eigen::Vector<T, 12> grad;
    grad.template segment<3>(0) = n;
    grad.template segment<3>(3) =
        jac_n.template leftCols<3>().transpose() * (p - t0) - n;
    grad.template segment<3>(6) =
        jac_n.template middleCols<3>(3).transpose() * (p - t0);
    grad.template segment<3>(9) =
        jac_n.template rightCols<3>().transpose() * (p - t0);
    return grad;
}

/// Compute the Hessian (matrix of second derivatives) of the signed distance.
///
/// The returned matrix is 12x12 with variables ordered as:
///   [ p (3), t0 (3), t1 (3), t2 (3) ].
///
/// @param p   The query point (3D).
/// @param t0  First vertex of the triangle (3D).
/// @param t1  Second vertex of the triangle (3D).
/// @param t2  Third vertex of the triangle (3D).
/// @return    A 12x12 Hessian matrix of the signed distance.
template <typename T>
Eigen::Matrix<T, 12, 12> point_plane_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2);

// --- EigenExpression wrappers ---

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_plane_signed_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2) -> typename DerivedP::Scalar
{
    using T = typename DerivedP::Scalar;
    return point_plane_signed_distance(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_plane_signed_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
    -> Eigen::Vector<typename DerivedP::Scalar, 12>
{
    using T = typename DerivedP::Scalar;
    return point_plane_signed_distance_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_plane_signed_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
    -> Eigen::Matrix<typename DerivedP::Scalar, 12, 12>
{
    using T = typename DerivedP::Scalar;
    return point_plane_signed_distance_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

} // namespace ipc
