#pragma once

#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

namespace detail {
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
        const Eigen::Matrix<T, 3, 9> jac_n =
            triangle_normal_jacobian(t0, t1, t2);
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

    /// Compute the Hessian (matrix of second derivatives) of the signed
    /// distance.
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
} // namespace detail

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_plane_signed_distance(
    const DerivedP& p,
    const DerivedT0& t0,
    const DerivedT1& t1,
    const DerivedT2& t2)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedP::Scalar;
    // NOTE: explicit <T>: the detail overload cannot deduce T from an
    // arbitrary Eigen expression.
    return detail::point_plane_signed_distance<T>(p, t0, t1, t2);
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_plane_signed_distance_gradient(
    const DerivedP& p,
    const DerivedT0& t0,
    const DerivedT1& t1,
    const DerivedT2& t2)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedP::Scalar;
    return detail::point_plane_signed_distance_gradient<T>(p, t0, t1, t2);
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_plane_signed_distance_hessian(
    const DerivedP& p,
    const DerivedT0& t0,
    const DerivedT1& t1,
    const DerivedT2& t2)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedP::Scalar;
    return detail::point_plane_signed_distance_hessian<T>(p, t0, t1, t2);
}

} // namespace ipc
