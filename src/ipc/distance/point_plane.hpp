#pragma once

#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    // clang-format off
    template <typename T>
    void point_plane_distance_gradient(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T g[12]);
    template <typename T>
    void point_plane_distance_hessian(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T H[144]);
    // clang-format on
} // namespace autogen

/// @brief Compute the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The distance between the point and plane.
template <typename T>
inline T point_plane_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> origin,
    Eigen::ConstRef<Eigen::Vector3<T>> normal)
{
    const T point_to_plane = (p - origin).dot(normal);
    return point_to_plane * point_to_plane / normal.squaredNorm();
}

/// @brief Compute the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance between the point and plane.
template <typename T>
inline T point_plane_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    // Inline the (p, origin, normal) overload to avoid two-phase lookup issues
    // with dependent argument types.
    const Eigen::Vector3<T> n = triangle_normal(t0, t1, t2);
    const T d = (p - t0).dot(n);
    return d * d / n.squaredNorm();
}

/// @brief Compute the gradient of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The gradient of the distance wrt p.
template <typename T>
inline Eigen::Vector3<T> point_plane_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> origin,
    Eigen::ConstRef<Eigen::Vector3<T>> normal)
{
    return (T(2) / normal.squaredNorm()) * (p - origin).dot(normal) * normal;
}

/// @brief Compute the gradient of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The gradient of the distance wrt p, t0, t1, and t2.
template <typename T>
inline Eigen::Vector<T, 12> point_plane_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    Eigen::Vector<T, 12> grad;
    autogen::point_plane_distance_gradient(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], grad.data());
    return grad;
}

/// @brief Compute the hessian of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The hessian of the distance wrt p.
template <typename T>
inline Eigen::Matrix3<T> point_plane_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> origin,
    Eigen::ConstRef<Eigen::Vector3<T>> normal)
{
    return T(2) / normal.squaredNorm() * normal * normal.transpose();
}

/// @brief Compute the hessian of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The hessian of the distance wrt p, t0, t1, and t2.
template <typename T>
inline Eigen::Matrix<T, 12, 12> point_plane_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    Eigen::Matrix<T, 12, 12> hess;
    autogen::point_plane_distance_hessian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], hess.data());
    return hess;
}

// --- EigenExpression wrappers ---

template <
    EigenExpression DerivedP,
    EigenExpression DerivedOrigin,
    EigenExpression DerivedNormal>
inline auto point_plane_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedOrigin>& origin,
    const Eigen::MatrixBase<DerivedNormal>& normal) -> typename DerivedP::Scalar
{
    using T = typename DerivedP::Scalar;
    return point_plane_distance(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(origin),
        Eigen::Ref<const Eigen::Vector3<T>>(normal));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedT0,
    EigenExpression DerivedT1,
    EigenExpression DerivedT2>
inline auto point_plane_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2) -> typename DerivedP::Scalar
{
    using T = typename DerivedP::Scalar;
    return point_plane_distance(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedOrigin,
    EigenExpression DerivedNormal>
inline auto point_plane_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedOrigin>& origin,
    const Eigen::MatrixBase<DerivedNormal>& normal)
    -> Eigen::Vector3<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_plane_distance_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(origin),
        Eigen::Ref<const Eigen::Vector3<T>>(normal));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedT0,
    EigenExpression DerivedT1,
    EigenExpression DerivedT2>
inline auto point_plane_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
    -> Eigen::Vector<typename DerivedP::Scalar, 12>
{
    using T = typename DerivedP::Scalar;
    return point_plane_distance_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedOrigin,
    EigenExpression DerivedNormal>
inline auto point_plane_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedOrigin>& origin,
    const Eigen::MatrixBase<DerivedNormal>& normal)
    -> Eigen::Matrix<typename DerivedP::Scalar, 3, 3>
{
    using T = typename DerivedP::Scalar;
    return point_plane_distance_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(origin),
        Eigen::Ref<const Eigen::Vector3<T>>(normal));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedT0,
    EigenExpression DerivedT1,
    EigenExpression DerivedT2>
inline auto point_plane_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
    -> Eigen::Matrix<typename DerivedP::Scalar, 12, 12>
{
    using T = typename DerivedP::Scalar;
    return point_plane_distance_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

} // namespace ipc
