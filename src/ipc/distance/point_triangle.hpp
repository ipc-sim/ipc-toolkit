#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

namespace detail {
    /// @brief Compute the distance between a points and a triangle.
    /// @note The distance is actually squared distance.
    /// @param p The point.
    /// @param t0 The first vertex of the triangle.
    /// @param t1 The second vertex of the triangle.
    /// @param t2 The third vertex of the triangle.
    /// @param dtype The point-triangle distance type to compute.
    /// @return The distance between the point and triangle.
    template <typename T>
    T point_triangle_distance(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2,
        PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

    /// @brief Compute the gradient of the distance between a points and a triangle.
    /// @note The distance is actually squared distance.
    /// @param p The point.
    /// @param t0 The first vertex of the triangle.
    /// @param t1 The second vertex of the triangle.
    /// @param t2 The third vertex of the triangle.
    /// @param dtype The point-triangle distance type to compute.
    /// @return The gradient of the distance wrt p, t0, t1, and t2.
    template <typename T>
    Eigen::Vector<T, 12> point_triangle_distance_gradient(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2,
        PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

    /// @brief Compute the hessian of the distance between a points and a triangle.
    /// @note The distance is actually squared distance.
    /// @param p The point.
    /// @param t0 The first vertex of the triangle.
    /// @param t1 The second vertex of the triangle.
    /// @param t2 The third vertex of the triangle.
    /// @param dtype The point-triangle distance type to compute.
    /// @return The hessian of the distance wrt p, t0, t1, and t2.
    template <typename T>
    Eigen::Matrix<T, 12, 12> point_triangle_distance_hessian(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2,
        PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);
} // namespace detail

// --- EigenExpression wrappers ---

/// @brief Compute the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @param dtype The point-triangle distance type to compute.
/// @return The distance between the point and triangle.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance<T>(p, t0, t1, t2, dtype);
}

/// @brief Compute the gradient of the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @param dtype The point-triangle distance type to compute.
/// @return The gradient of the distance wrt p, t0, t1, and t2.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance_gradient<T>(p, t0, t1, t2, dtype);
}

/// @brief Compute the hessian of the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @param dtype The point-triangle distance type to compute.
/// @return The hessian of the distance wrt p, t0, t1, and t2.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance_hessian<T>(p, t0, t1, t2, dtype);
}

} // namespace ipc
