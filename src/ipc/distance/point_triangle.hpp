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
        PointTriangleDistanceType dtype);

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
        PointTriangleDistanceType dtype);

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
        PointTriangleDistanceType dtype);
} // namespace detail

// --- EigenExpression wrappers ---

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_distance(
    const DerivedP& p,
    const DerivedT0& t0,
    const DerivedT1& t1,
    const DerivedT2& t2,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance<T>(p, t0, t1, t2, dtype);
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_distance_gradient(
    const DerivedP& p,
    const DerivedT0& t0,
    const DerivedT1& t1,
    const DerivedT2& t2,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance_gradient<T>(p, t0, t1, t2, dtype);
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_distance_hessian(
    const DerivedP& p,
    const DerivedT0& t0,
    const DerivedT1& t1,
    const DerivedT2& t2,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance_hessian<T>(p, t0, t1, t2, dtype);
}

} // namespace ipc
