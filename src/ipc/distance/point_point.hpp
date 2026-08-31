#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <cassert>

namespace ipc {

namespace detail {
    /// @brief Compute the distance between two points.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p0 The first point.
    /// @param p1 The second point.
    /// @return The distance between p0 and p1.
    template <typename T, int dim>
    inline T point_point_distance(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> p1)
    {
        static_assert(dim == 2 || dim == 3, "point-point is only 2D or 3D");
        return (p1 - p0).squaredNorm();
    }

    /// @brief Compute the gradient of the distance between two points.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p0 The first point.
    /// @param p1 The second point.
    /// @return The computed gradient.
    template <typename T, int dim>
    inline Eigen::Vector<T, 2 * dim> point_point_distance_gradient(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> p1)
    {
        static_assert(dim == 2 || dim == 3, "point-point is only 2D or 3D");
        Eigen::Vector<T, 2 * dim> grad;
        grad.template head<dim>() = T(2) * (p0 - p1);
        grad.template tail<dim>() = -grad.template head<dim>();
        return grad;
    }

    /// @brief Compute the hessian of the distance between two points.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p0 The first point.
    /// @param p1 The second point.
    /// @return The computed hessian.
    template <typename T, int dim>
    inline Eigen::Matrix<T, 2 * dim, 2 * dim> point_point_distance_hessian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> /*p0*/,
        Eigen::ConstRef<Eigen::Vector<T, dim>> /*p1*/)
    {
        static_assert(dim == 2 || dim == 3, "point-point is only 2D or 3D");
        const Eigen::Matrix<T, dim, dim> I2 =
            T(2) * Eigen::Matrix<T, dim, dim>::Identity();
        Eigen::Matrix<T, 2 * dim, 2 * dim> hess;
        hess.template topLeftCorner<dim, dim>() = I2;
        hess.template bottomRightCorner<dim, dim>() = I2;
        hess.template topRightCorner<dim, dim>() = -I2;
        hess.template bottomLeftCorner<dim, dim>() = -I2;
        return hess;
    }
} // namespace detail

/// @brief Compute the distance between two points.
///
/// @note The distance is actually squared distance.
///
/// Accepts any Eigen expression. The dimension is resolved at compile time when
/// the argument type knows it and with a single branch otherwise; each argument
/// is evaluated exactly once.
///
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The distance between p0 and p1.
template <typename DerivedP0, typename DerivedP1>
inline auto point_point_distance(const DerivedP0& p0, const DerivedP1& p1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP0, DerivedP1);
    using T = typename DerivedP0::Scalar;

    if constexpr (dim_v<DerivedP0> == 2) {
        return detail::point_point_distance<T, 2>(p0, p1);
    } else if constexpr (dim_v<DerivedP0> == 3) {
        return detail::point_point_distance<T, 3>(p0, p1);
    } else if (p0.size() == 2) {
        return detail::point_point_distance<T, 2>(p0, p1);
    } else {
        assert(p0.size() == 3);
        return detail::point_point_distance<T, 3>(p0, p1);
    }
}

/// @brief Compute the gradient of the distance between two points.
///
/// @note The distance is actually squared distance.
///
/// Accepts any Eigen expression. When the dimension is known at compile time
/// the return type is the fixed-size Eigen::Vector<T, 2 * dim>; otherwise it is
/// VectorMax6<T>.
///
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The computed gradient.
template <typename DerivedP0, typename DerivedP1>
inline auto
point_point_distance_gradient(const DerivedP0& p0, const DerivedP1& p1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP0, DerivedP1);
    using T = typename DerivedP0::Scalar;
    if constexpr (dim_v<DerivedP0> == 2) {
        return detail::point_point_distance_gradient<T, 2>(p0, p1);
    } else if constexpr (dim_v<DerivedP0> == 3) {
        return detail::point_point_distance_gradient<T, 3>(p0, p1);
    } else if (p0.size() == 2) {
        return VectorMax6<T>(
            detail::point_point_distance_gradient<T, 2>(p0, p1));
    } else {
        assert(p0.size() == 3);
        return VectorMax6<T>(
            detail::point_point_distance_gradient<T, 3>(p0, p1));
    }
}

/// @brief Compute the hessian of the distance between two points.
///
/// @note The distance is actually squared distance.
///
/// Accepts any Eigen expression. When the dimension is known at compile time
/// the return type is the fixed-size Eigen::Matrix<T, 2 * dim, 2 * dim>;
/// otherwise it is MatrixMax6<T>.
///
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The computed hessian.
template <typename DerivedP0, typename DerivedP1>
inline auto
point_point_distance_hessian(const DerivedP0& p0, const DerivedP1& p1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP0, DerivedP1);
    using T = typename DerivedP0::Scalar;
    if constexpr (dim_v<DerivedP0> == 2) {
        return detail::point_point_distance_hessian<T, 2>(p0, p1);
    } else if constexpr (dim_v<DerivedP0> == 3) {
        return detail::point_point_distance_hessian<T, 3>(p0, p1);
    } else if (p0.size() == 2) {
        return MatrixMax6<T>(
            detail::point_point_distance_hessian<T, 2>(p0, p1));
    } else {
        assert(p0.size() == 3);
        return MatrixMax6<T>(
            detail::point_point_distance_hessian<T, 3>(p0, p1));
    }
}

} // namespace ipc
