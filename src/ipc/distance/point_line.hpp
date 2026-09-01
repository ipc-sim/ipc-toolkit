#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    template <typename T>
    void point_line_distance_gradient_2D(
        T v01, T v02, T v11, T v12, T v21, T v22, T g[6]);
    template <typename T>
    void point_line_distance_gradient_3D(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T g[9]);
    template <typename T>
    void point_line_distance_hessian_2D(
        T v01, T v02, T v11, T v12, T v21, T v22, T H[36]);
    template <typename T>
    void point_line_distance_hessian_3D(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T H[81]);
} // namespace autogen

namespace detail {
    /// @brief Compute the distance between a point and line in 2D or 3D.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge defining the line.
    /// @param e1 The second vertex of the edge defining the line.
    /// @return The distance between the point and line.
    template <typename T, int dim>
    inline T point_line_distance(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "point-line is only 2D or 3D");
        if constexpr (dim == 2) {
            const Eigen::Vector2<T> e = e1 - e0;
            const T numerator =
                (e[1] * p[0] - e[0] * p[1] + e1[0] * e0[1] - e1[1] * e0[0]);
            return numerator * numerator / e.squaredNorm();
        } else {
            const Eigen::Vector3<T> p_to_e0 = e0 - p;
            const Eigen::Vector3<T> p_to_e1 = e1 - p;
            return p_to_e0.cross(p_to_e1).squaredNorm()
                / (e1 - e0).squaredNorm();
        }
    }

    /// @brief Compute the gradient of the distance between a point and line.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge defining the line.
    /// @param e1 The second vertex of the edge defining the line.
    /// @return The gradient of the distance wrt p, e0, and e1.
    template <typename T, int dim>
    inline Eigen::Vector<T, 3 * dim> point_line_distance_gradient(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "point-line is only 2D or 3D");
        Eigen::Vector<T, 3 * dim> grad;
        if constexpr (dim == 2) {
            autogen::point_line_distance_gradient_2D(
                p[0], p[1], e0[0], e0[1], e1[0], e1[1], grad.data());
        } else {
            autogen::point_line_distance_gradient_3D(
                p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
                grad.data());
        }
        return grad;
    }

    /// @brief Compute the hessian of the distance between a point and line.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge defining the line.
    /// @param e1 The second vertex of the edge defining the line.
    /// @return The hessian of the distance wrt p, e0, and e1.
    template <typename T, int dim>
    inline Eigen::Matrix<T, 3 * dim, 3 * dim> point_line_distance_hessian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "point-line is only 2D or 3D");
        Eigen::Matrix<T, 3 * dim, 3 * dim> hess;
        if constexpr (dim == 2) {
            autogen::point_line_distance_hessian_2D(
                p[0], p[1], e0[0], e0[1], e1[0], e1[1], hess.data());
        } else {
            autogen::point_line_distance_hessian_3D(
                p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
                hess.data());
        }
        return hess;
    }
} // namespace detail

/// @brief Compute the distance between a point and line in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The distance between the point and line.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto
point_line_distance(const DerivedP& p, const DerivedE0& e0, const DerivedE1& e1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedE0, DerivedE1);
    using T = typename DerivedP::Scalar;

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_line_distance<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_line_distance<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        return detail::point_line_distance<T, 2>(p, e0, e1);
    } else {
        assert(p.size() == 3);
        return detail::point_line_distance<T, 3>(p, e0, e1);
    }
}

/// @brief Compute the gradient of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The gradient of the distance wrt p, e0, and e1.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_distance_gradient(
    const DerivedP& p, const DerivedE0& e0, const DerivedE1& e1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedE0, DerivedE1);
    using T = typename DerivedP::Scalar;
    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_line_distance_gradient<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_line_distance_gradient<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        VectorMax9<T> grad(6);
        autogen::point_line_distance_gradient_2D(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], grad.data());
        return grad;
    } else {
        assert(p.size() == 3);
        VectorMax9<T> grad(9);
        autogen::point_line_distance_gradient_3D(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            grad.data());
        return grad;
    }
}

/// @brief Compute the hessian of the distance between a point and line.
///
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The hessian of the distance wrt p, e0, and e1.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_distance_hessian(
    const DerivedP& p, const DerivedE0& e0, const DerivedE1& e1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedE0, DerivedE1);
    using T = typename DerivedP::Scalar;
    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_line_distance_hessian<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_line_distance_hessian<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        MatrixMax9<T> hess(6, 6);
        autogen::point_line_distance_hessian_2D(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], hess.data());
        return hess;
    } else {
        assert(p.size() == 3);
        MatrixMax9<T> hess(9, 9);
        autogen::point_line_distance_hessian_3D(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            hess.data());
        return hess;
    }
}

} // namespace ipc
