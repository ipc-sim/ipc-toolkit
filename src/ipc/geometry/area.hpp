#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Geometry>

#include <cassert>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {

    // clang-format off
    /// dA is (9×1) flattened in column-major order
    template <typename T>
    void triangle_area_gradient(
        T t0_x, T t0_y, T t0_z, T t1_x, T t1_y, T t1_z, T t2_x, T t2_y, T t2_z, T dA[9]);
    // clang-format on

} // namespace autogen

namespace detail {
    /// @brief Compute the length of an edge.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param e0 The first vertex of the edge.
    /// @param e1 The second vertex of the edge.
    /// @return The length of the edge.
    template <typename T, int dim>
    inline T edge_length(
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "edges are only 2D or 3D");
        return (e1 - e0).norm();
    }

    /// @brief Compute the gradient of an edge's length.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param e0 The first vertex of the edge.
    /// @param e1 The second vertex of the edge.
    /// @return The gradient of the edge's length wrt e0, and e1.
    template <typename T, int dim>
    inline Eigen::Vector<T, 2 * dim> edge_length_gradient(
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "edges are only 2D or 3D");
        assert((e1 - e0).norm() != 0);

        // ∇ ‖e₁ - e₀‖
        Eigen::Vector<T, 2 * dim> grad;
        grad.template head<dim>() = (e0 - e1) / (e1 - e0).norm();
        grad.template tail<dim>() = -grad.template head<dim>();
        return grad;
    }

    /// @brief Compute the area of a triangle.
    /// @tparam T The scalar type.
    /// @param t0 The first vertex of the triangle.
    /// @param t1 The second vertex of the triangle.
    /// @param t2 The third vertex of the triangle.
    /// @return The area of the triangle.
    template <typename T>
    inline T triangle_area(
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2)
    {
        return T(0.5) * (t1 - t0).cross(t2 - t0).norm();
    }

    /// @brief Compute the gradient of the area of a triangle.
    /// @tparam T The scalar type.
    /// @param t0 The first vertex of the triangle.
    /// @param t1 The second vertex of the triangle.
    /// @param t2 The third vertex of the triangle.
    /// @return The gradient of the triangle's area t0, t1, and t2.
    template <typename T>
    inline Eigen::Vector<T, 9> triangle_area_gradient(
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2)
    {
        Eigen::Vector<T, 9> grad;
        autogen::triangle_area_gradient(
            t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0], t2[1], t2[2],
            grad.data());
        return grad;
    }
} // namespace detail

/// @brief Compute the length of an edge.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The length of the edge.
template <typename DerivedE0, typename DerivedE1>
inline auto edge_length(const DerivedE0& e0, const DerivedE1& e1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedE0, DerivedE1);
    using T = typename DerivedE0::Scalar;

    if constexpr (dim_v<DerivedE0> == 2) {
        return detail::edge_length<T, 2>(e0, e1);
    } else if constexpr (dim_v<DerivedE0> == 3) {
        return detail::edge_length<T, 3>(e0, e1);
    } else if (e0.size() == 2) {
        return detail::edge_length<T, 2>(e0, e1);
    } else {
        assert(e0.size() == 3);
        return detail::edge_length<T, 3>(e0, e1);
    }
}

/// @brief Compute the gradient of an edge's length.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The gradient of the edge's length wrt e0, and e1.
template <typename DerivedE0, typename DerivedE1>
inline auto edge_length_gradient(const DerivedE0& e0, const DerivedE1& e1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedE0, DerivedE1);
    using T = typename DerivedE0::Scalar;

    if constexpr (dim_v<DerivedE0> == 2) {
        return detail::edge_length_gradient<T, 2>(e0, e1);
    } else if constexpr (dim_v<DerivedE0> == 3) {
        return detail::edge_length_gradient<T, 3>(e0, e1);
    } else if (e0.size() == 2) {
        return VectorMax6<T>(detail::edge_length_gradient<T, 2>(e0, e1));
    } else {
        assert(e0.size() == 3);
        return VectorMax6<T>(detail::edge_length_gradient<T, 3>(e0, e1));
    }
}

/// @brief Compute the area of a triangle.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The area of the triangle.
template <typename DerivedT0, typename DerivedT1, typename DerivedT2>
inline auto
triangle_area(const DerivedT0& t0, const DerivedT1& t1, const DerivedT2& t2)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedT0, DerivedT1, DerivedT2);
    // NOTE: explicit <T>: the detail overload cannot deduce T from an
    // arbitrary Eigen expression.
    using T = typename DerivedT0::Scalar;
    return detail::triangle_area<T>(t0, t1, t2);
}

/// @brief Compute the gradient of the area of a triangle.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The gradient of the triangle's area t0, t1, and t2.
template <typename DerivedT0, typename DerivedT1, typename DerivedT2>
inline auto triangle_area_gradient(
    const DerivedT0& t0, const DerivedT1& t1, const DerivedT2& t2)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedT0, DerivedT1, DerivedT2);
    using T = typename DerivedT0::Scalar;
    return detail::triangle_area_gradient<T>(t0, t1, t2);
}

} // namespace ipc
