#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <cassert>

namespace ipc {

namespace detail {
    // ========================================================================
    // Point - Point

    /// @brief Compute the relative velocity of two points
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param dp0 Velocity of the first point
    /// @param dp1 Velocity of the second point
    /// @return The relative velocity of the two points
    template <typename T, int dim>
    inline Eigen::Vector<T, dim> point_point_relative_velocity(
        Eigen::ConstRef<Eigen::Vector<T, dim>> dp0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> dp1)
    {
        static_assert(
            dim == 2 || dim == 3, "point-point velocity is only 2D or 3D");
        return dp0 - dp1;
    }

    /// @brief Compute du/dx where u is the relative velocity of two points
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @return The relative velocity Jacobian du/dx
    template <typename T, int dim>
    inline Eigen::Matrix<T, dim, 2 * dim>
    point_point_relative_velocity_jacobian()
    {
        static_assert(
            dim == 2 || dim == 3, "point-point velocity is only 2D or 3D");
        Eigen::Matrix<T, dim, 2 * dim> J;
        J.template leftCols<dim>() = Eigen::Matrix<T, dim, dim>::Identity();
        J.template rightCols<dim>() = -Eigen::Matrix<T, dim, dim>::Identity();
        return J;
    }

    /// @brief Compute d²u/dxdβ where u is the relative velocity of two points
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @return The vectorized tensor of d²u/dxdβ
    template <typename T, int dim>
    inline Eigen::Vector<T, 2 * dim * dim>
    point_point_relative_velocity_dx_dbeta()
    {
        static_assert(
            dim == 2 || dim == 3, "point-point velocity is only 2D or 3D");
        // Γ is constant (does not depend on β), so the derivative is zero.
        return Eigen::Vector<T, 2 * dim * dim>::Zero();
    }

    // ========================================================================
    // Point - Edge

    /// @brief Compute the relative velocity of a point and an edge
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param dp Velocity of the point
    /// @param de0 Velocity of the first endpoint of the edge
    /// @param de1 Velocity of the second endpoint of the edge
    /// @param alpha Parametric coordinate of the closest point on the edge
    /// @return The relative velocity of the point and the edge
    template <typename T, int dim>
    inline Eigen::Vector<T, dim> point_edge_relative_velocity(
        Eigen::ConstRef<Eigen::Vector<T, dim>> dp,
        Eigen::ConstRef<Eigen::Vector<T, dim>> de0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> de1,
        const T alpha)
    {
        static_assert(
            dim == 2 || dim == 3, "point-edge velocity is only 2D or 3D");
        return dp - ((de1 - de0) * alpha + de0);
    }

    /// @brief Compute du/dx where u is the relative velocity of a point and an
    /// edge
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param alpha Parametric coordinate of the closest point on the edge
    /// @return The relative velocity Jacobian du/dx
    template <typename T, int dim>
    inline Eigen::Matrix<T, dim, 3 * dim>
    point_edge_relative_velocity_jacobian(const T alpha)
    {
        static_assert(
            dim == 2 || dim == 3, "point-edge velocity is only 2D or 3D");
        Eigen::Matrix<T, dim, 3 * dim> J =
            Eigen::Matrix<T, dim, 3 * dim>::Zero();
        J.template leftCols<dim>().diagonal().setOnes();
        J.template middleCols<dim>(dim).diagonal().setConstant(alpha - T(1));
        J.template rightCols<dim>().diagonal().setConstant(-alpha);
        return J;
    }

    /// @brief Compute d²u/dxdα where u is the relative velocity of a point and
    /// an edge
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param alpha Parametric coordinate of the closest point on the edge
    /// @return The vectorized tensor of d²u/dxdα
    template <typename T, int dim>
    inline Eigen::Vector<T, 3 * dim * dim>
    point_edge_relative_velocity_dx_dbeta(const T alpha)
    {
        static_assert(
            dim == 2 || dim == 3, "point-edge velocity is only 2D or 3D");
        Eigen::Vector<T, 3 * dim * dim> J =
            Eigen::Vector<T, 3 * dim * dim>::Zero();
        for (int i = 0; i < dim; ++i) {
            // I block at cols [dim, 2·dim)
            J[(dim + i) * dim + i] = T(1);
            // -I block at cols [2·dim, 3·dim)
            J[(2 * dim + i) * dim + i] = T(-1);
        }
        return J;
    }

    // ========================================================================
    // Edge - Edge

    /// @brief Compute the relative velocity of the edges.
    /// @tparam T The scalar type.
    /// @param dea0 Velocity of the first endpoint of the first edge
    /// @param dea1 Velocity of the second endpoint of the first edge
    /// @param deb0 Velocity of the first endpoint of the second edge
    /// @param deb1 Velocity of the second endpoint of the second edge
    /// @param coords Two parametric coordinates of the closest points
    /// @return The relative velocity of the edges
    template <typename T>
    inline Eigen::Vector3<T> edge_edge_relative_velocity(
        Eigen::ConstRef<Eigen::Vector3<T>> dea0,
        Eigen::ConstRef<Eigen::Vector3<T>> dea1,
        Eigen::ConstRef<Eigen::Vector3<T>> deb0,
        Eigen::ConstRef<Eigen::Vector3<T>> deb1,
        Eigen::ConstRef<Eigen::Vector2<T>> coords)
    {
        // closest_point_a_velocity - closest_point_b_velocity
        return ((dea1 - dea0) * coords[0] + dea0)
            - ((deb1 - deb0) * coords[1] + deb0);
    }

    /// @brief Compute du/dx where u is the relative velocity of two edges
    /// @tparam T The scalar type.
    /// @param coords Two parametric coordinates of the closest points
    /// @return The relative velocity Jacobian du/dx
    template <typename T>
    inline Eigen::Matrix<T, 3, 12> edge_edge_relative_velocity_jacobian(
        Eigen::ConstRef<Eigen::Vector2<T>> coords)
    {
        Eigen::Matrix<T, 3, 12> J = Eigen::Matrix<T, 3, 12>::Zero();
        J.template leftCols<3>().diagonal().setConstant(T(1) - coords[0]);
        J.template middleCols<3>(3).diagonal().setConstant(coords[0]);
        J.template middleCols<3>(6).diagonal().setConstant(coords[1] - T(1));
        J.template rightCols<3>().diagonal().setConstant(-coords[1]);
        return J;
    }

    /// @brief Compute d²u/dxdβ where u is the relative velocity of two edges
    /// @tparam T The scalar type.
    /// @param coords Two parametric coordinates of the closest points
    /// @return The vectorized tensor of d²u/dxdβ
    template <typename T>
    inline Eigen::Matrix<T, 36, 2> edge_edge_relative_velocity_dx_dbeta(
        Eigen::ConstRef<Eigen::Vector2<T>> coords)
    {
        constexpr int dim = 3;
        Eigen::Matrix<T, 36, 2> J = Eigen::Matrix<T, 36, 2>::Zero();
        for (int i = 0; i < dim; ++i) {
            // wrt β₁: -I at cols [0,3), I at cols [3,6)
            J((0 + i) * dim + i, 0) = T(-1);
            J((3 + i) * dim + i, 0) = T(1);
            // wrt β₂: I at cols [6,9), -I at cols [9,12)
            J((6 + i) * dim + i, 1) = T(1);
            J((9 + i) * dim + i, 1) = T(-1);
        }
        return J;
    }

    // ========================================================================
    // Point - Triangle

    /// @brief Compute the relative velocity of the point to the triangle.
    /// @tparam T The scalar type.
    /// @param dp Velocity of the point
    /// @param dt0 Velocity of the first vertex of the triangle
    /// @param dt1 Velocity of the second vertex of the triangle
    /// @param dt2 Velocity of the third vertex of the triangle
    /// @param coords Baricentric coordinates of the closest point
    /// @return The relative velocity of the point to the triangle
    template <typename T>
    inline Eigen::Vector3<T> point_triangle_relative_velocity(
        Eigen::ConstRef<Eigen::Vector3<T>> dp,
        Eigen::ConstRef<Eigen::Vector3<T>> dt0,
        Eigen::ConstRef<Eigen::Vector3<T>> dt1,
        Eigen::ConstRef<Eigen::Vector3<T>> dt2,
        Eigen::ConstRef<Eigen::Vector2<T>> coords)
    {
        // Compute the velocity of the closest point and subtract it from the
        // points velocity.
        return dp - (dt0 + coords[0] * (dt1 - dt0) + coords[1] * (dt2 - dt0));
    }

    /// @brief Compute du/dx where u is the relative velocity of a point and a
    /// triangle
    /// @tparam T The scalar type.
    /// @param coords Barycentric coordinates of the closest point
    /// @return The relative velocity Jacobian du/dx
    template <typename T>
    inline Eigen::Matrix<T, 3, 12> point_triangle_relative_velocity_jacobian(
        Eigen::ConstRef<Eigen::Vector2<T>> coords)
    {
        Eigen::Matrix<T, 3, 12> J = Eigen::Matrix<T, 3, 12>::Zero();
        J.template leftCols<3>().diagonal().setOnes();
        J.template middleCols<3>(3).diagonal().setConstant(
            coords[0] + coords[1] - T(1));
        J.template middleCols<3>(6).diagonal().setConstant(-coords[0]);
        J.template rightCols<3>().diagonal().setConstant(-coords[1]);
        return J;
    }

    /// @brief Compute d²u/dxdβ where u is the relative velocity of a point and
    /// a triangle
    /// @tparam T The scalar type.
    /// @param coords Baricentric coordinates of the closest point
    /// @return The vectorized tensor of d²u/dxdβ
    template <typename T>
    inline Eigen::Matrix<T, 36, 2> point_triangle_relative_velocity_dx_dbeta(
        Eigen::ConstRef<Eigen::Vector2<T>> coords)
    {
        constexpr int dim = 3;
        Eigen::Matrix<T, 36, 2> J = Eigen::Matrix<T, 36, 2>::Zero();
        for (int i = 0; i < dim; ++i) {
            // wrt β₁: I at cols [3,6), -I at cols [6,9)
            J((3 + i) * dim + i, 0) = T(1);
            J((6 + i) * dim + i, 0) = T(-1);
            // wrt β₂: I at cols [3,6), -I at cols [9,12)
            J((3 + i) * dim + i, 1) = T(1);
            J((9 + i) * dim + i, 1) = T(-1);
        }
        return J;
    }
} // namespace detail

// ============================================================================
// Point - Point

/// @brief Compute the relative velocity of two points
/// @param dp0 Velocity of the first point
/// @param dp1 Velocity of the second point
/// @return The relative velocity of the two points
template <typename DerivedDP0, typename DerivedDP1>
inline auto point_point_relative_velocity(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1)
{
    using T = typename DerivedDP0::Scalar;
    assert(dp0.size() == dp1.size());

    if constexpr (dim_v<DerivedDP0> == 2) {
        return detail::point_point_relative_velocity<T, 2>(dp0, dp1);
    } else if constexpr (dim_v<DerivedDP0> == 3) {
        return detail::point_point_relative_velocity<T, 3>(dp0, dp1);
    } else if (dp0.size() == 2) {
        return VectorMax3<T>(
            detail::point_point_relative_velocity<T, 2>(dp0, dp1));
    } else {
        assert(dp0.size() == 3);
        return VectorMax3<T>(
            detail::point_point_relative_velocity<T, 3>(dp0, dp1));
    }
}

/// @brief Compute du/dx where u is the relative velocity of two points
/// @tparam T The scalar type (defaults to double).
/// @param dim Dimension (2 or 3)
/// @return The relative velocity Jacobian du/dx
template <typename T = double>
inline MatrixMax<T, 3, 6> point_point_relative_velocity_jacobian(const int dim)
{
    if (dim == 2) {
        return MatrixMax<T, 3, 6>(
            detail::point_point_relative_velocity_jacobian<T, 2>());
    } else {
        assert(dim == 3);
        return MatrixMax<T, 3, 6>(
            detail::point_point_relative_velocity_jacobian<T, 3>());
    }
}

/// @brief Compute d²u/dxdβ where u is the relative velocity of two points
/// @tparam T The scalar type (defaults to double).
/// @param dim Dimension (2 or 3)
/// @return The vectorized tensor of d²u/dxdβ
template <typename T = double>
inline VectorMax<T, 18> point_point_relative_velocity_dx_dbeta(const int dim)
{
    if (dim == 2) {
        return VectorMax<T, 18>(
            detail::point_point_relative_velocity_dx_dbeta<T, 2>());
    } else {
        assert(dim == 3);
        return VectorMax<T, 18>(
            detail::point_point_relative_velocity_dx_dbeta<T, 3>());
    }
}

// ============================================================================
// Point - Edge

/// @brief Compute the relative velocity of a point and an edge
/// @param dp Velocity of the point
/// @param de0 Velocity of the first endpoint of the edge
/// @param de1 Velocity of the second endpoint of the edge
/// @param alpha Parametric coordinate of the closest point on the edge
/// @return The relative velocity of the point and the edge
template <typename DerivedDP, typename DerivedDE0, typename DerivedDE1>
inline auto point_edge_relative_velocity(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const typename DerivedDP::Scalar alpha)
{
    using T = typename DerivedDP::Scalar;
    assert(dp.size() == de0.size() && dp.size() == de1.size());

    if constexpr (dim_v<DerivedDP> == 2) {
        return detail::point_edge_relative_velocity<T, 2>(dp, de0, de1, alpha);
    } else if constexpr (dim_v<DerivedDP> == 3) {
        return detail::point_edge_relative_velocity<T, 3>(dp, de0, de1, alpha);
    } else if (dp.size() == 2) {
        return VectorMax3<T>(
            detail::point_edge_relative_velocity<T, 2>(dp, de0, de1, alpha));
    } else {
        assert(dp.size() == 3);
        return VectorMax3<T>(
            detail::point_edge_relative_velocity<T, 3>(dp, de0, de1, alpha));
    }
}

/// @brief Compute du/dx where u is the relative velocity of a point and an edge
/// @tparam T The scalar type (deduced from alpha).
/// @param dim Dimension (2 or 3)
/// @param alpha Parametric coordinate of the closest point on the edge
/// @return The relative velocity Jacobian du/dx
template <typename T = double>
inline MatrixMax<T, 3, 9>
point_edge_relative_velocity_jacobian(const int dim, const T alpha)
{
    if (dim == 2) {
        return MatrixMax<T, 3, 9>(
            detail::point_edge_relative_velocity_jacobian<T, 2>(alpha));
    } else {
        assert(dim == 3);
        return MatrixMax<T, 3, 9>(
            detail::point_edge_relative_velocity_jacobian<T, 3>(alpha));
    }
}

/// @brief Compute d²u/dxdα where u is the relative velocity of a point and an edge
/// @tparam T The scalar type (deduced from alpha).
///
/// Γ(α) = [I, (α-1)I, -αI]  (dim × 3·dim)
/// ∂Γ/∂α = [0, I, -I]
///
/// Stored as vec(∂Γ/∂α) in column-major order (3rd-order convention).
/// Result is a column vector of size dim × 3·dim = 3·dim².
/// @param dim Dimension (2 or 3)
/// @param alpha Parametric coordinate of the closest point on the edge
/// @return The vectorized tensor of d²u/dxdα
template <typename T = double>
inline VectorMax<T, 27>
point_edge_relative_velocity_dx_dbeta(const int dim, const T alpha)
{
    if (dim == 2) {
        return VectorMax<T, 27>(
            detail::point_edge_relative_velocity_dx_dbeta<T, 2>(alpha));
    } else {
        assert(dim == 3);
        return VectorMax<T, 27>(
            detail::point_edge_relative_velocity_dx_dbeta<T, 3>(alpha));
    }
}

// ============================================================================
// Edge - Edge

/// @brief Compute the relative velocity of the edges.
/// @param dea0 Velocity of the first endpoint of the first edge
/// @param dea1 Velocity of the second endpoint of the first edge
/// @param deb0 Velocity of the first endpoint of the second edge
/// @param deb1 Velocity of the second endpoint of the second edge
/// @param coords Two parametric coordinates of the closest points on the edges
/// @return The relative velocity of the edges
template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1,
    typename DerivedCoords>
inline auto edge_edge_relative_velocity(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    using T = typename DerivedDEA0::Scalar;
    // NOTE: explicit <T>: the detail overload cannot deduce T from an
    // arbitrary Eigen expression.
    return detail::edge_edge_relative_velocity<T>(
        dea0, dea1, deb0, deb1, coords);
}

/// @brief Compute du/dx where u is the relative velocity of two edges
/// @param coords Two parametric coordinates of the closest points on the edges
/// @return The relative velocity Jacobian du/dx
template <typename DerivedCoords>
inline auto edge_edge_relative_velocity_jacobian(
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    using T = typename DerivedCoords::Scalar;
    return detail::edge_edge_relative_velocity_jacobian<T>(coords);
}

/// @brief Compute d²u/dxdβ where u is the relative velocity of two edges
///
/// Γ(β₁,β₂) = [(1-β₁)I, β₁I, (β₂-1)I, -β₂I]  (3 × 12)
/// ∂Γ/∂β₁ = [-I, I, 0, 0]
/// ∂Γ/∂β₂ = [ 0, 0, I,-I]
///
/// Stored as [vec(∂Γ/∂β₁) | vec(∂Γ/∂β₂)] in column-major order (3rd-order
/// convention). Result shape: (36, 2). For a (3 × 12) matrix M, element
/// M(r,c) maps to vec index c·3 + r.
/// @param coords Two parametric coordinates of the closest points on the edges
/// @return The vectorized tensor of d²u/dxdβ
template <typename DerivedCoords>
inline auto edge_edge_relative_velocity_dx_dbeta(
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    using T = typename DerivedCoords::Scalar;
    return detail::edge_edge_relative_velocity_dx_dbeta<T>(coords);
}

// ============================================================================
// Point - Triangle

/// @brief Compute the relative velocity of the point to the triangle.
/// @param dp Velocity of the point
/// @param dt0 Velocity of the first vertex of the triangle
/// @param dt1 Velocity of the second vertex of the triangle
/// @param dt2 Velocity of the third vertex of the triangle
/// @param coords Baricentric coordinates of the closest point on the triangle
/// @return The relative velocity of the point to the triangle
template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2,
    typename DerivedCoords>
inline auto point_triangle_relative_velocity(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    using T = typename DerivedDP::Scalar;
    return detail::point_triangle_relative_velocity<T>(
        dp, dt0, dt1, dt2, coords);
}

/// @brief Compute du/dx where u is the relative velocity of a point and a triangle
/// @param coords Barycentric coordinates of the closest point on the triangle
/// @return The relative velocity Jacobian du/dx
template <typename DerivedCoords>
inline auto point_triangle_relative_velocity_jacobian(
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    using T = typename DerivedCoords::Scalar;
    return detail::point_triangle_relative_velocity_jacobian<T>(coords);
}

/// @brief Compute d²u/dxdβ where u is the relative velocity of a point and a triangle
///
/// Γ(β₁,β₂) = [I, (β₁+β₂-1)I, -β₁I, -β₂I]  (3 × 12)
/// ∂Γ/∂β₁ = [0, I, -I, 0]
/// ∂Γ/∂β₂ = [0, I,  0,-I]
///
/// Stored as [vec(∂Γ/∂β₁) | vec(∂Γ/∂β₂)] in column-major order (3rd-order
/// convention). Result shape: (36, 2). For a (3 × 12) matrix M, element
/// M(r,c) maps to vec index c·3 + r.
/// @param coords Baricentric coordinates of the closest point on the triangle
/// @return The vectorized tensor of d²u/dxdβ
template <typename DerivedCoords>
inline auto point_triangle_relative_velocity_dx_dbeta(
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    using T = typename DerivedCoords::Scalar;
    return detail::point_triangle_relative_velocity_dx_dbeta<T>(coords);
}

} // namespace ipc
