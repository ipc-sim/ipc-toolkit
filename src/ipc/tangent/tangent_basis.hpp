#pragma once

#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/simd.hpp>

#include <cassert>

namespace ipc {

namespace detail {
    // ========================================================================
    // Point - Point

    /// @brief Compute a basis for the space tangent to the point-point pair.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p0 First point
    /// @param p1 Second point
    /// @return A dim×(dim-1) matrix whose columns are the basis vectors.
    template <typename T, int dim>
    inline Eigen::Matrix<T, dim, dim - 1> point_point_tangent_basis(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> p1)
    {
        static_assert(
            dim == 2 || dim == 3, "point-point tangent basis is only 2D or 3D");

        if constexpr (dim == 2) {
            const Eigen::Vector2<T> p0_to_p1 = normalized(p1 - p0);
            return Eigen::Vector2<T>(-p0_to_p1.y(), p0_to_p1.x());
        } else {
            const Eigen::Vector3<T> p0_to_p1 = p1 - p0;

            const Eigen::Vector3<T> cross_x =
                Eigen::Vector3<T>::UnitX().cross(p0_to_p1);
            const Eigen::Vector3<T> cross_y =
                Eigen::Vector3<T>::UnitY().cross(p0_to_p1);

            // Prefer whichever reference axis is least parallel to the pair.
            // A batch cannot answer that with one bool, so it builds both
            // bases and blends them per-lane.
            Eigen::Matrix<T, 3, 2> basis;
            if constexpr (is_simd_batch_v<T>) {
                Eigen::Matrix<T, 3, 2> basis_x, basis_y;
                basis_x.col(0) = normalized(cross_x);
                basis_x.col(1) = normalized(p0_to_p1.cross(cross_x));
                basis_y.col(0) = normalized(cross_y);
                basis_y.col(1) = normalized(p0_to_p1.cross(cross_y));

                const auto prefer_x =
                    cross_x.squaredNorm() > cross_y.squaredNorm();
                for (Eigen::Index i = 0; i < basis.size(); ++i) {
                    basis(i) = select(prefer_x, basis_x(i), basis_y(i));
                }
            } else if (cross_x.squaredNorm() > cross_y.squaredNorm()) {
                basis.col(0) = cross_x.normalized();
                basis.col(1) = p0_to_p1.cross(cross_x).normalized();
            } else {
                basis.col(0) = cross_y.normalized();
                basis.col(1) = p0_to_p1.cross(cross_y).normalized();
            }

            return basis;
        }
    }

    /// @brief Compute the Jacobian of the tangent basis for the point-point
    /// pair.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p0 First point
    /// @param p1 Second point
    /// @return A (dim*(dim-1))×(2*dim) matrix.
    template <typename T, int dim>
    Eigen::Matrix<T, dim*(dim - 1), 2 * dim> point_point_tangent_basis_jacobian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> p1);

    // ========================================================================
    // Point - Edge

    /// @brief Compute a basis for the space tangent to the point-edge pair.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p Point
    /// @param e0 First edge point
    /// @param e1 Second edge point
    /// @return A dim×(dim-1) matrix whose columns are the basis vectors.
    template <typename T, int dim>
    inline Eigen::Matrix<T, dim, dim - 1> point_edge_tangent_basis(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(
            dim == 2 || dim == 3, "point-edge tangent basis is only 2D or 3D");

        if constexpr (dim == 2) {
            return normalized(e1 - e0);
        } else {
            const Eigen::Vector3<T> e = e1 - e0;

            Eigen::Matrix<T, 3, 2> basis;
            basis.col(0) = normalized(e);
            basis.col(1) = normalized(e.cross(p - e0));
            return basis;
        }
    }

    /// @brief Compute the Jacobian of the tangent basis for the point-edge pair.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p Point
    /// @param e0 First edge point
    /// @param e1 Second edge point
    /// @return A (dim*(dim-1))×(3*dim) matrix.
    template <typename T, int dim>
    Eigen::Matrix<T, dim*(dim - 1), 3 * dim> point_edge_tangent_basis_jacobian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1);

    // ========================================================================
    // Edge - Edge

    /// @brief Compute a basis for the space tangent to the edge-edge pair.
    /// @tparam T The scalar type.
    /// @param ea0 First point of the first edge
    /// @param ea1 Second point of the first edge
    /// @param eb0 First point of the second edge
    /// @param eb1 Second point of the second edge
    /// @return A 3x2 matrix whose columns are the basis vectors.
    template <typename T>
    inline Eigen::Matrix<T, 3, 2> edge_edge_tangent_basis(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        const Eigen::Vector3<T> ea = ea1 - ea0; // Edge A direction
        const Eigen::Vector3<T> normal = ea.cross(eb1 - eb0);
        // The normal will be zero if the edges are parallel (i.e. coplanar).
        if constexpr (std::is_floating_point_v<T>) {
            assert(normal.norm() != 0);
        }

        Eigen::Matrix<T, 3, 2> basis;
        // The first basis vector is along edge A.
        basis.col(0) = normalized(ea);
        // The second basis vector is orthogonal to the first and the edge-edge
        // normal.
        basis.col(1) = normalized(normal.cross(ea));
        return basis;
    }

    /// @brief Compute the Jacobian of the tangent basis for the edge-edge pair.
    /// @tparam T The scalar type.
    /// @param ea0 First point of the first edge
    /// @param ea1 Second point of the first edge
    /// @param eb0 First point of the second edge
    /// @param eb1 Second point of the second edge
    /// @return A (3*2)x12 matrix whose columns are the basis vectors.
    template <typename T>
    Eigen::Matrix<T, 6, 12> edge_edge_tangent_basis_jacobian(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1);

    // ========================================================================
    // Point - Triangle

    /// @brief Compute a basis for the space tangent to the point-triangle pair.
    /// @tparam T The scalar type.
    /// @param p Point
    /// @param t0 Triangle's first vertex
    /// @param t1 Triangle's second vertex
    /// @param t2 Triangle's third vertex
    /// @return A 3x2 matrix whose columns are the basis vectors.
    template <typename T>
    inline Eigen::Matrix<T, 3, 2> point_triangle_tangent_basis(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2)
    {
        const Eigen::Vector3<T> e0 = t1 - t0;
        const Eigen::Vector3<T> normal = e0.cross(t2 - t0);
        if constexpr (std::is_floating_point_v<T>) {
            assert(normal.norm() != 0);
        }

        Eigen::Matrix<T, 3, 2> basis;

        // The first basis vector is along first edge of the triangle.
        basis.col(0) = normalized(e0);
        // The second basis vector is orthogonal to the first and the triangle
        // normal.
        basis.col(1) = normalized(normal.cross(e0));

        return basis;
    }

    /// @brief Compute the Jacobian of the tangent basis for the point-triangle pair.
    /// @tparam T The scalar type.
    /// @param p Point
    /// @param t0 Triangle's first vertex
    /// @param t1 Triangle's second vertex
    /// @param t2 Triangle's third vertex
    /// @return A (3*2)x12 matrix whose columns are the basis vectors.
    template <typename T>
    Eigen::Matrix<T, 6, 12> point_triangle_tangent_basis_jacobian(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2);
} // namespace detail

// ============================================================================
// Point - Point

/// @brief Compute a basis for the space tangent to the point-point pair.
/// @param p0 First point
/// @param p1 Second point
/// @return A 3x2 matrix whose columns are the basis vectors.
template <typename DerivedP0, typename DerivedP1>
inline auto point_point_tangent_basis(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    using T = typename DerivedP0::Scalar;
    assert(p0.size() == p1.size());

    if constexpr (dim_v<DerivedP0> == 2) {
        return detail::point_point_tangent_basis<T, 2>(p0, p1);
    } else if constexpr (dim_v<DerivedP0> == 3) {
        return detail::point_point_tangent_basis<T, 3>(p0, p1);
    } else if (p0.size() == 2) {
        return MatrixMax<T, 3, 2>(
            detail::point_point_tangent_basis<T, 2>(p0, p1));
    } else {
        assert(p0.size() == 3);
        return MatrixMax<T, 3, 2>(
            detail::point_point_tangent_basis<T, 3>(p0, p1));
    }
}

/// @brief Compute the Jacobian of the tangent basis for the point-point pair.
/// @param p0 First point
/// @param p1 Second point
/// @return A (3*2)x6 matrix whose columns are the basis vectors.
template <typename DerivedP0, typename DerivedP1>
inline auto point_point_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    using T = typename DerivedP0::Scalar;
    assert(p0.size() == p1.size());

    if constexpr (dim_v<DerivedP0> == 2) {
        return detail::point_point_tangent_basis_jacobian<T, 2>(p0, p1);
    } else if constexpr (dim_v<DerivedP0> == 3) {
        return detail::point_point_tangent_basis_jacobian<T, 3>(p0, p1);
    } else if (p0.size() == 2) {
        return MatrixMax<T, 6, 6>(
            detail::point_point_tangent_basis_jacobian<T, 2>(p0, p1));
    } else {
        assert(p0.size() == 3);
        return MatrixMax<T, 6, 6>(
            detail::point_point_tangent_basis_jacobian<T, 3>(p0, p1));
    }
}

// ============================================================================
// Point - Edge

/// @brief Compute a basis for the space tangent to the point-edge pair.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return A 3x2 matrix whose columns are the basis vectors.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    using T = typename DerivedP::Scalar;
    assert(p.size() == e0.size() && p.size() == e1.size());

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_tangent_basis<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_tangent_basis<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        return MatrixMax<T, 3, 2>(
            detail::point_edge_tangent_basis<T, 2>(p, e0, e1));
    } else {
        assert(p.size() == 3);
        return MatrixMax<T, 3, 2>(
            detail::point_edge_tangent_basis<T, 3>(p, e0, e1));
    }
}

/// @brief Compute the Jacobian of the tangent basis for the point-edge pair.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return A (3*2)x9 matrix whose columns are the basis vectors.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    using T = typename DerivedP::Scalar;
    assert(p.size() == e0.size() && p.size() == e1.size());

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_tangent_basis_jacobian<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_tangent_basis_jacobian<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        return MatrixMax<T, 6, 9>(
            detail::point_edge_tangent_basis_jacobian<T, 2>(p, e0, e1));
    } else {
        assert(p.size() == 3);
        return MatrixMax<T, 6, 9>(
            detail::point_edge_tangent_basis_jacobian<T, 3>(p, e0, e1));
    }
}

// ============================================================================
// Edge - Edge

/// @brief Compute a basis for the space tangent to the edge-edge pair.
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return A 3x2 matrix whose columns are the basis vectors.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    // NOTE: explicit <T>: the detail overload cannot deduce T from an
    // arbitrary Eigen expression.
    return detail::edge_edge_tangent_basis<T>(ea0, ea1, eb0, eb1);
}

/// @brief Compute the Jacobian of the tangent basis for the edge-edge pair.
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return A (3*2)x12 matrix whose columns are the basis vectors.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_tangent_basis_jacobian<T>(ea0, ea1, eb0, eb1);
}

// ============================================================================
// Point - Triangle

/// @brief Compute a basis for the space tangent to the point-triangle pair.
///
/// \f\[
///     \begin{bmatrix}
///     \frac{t_1 - t_0}{\|t_1 - t_0\|} & \frac{((t_1 - t_0)\times(t_2 - t_0))
///     \times(t_1 - t_0)}{\|((t_1 - t_0)\times(t_2 - t_0))\times(t_1 - t_0)\|}
///     \end{bmatrix}
/// \f\]
///
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return A 3x2 matrix whose columns are the basis vectors.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_tangent_basis<T>(p, t0, t1, t2);
}

/// @brief Compute the Jacobian of the tangent basis for the point-triangle pair.
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return A (3*2)x12 matrix whose columns are the basis vectors.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_tangent_basis_jacobian<T>(p, t0, t1, t2);
}

// ============================================================================
// Symbolically generated derivatives

namespace autogen {

    // J is (2×4) flattened in column-major order
    template <typename T>
    void point_point_tangent_basis_2D_jacobian(
        T p0_x, T p0_y, T p1_x, T p1_y, T J[8]);

    // J is (6×6) flattened in column-major order
    template <typename T>
    void point_point_tangent_basis_3D_jacobian(
        T p0_x, T p0_y, T p0_z, T p1_x, T p1_y, T p1_z, T J[36]);

    // J is (2×6) flattened in column-major order
    template <typename T>
    void point_edge_tangent_basis_2D_jacobian(
        T p_x, T p_y, T e0_x, T e0_y, T e1_x, T e1_y, T J[12]);

    // J is (6×9) flattened in column-major order
    template <typename T>
    void point_edge_tangent_basis_3D_jacobian(
        T p_x,
        T p_y,
        T p_z,
        T e0_x,
        T e0_y,
        T e0_z,
        T e1_x,
        T e1_y,
        T e1_z,
        T J[54]);

    // J is (6×12) flattened in column-major order
    template <typename T>
    void edge_edge_tangent_basis_jacobian(
        T ea0_x,
        T ea0_y,
        T ea0_z,
        T ea1_x,
        T ea1_y,
        T ea1_z,
        T eb0_x,
        T eb0_y,
        T eb0_z,
        T eb1_x,
        T eb1_y,
        T eb1_z,
        T J[72]);

    // J is (6×12) flattened in column-major order
    template <typename T>
    void point_triangle_tangent_basis_jacobian(
        T p_x,
        T p_y,
        T p_z,
        T t0_x,
        T t0_y,
        T t0_z,
        T t1_x,
        T t1_y,
        T t1_z,
        T t2_x,
        T t2_y,
        T t2_z,
        T J[72]);
} // namespace autogen

} // namespace ipc
