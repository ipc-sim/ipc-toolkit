#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Cholesky>

#include <array>
#include <cassert>
#include <limits>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    /// hess is (6×6) flattened in column-major order
    template <typename T>
    void point_edge_closest_point_2D_hessian(
        T p_x, T p_y, T e0_x, T e0_y, T e1_x, T e1_y, T hess[36]);

    /// hess is (9×9) flattened in column-major order
    template <typename T>
    void point_edge_closest_point_3D_hessian(
        T p_x,
        T p_y,
        T p_z,
        T e0_x,
        T e0_y,
        T e0_z,
        T e1_x,
        T e1_y,
        T e1_z,
        T hess[81]);

    /// J is (2×12) flattened in column-major order
    template <typename T>
    void edge_edge_closest_point_jacobian(
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
        T J[24]);

    /// hess is (144×1) flattened in column-major order
    template <typename T>
    void edge_edge_closest_point_hessian_a(
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
        T hess[144]);

    /// hess is (144×1) flattened in column-major order
    template <typename T>
    void edge_edge_closest_point_hessian_b(
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
        T hess[144]);

    /// J is (2×12) flattened in column-major order
    template <typename T>
    void point_triangle_closest_point_jacobian(
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
        T J[24]);

    /// hess is (144×1) flattened in column-major order
    template <typename T>
    void point_triangle_closest_point_hessian_0(
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
        T hess[144]);

    /// hess is (144×1) flattened in column-major order
    template <typename T>
    void point_triangle_closest_point_hessian_1(
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
        T hess[144]);
} // namespace autogen

namespace detail {
    /// @brief Residual bound for the 2×2 solves below.
    ///
    /// The 1e-10 bound is tuned for double, so scale it by the relative
    /// precision of T; otherwise the float instantiation would be held to a
    /// double-precision bound.
    template <typename T>
    inline constexpr double CLOSEST_POINT_RESIDUAL_TOL = 1e-10
        * (static_cast<double>(std::numeric_limits<T>::epsilon())
           / std::numeric_limits<double>::epsilon());

    // ========================================================================
    // Point - Edge

    /// @brief Compute the barycentric coordinate of the closest point on the edge.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p Point
    /// @param e0 First edge point
    /// @param e1 Second edge point
    /// @return Barycentric coordinate of the closest point
    template <typename T, int dim>
    inline T point_edge_closest_point(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");
        const Eigen::Vector<T, dim> e = e1 - e0;
        return (p - e0).dot(e) / e.squaredNorm();
    }

    /// @brief Compute the Jacobian of the closest point on the edge.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p Point
    /// @param e0 First edge point
    /// @param e1 Second edge point
    /// @return Jacobian of the closest point
    template <typename T, int dim>
    inline Eigen::Vector<T, 3 * dim> point_edge_closest_point_jacobian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");

        const Eigen::Vector<T, dim> e = e1 - e0;
        const Eigen::Vector<T, dim> e2p = p - e0;
        const T e_sqnorm = e.squaredNorm();

        Eigen::Vector<T, 3 * dim> J;
        J.template head<dim>() = e / e_sqnorm;
        J.template segment<dim>(dim) =
            (T(2) / e_sqnorm * e.dot(e2p) * e - e - e2p) / e_sqnorm;
        J.template tail<dim>() =
            (e2p - T(2) / e_sqnorm * e.dot(e2p) * e) / e_sqnorm;
        return J;
    }

    /// @brief Compute the Hessian of the closest point on the edge.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p Point
    /// @param e0 First edge point
    /// @param e1 Second edge point
    /// @return Hessian of the closest point
    template <typename T, int dim>
    inline Eigen::Matrix<T, 3 * dim, 3 * dim> point_edge_closest_point_hessian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");

        Eigen::Matrix<T, 3 * dim, 3 * dim> H;
        if constexpr (dim == 2) {
            autogen::point_edge_closest_point_2D_hessian(
                p(0), p(1), e0(0), e0(1), e1(0), e1(1), H.data());
        } else {
            autogen::point_edge_closest_point_3D_hessian(
                p(0), p(1), p(2), e0(0), e0(1), e0(2), e1(0), e1(1), e1(2),
                H.data());
        }
        return H;
    }

    // ========================================================================
    // Edge - Edge

    /// @brief Compute the barycentric coordinates of the closest points between two edges.
    /// @tparam T The scalar type.
    /// @param ea0 First point of the first edge
    /// @param ea1 Second point of the first edge
    /// @param eb0 First point of the second edge
    /// @param eb1 Second point of the second edge
    /// @return Barycentric coordinates of the closest points
    template <typename T>
    inline Eigen::Vector2<T> edge_edge_closest_point(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        const Eigen::Vector3<T> eb_to_ea = ea0 - eb0;
        const Eigen::Vector3<T> ea = ea1 - ea0;
        const Eigen::Vector3<T> eb = eb1 - eb0;

        Eigen::Matrix2<T> A;
        A(0, 0) = ea.squaredNorm();
        A(0, 1) = A(1, 0) = -eb.dot(ea);
        A(1, 1) = eb.squaredNorm();

        Eigen::Vector2<T> rhs;
        rhs[0] = -eb_to_ea.dot(ea);
        rhs[1] = eb_to_ea.dot(eb);

        const Eigen::Vector2<T> x = A.ldlt().solve(rhs);
        assert((A * x - rhs).norm() < CLOSEST_POINT_RESIDUAL_TOL<T>);
        return x;
    }

    /// @brief Compute the Jacobian of the closest points between two edges.
    /// @tparam T The scalar type.
    /// @param ea0 First point of the first edge
    /// @param ea1 Second point of the first edge
    /// @param eb0 First point of the second edge
    /// @param eb1 Second point of the second edge
    /// @return Jacobian of the closest points
    template <typename T>
    inline Eigen::Matrix<T, 2, 12> edge_edge_closest_point_jacobian(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        Eigen::Matrix<T, 2, 12> J;
        autogen::edge_edge_closest_point_jacobian(
            ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1],
            eb0[2], eb1[0], eb1[1], eb1[2], J.data());
        return J;
    }

    /// @brief Compute the Hessian of the closest points between two edges.
    /// @tparam T The scalar type.
    /// @param ea0 First point of the first edge
    /// @param ea1 Second point of the first edge
    /// @param eb0 First point of the second edge
    /// @param eb1 Second point of the second edge
    /// @return Hessian of the closest points (2x12x12 tensor)
    template <typename T>
    inline std::array<Eigen::Matrix<T, 12, 12>, 2>
    edge_edge_closest_point_hessian(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        std::array<Eigen::Matrix<T, 12, 12>, 2> H;
        autogen::edge_edge_closest_point_hessian_a(
            ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1],
            eb0[2], eb1[0], eb1[1], eb1[2], H[0].data());
        autogen::edge_edge_closest_point_hessian_b(
            ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1],
            eb0[2], eb1[0], eb1[1], eb1[2], H[1].data());
        return H;
    }

    // ========================================================================
    // Point - Triangle

    /// @brief Compute the barycentric coordinates of the closest point on the triangle.
    /// @tparam T The scalar type.
    /// @param p Point
    /// @param t0 Triangle's first vertex
    /// @param t1 Triangle's second vertex
    /// @param t2 Triangle's third vertex
    /// @return Barycentric coordinates of the closest point
    template <typename T>
    inline Eigen::Vector2<T> point_triangle_closest_point(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2)
    {
        Eigen::Matrix<T, 2, 3> basis;
        basis.row(0) = Eigen::RowVector3<T>(t1 - t0); // edge 0
        basis.row(1) = Eigen::RowVector3<T>(t2 - t0); // edge 1
        const Eigen::Matrix2<T> A = basis * basis.transpose();
        const Eigen::Vector2<T> b = basis * (p - t0);
        const Eigen::Vector2<T> x = A.ldlt().solve(b);
        assert((A * x - b).norm() < CLOSEST_POINT_RESIDUAL_TOL<T>);
        return x;
    }

    /// @brief Compute the Jacobian of the closest point on the triangle.
    /// @tparam T The scalar type.
    /// @param p Point
    /// @param t0 Triangle's first vertex
    /// @param t1 Triangle's second vertex
    /// @param t2 Triangle's third vertex
    /// @return Jacobian of the closest point
    template <typename T>
    inline Eigen::Matrix<T, 2, 12> point_triangle_closest_point_jacobian(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2)
    {
        Eigen::Matrix<T, 2, 12> J;
        autogen::point_triangle_closest_point_jacobian(
            p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
            t2[1], t2[2], J.data());
        return J;
    }

    /// @brief Compute the Hessian of the closest point on the triangle.
    /// @tparam T The scalar type.
    /// @param p Point
    /// @param t0 Triangle's first vertex
    /// @param t1 Triangle's second vertex
    /// @param t2 Triangle's third vertex
    /// @return Hessian of the closest point (2x12x12 tensor)
    template <typename T>
    inline std::array<Eigen::Matrix<T, 12, 12>, 2>
    point_triangle_closest_point_hessian(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2)
    {
        std::array<Eigen::Matrix<T, 12, 12>, 2> H;
        autogen::point_triangle_closest_point_hessian_0(
            p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
            t2[1], t2[2], H[0].data());
        autogen::point_triangle_closest_point_hessian_1(
            p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
            t2[1], t2[2], H[1].data());
        return H;
    }
} // namespace detail

// ============================================================================
// Point - Edge

/// @brief Compute the baricentric coordinate of the closest point on the edge.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return barycentric coordinates of the closest point
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_closest_point(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    using T = typename DerivedP::Scalar;

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_closest_point<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_closest_point<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        return detail::point_edge_closest_point<T, 2>(p, e0, e1);
    } else {
        assert(p.size() == 3);
        return detail::point_edge_closest_point<T, 3>(p, e0, e1);
    }
}

/// @brief Compute the Jacobian of the closest point on the edge.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return Jacobian of the closest point
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_closest_point_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    using T = typename DerivedP::Scalar;

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_closest_point_jacobian<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_closest_point_jacobian<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        return VectorMax9<T>(
            detail::point_edge_closest_point_jacobian<T, 2>(p, e0, e1));
    } else {
        assert(p.size() == 3);
        return VectorMax9<T>(
            detail::point_edge_closest_point_jacobian<T, 3>(p, e0, e1));
    }
}

/// @brief Compute the hessian of the closest point on the edge
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return Hessian of the closest point
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_closest_point_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    using T = typename DerivedP::Scalar;

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_closest_point_hessian<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_closest_point_hessian<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        MatrixMax9<T> H(6, 6);
        autogen::point_edge_closest_point_2D_hessian(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], H.data());
        return H;
    } else {
        assert(p.size() == 3);
        MatrixMax9<T> H(9, 9);
        autogen::point_edge_closest_point_3D_hessian(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            H.data());
        return H;
    }
}

// ============================================================================
// Edge - Edge

/// @brief Compute the barycentric coordinates of the closest points between two edges.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return Barycentric coordinates of the closest points
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_closest_point(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_closest_point<T>(ea0, ea1, eb0, eb1);
}

/// @brief Compute the Jacobian of the closest points between two edges.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return Jacobian of the closest points
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_closest_point_jacobian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_closest_point_jacobian<T>(ea0, ea1, eb0, eb1);
}

/// @brief Compute the Hessian of the closest points between two edges.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return Hessian of the closest points (2x12x12 tensor)
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_closest_point_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_closest_point_hessian<T>(ea0, ea1, eb0, eb1);
}

// ============================================================================
// Point - Triangle

/// @brief Compute the barycentric coordinates of the closest point on the triangle.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return Barycentric coordinates of the closest point
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_closest_point(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_closest_point<T>(p, t0, t1, t2);
}

/// @brief Compute the Jacobian of the closest point on the triangle.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return Jacobian of the closest point
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_closest_point_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_closest_point_jacobian<T>(p, t0, t1, t2);
}

/// @brief Compute the Hessian of the closest point on the triangle.
///
/// Accepts any 3D Eigen expression; the scalar type is deduced.
///
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return Hessian of the closest point (2x12x12 tensor)
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_closest_point_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_closest_point_hessian<T>(p, t0, t1, t2);
}

} // namespace ipc
