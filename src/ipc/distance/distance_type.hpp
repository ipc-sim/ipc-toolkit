#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <cassert>
#include <type_traits>

namespace ipc {

/// @brief Closest pair between a point and point.
enum class PointPointDistanceType : uint8_t {
    P_P = 0, ///< The points are closest to each other.
    AUTO = 0 ///< Automatically determine the closest pair. (Same as P_P)
};

/// @brief Closest pair between a point and edge.
enum class PointEdgeDistanceType : uint8_t {
    P_E0, ///< The point is closest to edge vertex zero.
    P_E1, ///< The point is closest to edge vertex one.
    P_E,  ///< The point is closest to the interior of the edge.
    AUTO  ///< Automatically determine the closest pair.
};

/// @brief Closest pair between a point and triangle.
enum class PointTriangleDistanceType : uint8_t {
    P_T0, ///< The point is closest to triangle vertex zero.
    P_T1, ///< The point is closest to triangle vertex one.
    P_T2, ///< The point is closest to triangle vertex two.
    P_E0, ///< The point is closest to triangle edge zero (vertex zero to one).
    P_E1, ///< The point is closest to triangle edge one (vertex one to two).
    P_E2, ///< The point is closest to triangle edge two (vertex two to zero).
    P_T,  ///< The point is closest to the interior of the triangle.
    AUTO  ///< Automatically determine the closest pair.
};

/// @brief Closest pair between two edges.
enum class EdgeEdgeDistanceType : uint8_t {
    /// The edges are closest at vertex 0 of edge A and 0 of edge B.
    EA0_EB0,
    /// The edges are closest at vertex 0 of edge A and 1 of edge B.
    EA0_EB1,
    /// The edges are closest at vertex 1 of edge A and 0 of edge B.
    EA1_EB0,
    /// The edges are closest at vertex 1 of edge A and 1 of edge B.
    EA1_EB1,
    /// The edges are closest at the interior of edge A and vertex 0 of edge B.
    EA_EB0,
    /// The edges are closest at the interior of edge A and vertex 1 of edge B.
    EA_EB1,
    /// The edges are closest at vertex 0 of edge A and the interior of edge B.
    EA0_EB,
    /// The edges are closest at vertex 1 of edge A and the interior of edge B.
    EA1_EB,
    /// The edges are closest at an interior point of edge A and B.
    EA_EB,
    /// Automatically determine the closest pair.
    AUTO
};

namespace detail {
    /// @brief Warn about a degenerate edge.
    /// @note Out of line to keep the logger (and spdlog) out of this header.
    void warn_degenerate_point_edge() noexcept;

    /// @brief Throw for an invalid distance type.
    /// @note Out of line and [[noreturn]] so that constructing the exception does not consume the caller's inlining budget on the hot path.
    [[noreturn]] void throw_invalid_distance_type(const char* function);

    /// @brief Throw when AUTO is requested for a scalar that cannot resolve it.
    [[noreturn]] void throw_auto_requires_explicit_dtype(const char* function);

    /// @brief Determine the closest pair between a point and edge.
    /// @note Prefer the ipc::point_edge_distance_type front end below, which deduces both the scalar type and the dimension.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge.
    /// @param e1 The second vertex of the edge.
    /// @return The distance type of the point-edge pair.
    template <typename T, int dim>
    inline PointEdgeDistanceType point_edge_distance_type(
        const Eigen::Vector<T, dim>& p,
        const Eigen::Vector<T, dim>& e0,
        const Eigen::Vector<T, dim>& e1)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");

        const Eigen::Vector<T, dim> e = e1 - e0;
        const T e_length_sqr = e.squaredNorm();
        if (e_length_sqr == 0) {
            warn_degenerate_point_edge();
            // WARNING: use an arbitrary end-point
            return PointEdgeDistanceType::P_E0;
        }

        const T ratio = e.dot(p - e0) / e_length_sqr;
        if (ratio < 0) {
            return PointEdgeDistanceType::P_E0; // PP (p-e0)
        } else if (ratio > 1) {
            return PointEdgeDistanceType::P_E1; // PP (p-e1)
        } else {
            return PointEdgeDistanceType::P_E; // PE
        }
    }

    /// @brief Determine the closest pair between a point and triangle.
    /// @param p The point.
    /// @param t0 The first vertex of the triangle.
    /// @param t1 The second vertex of the triangle.
    /// @param t2 The third vertex of the triangle.
    /// @return The distance type of the point-triangle pair.
    template <typename T>
    PointTriangleDistanceType point_triangle_distance_type(
        Eigen::ConstRef<Eigen::Vector3<T>> p,
        Eigen::ConstRef<Eigen::Vector3<T>> t0,
        Eigen::ConstRef<Eigen::Vector3<T>> t1,
        Eigen::ConstRef<Eigen::Vector3<T>> t2);

    /// @brief Determine the closest pair between two edges.
    /// @param ea0 The first vertex of the first edge.
    /// @param ea1 The second vertex of the first edge.
    /// @param eb0 The first vertex of the second edge.
    /// @param eb1 The second vertex of the second edge.
    /// @return The distance type of the edge-edge pair.
    template <typename T>
    EdgeEdgeDistanceType edge_edge_distance_type(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1);

    /// @brief Determine the closest pair between two parallel edges.
    /// @param ea0 The first vertex of the first edge.
    /// @param ea1 The second vertex of the first edge.
    /// @param eb0 The first vertex of the second edge.
    /// @param eb1 The second vertex of the second edge.
    /// @return The distance type of the edge-edge pair.
    template <typename T>
    EdgeEdgeDistanceType edge_edge_parallel_distance_type(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1);
} // namespace detail

// --- EigenExpression wrappers ---
//
// NOTE: The std::is_class_v guard is required. Unlike the wrappers whose
// trailing return type mentions typename DerivedX::Scalar, these return a
// non-dependent enum, so nothing would reject the candidate when the first
// template argument is given explicitly as a scalar (e.g.
// edge_edge_distance_type<double>(...)). Substitution would then form
// Eigen::MatrixBase<double>, a hard error inside Eigen rather than a SFINAE
// failure.

/// @brief Determine the closest pair between a point and edge.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The distance type of the point-edge pair.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    using T = typename DerivedP::Scalar;

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_distance_type<T, 2>(p, e0, e1);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_distance_type<T, 3>(p, e0, e1);
    } else if (p.size() == 2) {
        return detail::point_edge_distance_type<T, 2>(p, e0, e1);
    } else {
        assert(p.size() == 3);
        return detail::point_edge_distance_type<T, 3>(p, e0, e1);
    }
}

/// @brief Determine the closest pair between a point and triangle.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance type of the point-triangle pair.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    std::enable_if_t<std::is_class_v<DerivedP>, int> = 0>
inline PointTriangleDistanceType point_triangle_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    using T = typename DerivedP::Scalar;
    return detail::point_triangle_distance_type<T>(p, t0, t1, t2);
}

/// @brief Determine the closest pair between two edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    std::enable_if_t<std::is_class_v<DerivedEA0>, int> = 0>
inline EdgeEdgeDistanceType edge_edge_distance_type(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_distance_type<T>(ea0, ea1, eb0, eb1);
}

/// @brief Determine the closest pair between two parallel edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    std::enable_if_t<std::is_class_v<DerivedEA0>, int> = 0>
inline EdgeEdgeDistanceType edge_edge_parallel_distance_type(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_parallel_distance_type<T>(ea0, ea1, eb0, eb1);
}

} // namespace ipc
