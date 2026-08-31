#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <cassert>
#include <type_traits>

namespace ipc {

namespace detail {
    /// @brief Compute the distance between a point and edge in 2D or 3D.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge.
    /// @param e1 The second vertex of the edge.
    /// @param dtype The point edge distance type to compute.
    /// @return The distance between the point and edge.
    template <typename T, int dim>
    inline T point_edge_distance(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1,
        PointEdgeDistanceType dtype)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");

        if constexpr (std::is_floating_point_v<T>) {
            if (dtype == PointEdgeDistanceType::AUTO) {
                dtype = point_edge_distance_type<T, dim>(p, e0, e1);
            }
        } else if (dtype == PointEdgeDistanceType::AUTO) {
            throw_auto_requires_explicit_dtype("point_edge_distance");
        }

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            return point_point_distance(p, e0);

        case PointEdgeDistanceType::P_E1:
            return point_point_distance(p, e1);

        case PointEdgeDistanceType::P_E:
            return point_line_distance(p, e0, e1);

        default:
            throw_invalid_distance_type("point_edge_distance");
        }
    }

    /// @brief Compute the gradient of the distance between a point and edge.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge.
    /// @param e1 The second vertex of the edge.
    /// @param dtype The point edge distance type to compute.
    /// @return The gradient of the distance wrt p, e0, and e1.
    template <typename T, int dim>
    inline Eigen::Vector<T, 3 * dim> point_edge_distance_gradient(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1,
        PointEdgeDistanceType dtype)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");

        if constexpr (std::is_floating_point_v<T>) {
            if (dtype == PointEdgeDistanceType::AUTO) {
                dtype = point_edge_distance_type<T, dim>(p, e0, e1);
            }
        } else if (dtype == PointEdgeDistanceType::AUTO) {
            throw_auto_requires_explicit_dtype("point_edge_distance_gradient");
        }

        Eigen::Vector<T, 3 * dim> grad = Eigen::Vector<T, 3 * dim>::Zero();

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            grad.template head<2 * dim>() =
                point_point_distance_gradient<T, dim>(p, e0);
            break;

        case PointEdgeDistanceType::P_E1: {
            const Eigen::Vector<T, 2 * dim> local_grad =
                point_point_distance_gradient<T, dim>(p, e1);
            grad.template head<dim>() = local_grad.template head<dim>();
            grad.template tail<dim>() = local_grad.template tail<dim>();
            break;
        }

        case PointEdgeDistanceType::P_E:
            grad = point_line_distance_gradient<T, dim>(p, e0, e1);
            break;

        default:
            throw_invalid_distance_type("point_edge_distance_gradient");
        }
        return grad;
    }

    /// @brief Compute the hessian of the distance between a point and edge.
    /// @note The distance is actually squared distance.
    /// @tparam T The scalar type.
    /// @tparam dim The dimension (2 or 3).
    /// @param p The point.
    /// @param e0 The first vertex of the edge.
    /// @param e1 The second vertex of the edge.
    /// @param dtype The point edge distance type to compute.
    /// @return The hessian of the distance wrt p, e0, and e1.
    template <typename T, int dim>
    Eigen::Matrix<T, 3 * dim, 3 * dim> point_edge_distance_hessian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1,
        PointEdgeDistanceType dtype);
} // namespace detail

/// @brief Compute the distance between a point and edge in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The distance between the point and edge.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_distance(
    const DerivedP& p,
    const DerivedE0& e0,
    const DerivedE1& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedE0, DerivedE1);
    using T = typename DerivedP::Scalar;

    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_distance<T, 2>(p, e0, e1, dtype);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_distance<T, 3>(p, e0, e1, dtype);
    } else if (p.size() == 2) {
        return detail::point_edge_distance<T, 2>(p, e0, e1, dtype);
    } else {
        assert(p.size() == 3);
        return detail::point_edge_distance<T, 3>(p, e0, e1, dtype);
    }
}

/// @brief Compute the gradient of the distance between a point and edge.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The gradient of the distance wrt p, e0, and e1.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_distance_gradient(
    const DerivedP& p,
    const DerivedE0& e0,
    const DerivedE1& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedE0, DerivedE1);
    using T = typename DerivedP::Scalar;
    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_distance_gradient<T, 2>(p, e0, e1, dtype);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_distance_gradient<T, 3>(p, e0, e1, dtype);
    } else if (p.size() == 2) {
        return VectorMax9<T>(
            detail::point_edge_distance_gradient<T, 2>(p, e0, e1, dtype));
    } else {
        assert(p.size() == 3);
        return VectorMax9<T>(
            detail::point_edge_distance_gradient<T, 3>(p, e0, e1, dtype));
    }
}

/// @brief Compute the hessian of the distance between a point and edge.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The hessian of the distance wrt p, e0, and e1.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_distance_hessian(
    const DerivedP& p,
    const DerivedE0& e0,
    const DerivedE1& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedP, DerivedE0, DerivedE1);
    using T = typename DerivedP::Scalar;
    if constexpr (dim_v<DerivedP> == 2) {
        return detail::point_edge_distance_hessian<T, 2>(p, e0, e1, dtype);
    } else if constexpr (dim_v<DerivedP> == 3) {
        return detail::point_edge_distance_hessian<T, 3>(p, e0, e1, dtype);
    } else if (p.size() == 2) {
        return MatrixMax9<T>(
            detail::point_edge_distance_hessian<T, 2>(p, e0, e1, dtype));
    } else {
        assert(p.size() == 3);
        return MatrixMax9<T>(
            detail::point_edge_distance_hessian<T, 3>(p, e0, e1, dtype));
    }
}

} // namespace ipc
