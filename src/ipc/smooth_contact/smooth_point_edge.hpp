#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>

namespace ipc {

    /// @brief Compute pointwise potential for a point p and a point specified by uv on edge [e0, e1]
    /// @param p One point outside of edge
    /// @param e0 One end point of the edge
    /// @param e1 One end point of the edge
    /// @param uv Barycentric coordinate
    /// @param dhat The effective distance of barrier
    /// @param alpha The effective angle of barrier
    template <typename scalar>
    scalar smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const Vector2<scalar>>& p,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const double &uv,
        const ParameterType &params);

    /// @brief Compute potential for a point p and an edge [e0, e1], integrated over the edge with high order quadrature
    template <typename scalar>
    scalar smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const Vector2<scalar>>& p,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const ParameterType &params);

    /// @brief Compute potential for a list of points and an edge [e0, e1], integrated over the edge with high order quadrature
    template <typename scalar>
    Vector<scalar, -1, -1> smooth_point_edge_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<scalar, -1, -1, Eigen::RowMajor, -1, 3>>& points,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const ParameterType &params);

    /// @brief Compute potential for a point p and an edge [e0, e1], using the smooth closest point
    template <typename scalar>
    scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const Vector2<scalar>>& p,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const ParameterType &params);
}