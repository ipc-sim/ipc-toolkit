#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>
#include "smooth_edge_edge.hpp"
#include <iostream>

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

    /// @brief 
    /// @tparam scalar 
    /// @param p 
    /// @param e0 
    /// @param e1 edge [e0, e1]
    /// @param t0 edge [p, t0]
    /// @param t1 edge [t1, p]
    /// @param params 
    /// @return 
    template <typename scalar>
    scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const Vector2<scalar>>& p,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const Eigen::Ref<const Vector2<scalar>>& x0,
        const Eigen::Ref<const Vector2<scalar>>& x1,
        const ParameterType &params)
    {
        Vector2<scalar> tangent = e1 - e0;
        const scalar len = tangent.norm();
        tangent = tangent / len;

        Vector2<scalar> direc = point_edge_closest_point_direction<scalar>(p, e0, e1, PointEdgeDistanceType::AUTO);
        const scalar dist = direc.norm();
        direc = direc / dist;
        const scalar Phi = 1 - cross2<scalar>(p - e0, tangent) / dist; // intpow(diff.dot(tangent), 2) / dist_sqr;

        Vector2<scalar> t0 = x0 - p, t1 = p - x1;
        scalar l0 = t0.norm(), l1 = t1.norm();
        scalar tangent_term = smooth_heaviside<scalar>(t0.dot(direc) / l0 / params.alpha) *
                            smooth_heaviside<scalar>(t1.dot(direc) / l1 / params.alpha);

        const scalar mollifier_val = edge_mollifier<scalar>(p, e0, e1, intpow(dist, 2));

        return 0.5 * (l0 + l1) * len * cubic_spline(Phi * (2. / params.alpha)) * inv_barrier(intpow(dist, 2) / params.eps, params.r) * tangent_term * mollifier_val;
    }

    template <typename scalar>
    scalar smooth_point_edge_potential_single_point_3d(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const ParameterType &params)
    {
        const Vector3<scalar> u = e1 - e0;
        const scalar len = u.squaredNorm();

        Vector3<scalar> direc = point_edge_closest_point_direction<scalar>(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const scalar dist_sqr = direc.squaredNorm();
        direc = direc / sqrt(dist_sqr);

        scalar out = smooth_edge_edge_potential_tangent_term<scalar>(f0, e0, e1, -direc, params.alpha, HEAVISIDE_TYPE::VARIANT) *
                smooth_edge_edge_potential_tangent_term<scalar>(f1, e1, e0, -direc, params.alpha, HEAVISIDE_TYPE::VARIANT) *
                inv_barrier<scalar>(dist_sqr / params.eps, params.r);

        const scalar normal_penalty = smooth_edge_edge_potential_normal_term<scalar>(f0, e0, e1, direc, params.alpha, HEAVISIDE_TYPE::VARIANT) + 
                                    smooth_edge_edge_potential_normal_term<scalar>(f1, e1, e0, direc, params.alpha, HEAVISIDE_TYPE::VARIANT);

        const scalar mollifier_val = edge_mollifier<scalar>(p, e0, e1, dist_sqr);

        return 0.5 * sqrt(len) * out * normal_penalty * mollifier_val;
    }
}