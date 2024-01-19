#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>
#include "smooth_edge_edge.hpp"
#include "smooth_point_point.hpp"

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
                            smooth_heaviside<scalar>(-t1.dot(direc) / l1 / params.alpha);

        scalar normal_term = smooth_heaviside<scalar>(-cross2<scalar>(t0, direc) / l0 / params.alpha) + 
                            smooth_heaviside<scalar>(-cross2<scalar>(t1, direc) / l1 / params.alpha);

        const scalar mollifier_val = edge_mollifier<scalar>(p, e0, e1, intpow(dist, 2));

        return 0.5 * (l0 + l1) * len * cubic_spline(Phi * (2. / params.alpha)) * inv_barrier(intpow(dist, 2) / params.eps, params.r) * tangent_term * normal_term * mollifier_val;
    }

    inline bool smooth_point_edge_potential_single_point_3d_type(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Matrix<double, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const ParameterType &params)
    {
        Vector3<double> direc = point_edge_closest_point_direction<double>(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const double dist_sqr = direc.squaredNorm();
        if (dist_sqr >= params.eps)
            return false;
        
        direc = direc / sqrt(dist_sqr);

        const bool edge_term = smooth_edge3_term_type(direc, e0, e1, f0, f1, params.alpha);
        const bool vert_term = smooth_point3_term_type(p, direc, neighbors, params.alpha);

        return edge_term && vert_term;
    }

    template <typename scalar>
    scalar smooth_point_edge_potential_single_point_3d(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Matrix<scalar, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const ParameterType &params)
    {
        Vector3<scalar> direc = point_edge_closest_point_direction<scalar>(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const scalar dist_sqr = direc.squaredNorm();
        direc = direc / sqrt(dist_sqr);

        const scalar edge_term = smooth_edge3_term<scalar>(direc, e0, e1, f0, f1, params.alpha);
        const scalar barrier = inv_barrier<scalar>(dist_sqr / params.eps, params.r);
        const scalar mollifier_val = edge_mollifier<scalar>(p, e0, e1, dist_sqr);
        const scalar vert_term = smooth_point3_term<scalar>(p, direc, neighbors, params.alpha);

        return (e1 - e0).norm() * edge_term * mollifier_val * vert_term * barrier;
    }
}