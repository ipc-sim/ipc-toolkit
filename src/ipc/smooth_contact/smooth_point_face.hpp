#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>
#include "smooth_point_point.hpp"
#include <ipc/utils/distance_autodiff.hpp>

namespace ipc {

    /// @brief Compute pointwise potential for a point p and a point specified by uv on face [v0, v1, v2]
    /// @param p One point outside of face
    /// @param uv Barycentric coordinate
    /// @param dhat The effective distance of barrier
    /// @param alpha The effective angle of barrier
    template <typename scalar>
    scalar smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params);

    /// @brief Compute potential for a point p and an face [v0, v1, v2], integrated over the face with high order quadrature
    template <typename scalar>
    scalar smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const ParameterType &params);

    /// @brief Compute potential for a list of points and an face [v0, v1, v2], integrated over the face with high order quadrature
    template <typename scalar>
    Vector<scalar, -1, -1> smooth_point_face_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<scalar, -1, 3>>& points,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const ParameterType &params);

    /// @brief Compute potential for a point p and an face [v0, v1, v2], using the smooth closest point
    template <typename scalar>
    scalar smooth_face_term(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const scalar &dist_sqr,
        const double &alpha)
    {
        const Vector3<scalar> normal = (v1 - v0).cross(v2 - v0);
        const scalar sqr_norm = normal.squaredNorm();
        const scalar Phi = 1 - (p - v0).dot(normal) / sqrt(dist_sqr * sqr_norm);

        return 0.5 * sqr_norm * cubic_spline(Phi / alpha);
    }

    template <typename scalar>
    scalar smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Matrix<scalar, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype)
    {
        Vector3<scalar> direc = point_triangle_closest_point_direction<scalar>(p, v0, v1, v2, dtype);
        const scalar dist_sqr = direc.squaredNorm();

        auto b = inv_barrier(dist_sqr / params.eps, params.r);
        auto ff = smooth_face_term<scalar>(p, v0, v1, v2, dist_sqr, params.alpha);
        auto pp = smooth_point3_term<scalar>(p, direc / direc.norm(), neighbors, params.alpha, params.beta);
        auto tt = triangle_mollifier<scalar>(p - direc, v0, v1, v2, dist_sqr);

        // if constexpr (std::is_same<double,scalar>::value)
        // {
        //     if (dist_sqr < 1e-20)
        //         logger().error("barrier {}, face {}, point {}, mollifier {}", b, ff, pp, tt);
        // }

        return b * ff * pp * tt;
    }
}