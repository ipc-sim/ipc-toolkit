#pragma once

#include <ipc/utils/distance_autodiff.hpp>

namespace ipc {
    /// @brief Compute potential for a point p and an face [v0, v1, v2], using the smooth closest point
    // template <typename scalar>
    // scalar smooth_face_term(
    //     const Eigen::Ref<const Vector3<scalar>>& p,
    //     const Eigen::Ref<const Vector3<scalar>>& v0,
    //     const Eigen::Ref<const Vector3<scalar>>& v1,
    //     const Eigen::Ref<const Vector3<scalar>>& v2,
    //     const scalar &dist_sqr,
    //     const double &alpha)
    // {
    //     const Vector3<scalar> normal = (v1 - v0).cross(v2 - v0);
    //     const scalar sqr_norm = normal.squaredNorm();
    //     const scalar Phi = 1 - (p - v0).dot(normal) / sqrt(dist_sqr * sqr_norm);
    //     const scalar tri_area = 0.5 * sqr_norm;

    //     return tri_area * cubic_spline(Phi / alpha);
    // }

    /// @brief d points from triangle to the point
    template <typename scalar>
    scalar smooth_face_term(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2)
    {
        const Vector3<scalar> normal = (v1 - v0).cross(v2 - v0);

        if (normal.dot(p - v0) < 0)
            return scalar(0.);
        else
            return 0.5 * normal.norm(); // area of triangle
    }
}