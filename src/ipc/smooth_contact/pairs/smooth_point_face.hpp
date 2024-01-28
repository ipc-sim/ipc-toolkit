#pragma once

#include <ipc/smooth_contact/primitives/point.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>
#include <ipc/smooth_contact/common.hpp>

namespace ipc {
    template <typename scalar>
    scalar smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Matrix<scalar, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype,
        const ORIENTATION_TYPES &otypes)
    {
        Vector3<scalar> direc = point_triangle_closest_point_direction<scalar>(p, v0, v1, v2, dtype);
        const scalar dist_sqr = direc.squaredNorm();

        auto b = inv_barrier(sqrt(dist_sqr) / params.dhat, params.r);
        auto ff = smooth_face_term<scalar>(p, v0, v1, v2);
        auto pp = smooth_point3_term<scalar>(p, direc / direc.norm(), neighbors, params.alpha, params.beta, otypes);
        auto tt = triangle_mollifier<scalar>(p - direc, v0, v1, v2, dist_sqr);

        // if constexpr (std::is_same<double,scalar>::value)
        // {
        //     if (dist_sqr < 1e-20)
        //         logger().error("barrier {}, face {}, point {}, mollifier {}", b, ff, pp, tt);
        // }

        return b * ff * pp * tt;
    }
}