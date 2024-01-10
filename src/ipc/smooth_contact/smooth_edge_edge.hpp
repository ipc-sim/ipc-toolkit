#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>

namespace ipc {
    /// @param p is one point outside of edge
    /// @param e0, e1 form one edge
    /// @param (e0, e1, f0) and (e0, e1, f1) are two triangle faces (orientation doesn't matter)
    /// @param uv is the barycentric coordinate on (e0, e1) that is closest to p
    template <typename scalar>
    scalar smooth_point_edge_potential_single_point_3d(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const double &uv,
        const ParameterType &params);
}