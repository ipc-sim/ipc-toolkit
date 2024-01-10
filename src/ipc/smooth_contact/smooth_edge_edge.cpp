#include "smooth_edge_edge.hpp"

namespace ipc {
    template <typename scalar>
    scalar smooth_point_edge_potential_single_point_3d(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const double &uv,
        const ParameterType &params)
    {
        const Vector3<scalar> tangent = e1 - e0;
        const Vector3<scalar> q = e0 + tangent * uv;
        const Vector3<scalar> tmp = f0 - q;
        const Vector3<scalar> t0 = tmp - (tmp.dot(tangent) / tangent.squaredNorm()) * tangent;
    }
}