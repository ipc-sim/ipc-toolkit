#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ipc {

/// @brief Compute the distance between a point and a plane (defined by a
/// triangle).
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0,t1,t2 The points of the triangle defining the plane.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_plane_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    auto normal = (t1 - t0).cross(t2 - t0);
    auto point_to_plane = (p - t0).dot(normal);
    return point_to_plane * point_to_plane / normal.squaredNorm();
}

// Symbolically generated derivatives;
namespace autogen {
    void point_plane_distance_gradient(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double g[12]);

    void point_plane_distance_hessian(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double H[144]);
} // namespace autogen

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_plane_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    Eigen::VectorXd& grad)
{
    grad.resize(p.size() + t0.size() + t1.size() + t2.size());
    autogen::point_plane_distance_gradient(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], grad.data());
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_plane_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    Eigen::MatrixXd& hess)
{
    hess.resize(
        p.size() + t0.size() + t1.size() + t2.size(),
        p.size() + t0.size() + t1.size() + t2.size());
    autogen::point_plane_distance_hessian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], hess.data());
}

} // namespace ipc
