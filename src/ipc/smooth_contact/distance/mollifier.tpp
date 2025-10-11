#pragma once

#include "mollifier.hpp"

namespace ipc {
template <typename scalar, int dim>
scalar point_edge_mollifier(
    Eigen::ConstRef<Vector<scalar, dim>> p,
    Eigen::ConstRef<Vector<scalar, dim>> e0,
    Eigen::ConstRef<Vector<scalar, dim>> e1,
    const scalar& dist_sqr)
{
    const scalar denominator = dist_sqr * MOLLIFIER_THRESHOLD_EPS;
    return Math<scalar>::mollifier(
               ((p - e0).squaredNorm() - dist_sqr) / denominator)
        * Math<scalar>::mollifier(
               ((p - e1).squaredNorm() - dist_sqr) / denominator);
}

template <typename scalar>
scalar edge_edge_mollifier(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1,
    const std::array<HeavisideType, 4>& mtypes,
    const scalar& dist_sqr)
{
    const scalar da = dist_sqr * MOLLIFIER_THRESHOLD_EPS;
    const scalar db = dist_sqr * MOLLIFIER_THRESHOLD_EPS;
    scalar a = (mtypes[0] == HeavisideType::VARIANT)
        ? Math<scalar>::mollifier(
              (PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(
                   ea0, eb0, eb1)
               - dist_sqr)
              / db)
        : scalar(1.);
    scalar b = (mtypes[1] == HeavisideType::VARIANT)
        ? Math<scalar>::mollifier(
              (PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(
                   ea1, eb0, eb1)
               - dist_sqr)
              / db)
        : scalar(1.);
    scalar c = (mtypes[2] == HeavisideType::VARIANT)
        ? Math<scalar>::mollifier(
              (PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(
                   eb0, ea0, ea1)
               - dist_sqr)
              / da)
        : scalar(1.);
    scalar d = (mtypes[3] == HeavisideType::VARIANT)
        ? Math<scalar>::mollifier(
              (PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(
                   eb1, ea0, ea1)
               - dist_sqr)
              / da)
        : scalar(1.);

    return a * b * c * d;
}

template <typename scalar>
scalar point_face_mollifier(
    Eigen::ConstRef<Vector3<scalar>> p,
    Eigen::ConstRef<Vector3<scalar>> e0,
    Eigen::ConstRef<Vector3<scalar>> e1,
    Eigen::ConstRef<Vector3<scalar>> e2,
    const scalar& dist_sqr)
{
    // use point-line distance instead of point-edge distance because
    // this function vanishes if the point is outside the triangle, so
    // whenever this function is nonzero the point-edge distance equals
    // the point-line distance
    return Math<scalar>::mollifier(
               (PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, e0, e1)
                - dist_sqr)
               / MOLLIFIER_THRESHOLD_EPS / dist_sqr)
        * Math<scalar>::mollifier(
               (PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, e2, e1)
                - dist_sqr)
               / MOLLIFIER_THRESHOLD_EPS / dist_sqr)
        * Math<scalar>::mollifier(
               (PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, e0, e2)
                - dist_sqr)
               / MOLLIFIER_THRESHOLD_EPS / dist_sqr);
}
} // namespace ipc