#pragma once

#include "mollifier.hpp"

#include <ipc/distance/point_line.hpp>
#include <ipc/gcp/distance/mollifier_config.hpp>

namespace ipc {
template <typename scalar, int dim>
scalar point_edge_mollifier(
    Eigen::ConstRef<Eigen::Vector<scalar, dim>> p,
    Eigen::ConstRef<Eigen::Vector<scalar, dim>> e0,
    Eigen::ConstRef<Eigen::Vector<scalar, dim>> e1,
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
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1,
    const std::array<HeavisideType, 4>& mtypes,
    const scalar& dist_sqr)
{
    if (ee_mollifier_config().smoothstep) {
        // PROTOTYPE: H(c / l) with H(t) = 10t^3 - 15t^4 + 6t^5 (C2), where
        // c = sqrt(d_PE^2 - d^2) is the endpoint clearance and l is fixed.
        using std::sqrt;
        const scalar ell2 = ee_mollifier_config().rel_width
            ? scalar(dist_sqr * ee_mollifier_config().rel_eps)
            : ee_mollifier_config().sqrt_width
            ? scalar(sqrt(dist_sqr)
                     * (MOLLIFIER_THRESHOLD_EPS * ee_mollifier_config().L))
            : scalar(ee_mollifier_config().ell * ee_mollifier_config().ell);
        const auto H = [&](const scalar& c2) -> scalar {
            if (c2 <= 0) {
                return scalar(0.);
            }
            if (c2 >= ell2) {
                return scalar(1.);
            }
            const scalar t = sqrt(c2 / ell2);
            return t * t * t * (10. + t * (-15. + 6. * t));
        };
        using PE = PointEdgeDistance<scalar, 3>;
        scalar w(1.);
        if (mtypes[0] == HeavisideType::VARIANT)
            w = w * H(PE::point_edge_sqr_distance(ea0, eb0, eb1) - dist_sqr);
        if (mtypes[1] == HeavisideType::VARIANT)
            w = w * H(PE::point_edge_sqr_distance(ea1, eb0, eb1) - dist_sqr);
        if (mtypes[2] == HeavisideType::VARIANT)
            w = w * H(PE::point_edge_sqr_distance(eb0, ea0, ea1) - dist_sqr);
        if (mtypes[3] == HeavisideType::VARIANT)
            w = w * H(PE::point_edge_sqr_distance(eb1, ea0, ea1) - dist_sqr);
        return w;
    }
    scalar denom = dist_sqr * MOLLIFIER_THRESHOLD_EPS;
    if (ee_mollifier_config().fixed) { // PROTOTYPE
        denom = scalar(ee_mollifier_config().ell * ee_mollifier_config().ell);
    } else if (ee_mollifier_config().sqrt_width) { // PROTOTYPE
        using std::sqrt;
        denom = sqrt(dist_sqr)
            * (MOLLIFIER_THRESHOLD_EPS * ee_mollifier_config().L);
    }
    const scalar da = denom;
    const scalar db = denom;
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
    Eigen::ConstRef<Eigen::Vector3<scalar>> p,
    Eigen::ConstRef<Eigen::Vector3<scalar>> e0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> e1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> e2,
    const scalar& dist_sqr)
{
    // use point-line distance instead of point-edge distance because
    // this function vanishes if the point is outside the triangle, so
    // whenever this function is nonzero the point-edge distance equals
    // the point-line distance
    return Math<scalar>::mollifier(
               (point_line_distance(p, e0, e1) - dist_sqr)
               / MOLLIFIER_THRESHOLD_EPS / dist_sqr)
        * Math<scalar>::mollifier(
               (point_line_distance(p, e2, e1) - dist_sqr)
               / MOLLIFIER_THRESHOLD_EPS / dist_sqr)
        * Math<scalar>::mollifier(
               (point_line_distance(p, e0, e2) - dist_sqr)
               / MOLLIFIER_THRESHOLD_EPS / dist_sqr);
}
} // namespace ipc