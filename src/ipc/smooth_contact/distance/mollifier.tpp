#pragma once

#include "mollifier.hpp"

namespace ipc {
    template <typename scalar>
    scalar point_edge_mollifier(const VectorMax3<scalar> &p, const VectorMax3<scalar> &e0, const VectorMax3<scalar> &e1, const scalar &dist_sqr)
    {
        const scalar denominator = (e1 - e0).squaredNorm() * mollifier_threshold_eps;
        return Math<scalar>::mollifier(((p - e0).squaredNorm() - dist_sqr) / denominator) *
            Math<scalar>::mollifier(((p - e1).squaredNorm() - dist_sqr) / denominator);
    }

    inline std::array<HEAVISIDE_TYPE, 4> edge_edge_mollifier_type(
        const Vector3<double> &ea0, const Vector3<double> &ea1,
        const Vector3<double> &eb0, const Vector3<double> &eb1,
        const double &dist_sqr)
    {
        std::array<HEAVISIDE_TYPE, 4> mtypes;
        mtypes[0] = (point_edge_distance(ea0, eb0, eb1) - dist_sqr) >= (eb1 - eb0).squaredNorm() * mollifier_threshold_eps ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        mtypes[1] = (point_edge_distance(ea1, eb0, eb1) - dist_sqr) >= (eb1 - eb0).squaredNorm() * mollifier_threshold_eps ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        mtypes[2] = (point_edge_distance(eb0, ea0, ea1) - dist_sqr) >= (ea1 - ea0).squaredNorm() * mollifier_threshold_eps ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        mtypes[3] = (point_edge_distance(eb1, ea0, ea1) - dist_sqr) >= (ea1 - ea0).squaredNorm() * mollifier_threshold_eps ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        return mtypes;
    }

    template <typename scalar>
    scalar edge_edge_mollifier(
        const Vector3<scalar> &ea0, const Vector3<scalar> &ea1,
        const Vector3<scalar> &eb0, const Vector3<scalar> &eb1, 
        const std::array<HEAVISIDE_TYPE, 4> &mtypes,
        const scalar &dist_sqr)
    {
        const scalar da = (ea1 - ea0).squaredNorm() * mollifier_threshold_eps;
        const scalar db = (eb1 - eb0).squaredNorm() * mollifier_threshold_eps;
        scalar a = (mtypes[0] == HEAVISIDE_TYPE::VARIANT) ? Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(ea0, eb0, eb1) - dist_sqr) / db) : scalar(1.);
        scalar b = (mtypes[1] == HEAVISIDE_TYPE::VARIANT) ? Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(ea1, eb0, eb1) - dist_sqr) / db) : scalar(1.);
        scalar c = (mtypes[2] == HEAVISIDE_TYPE::VARIANT) ? Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(eb0, ea0, ea1) - dist_sqr) / da) : scalar(1.);
        scalar d = (mtypes[3] == HEAVISIDE_TYPE::VARIANT) ? Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(eb1, ea0, ea1) - dist_sqr) / da) : scalar(1.);
        
        return a * b * c * d;
    }

    template <typename scalar>
    scalar point_face_mollifier(
        const VectorMax3<scalar> &p, 
        const VectorMax3<scalar> &e0,
        const VectorMax3<scalar> &e1,
        const VectorMax3<scalar> &e2,
        const scalar &dist_sqr)
    {
        return Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(p, e0, e1) - dist_sqr) / mollifier_threshold_eps / PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(e2, e0, e1)) *
            Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(p, e2, e1) - dist_sqr) / mollifier_threshold_eps / PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(e0, e2, e1)) *
            Math<scalar>::mollifier((PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(p, e0, e2) - dist_sqr) / mollifier_threshold_eps / PointEdgeDistance<scalar, 3>::point_edge_sqr_distance(e1, e0, e2));
    }
}