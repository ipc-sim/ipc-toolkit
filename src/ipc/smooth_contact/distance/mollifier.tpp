#pragma once

#include "mollifier.hpp"

namespace ipc {
    template <typename scalar>
    scalar edge_mollifier(const VectorMax3<scalar> &p, const VectorMax3<scalar> &e0, const VectorMax3<scalar> &e1, const scalar &dist)
    {
        const scalar denominator = dist * mollifier_threshold_eps;
        return mollifier<scalar>(((p - e0).norm() - dist) / denominator) *
            mollifier<scalar>(((p - e1).norm() - dist) / denominator);
    }

    inline std::array<HEAVISIDE_TYPE, 4> edge_edge_mollifier_type(
        const Vector3<double> &ea0, const Vector3<double> &ea1,
        const Vector3<double> &eb0, const Vector3<double> &eb1, 
        const double &dist)
    {
        const double dist_sqr = dist * dist;
        const double denominator = dist_sqr * mollifier_threshold_eps;
        std::array<HEAVISIDE_TYPE, 4> mtypes;
        mtypes[0] = (point_edge_distance(ea0, eb0, eb1) - dist_sqr) >= denominator ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        mtypes[1] = (point_edge_distance(ea1, eb0, eb1) - dist_sqr) >= denominator ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        mtypes[2] = (point_edge_distance(eb0, ea0, ea1) - dist_sqr) >= denominator ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        mtypes[3] = (point_edge_distance(eb1, ea0, ea1) - dist_sqr) >= denominator ? HEAVISIDE_TYPE::ONE : HEAVISIDE_TYPE::VARIANT;
        return mtypes;
    }

    template <typename scalar>
    scalar edge_edge_mollifier(
        const Vector3<scalar> &ea0, const Vector3<scalar> &ea1,
        const Vector3<scalar> &eb0, const Vector3<scalar> &eb1, 
        const std::array<HEAVISIDE_TYPE, 4> mtypes,
        const scalar &dist)
    {
        const scalar dist_sqr = dist * dist;
        const scalar denominator = dist_sqr * mollifier_threshold_eps;
        scalar a = mtypes[0] == HEAVISIDE_TYPE::VARIANT ? mollifier<scalar>((point_edge_sqr_distance<scalar>(ea0, eb0, eb1) - dist_sqr) / denominator) : (mtypes[0] == HEAVISIDE_TYPE::ZERO ? scalar(0.) : scalar(1.));
        scalar b = mtypes[1] == HEAVISIDE_TYPE::VARIANT ? mollifier<scalar>((point_edge_sqr_distance<scalar>(ea1, eb0, eb1) - dist_sqr) / denominator) : (mtypes[1] == HEAVISIDE_TYPE::ZERO ? scalar(0.) : scalar(1.));
        scalar c = mtypes[2] == HEAVISIDE_TYPE::VARIANT ? mollifier<scalar>((point_edge_sqr_distance<scalar>(eb0, ea0, ea1) - dist_sqr) / denominator) : (mtypes[2] == HEAVISIDE_TYPE::ZERO ? scalar(0.) : scalar(1.));
        scalar d = mtypes[3] == HEAVISIDE_TYPE::VARIANT ? mollifier<scalar>((point_edge_sqr_distance<scalar>(eb1, ea0, ea1) - dist_sqr) / denominator) : (mtypes[3] == HEAVISIDE_TYPE::ZERO ? scalar(0.) : scalar(1.));
        
        return a * b * c * d;
    }

    template <typename scalar>
    scalar triangle_mollifier(
        const VectorMax3<scalar> &p, 
        const VectorMax3<scalar> &e0,
        const VectorMax3<scalar> &e1,
        const VectorMax3<scalar> &e2,
        const scalar &dist)
    {
        const scalar denominator = dist*dist * mollifier_threshold_eps;
        return mollifier<scalar>(point_edge_sqr_distance<scalar>(p, e0, e1) / denominator) *
            mollifier<scalar>(point_edge_sqr_distance<scalar>(p, e2, e1) / denominator) *
            mollifier<scalar>(point_edge_sqr_distance<scalar>(p, e0, e2) / denominator);
    }
}