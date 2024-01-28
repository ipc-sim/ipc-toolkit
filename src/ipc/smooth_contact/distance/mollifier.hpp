#pragma once
#include "point_face.hpp"
#include "point_edge.hpp"
#include "edge_edge.hpp"

namespace ipc {
    template <typename scalar>
    scalar edge_mollifier(const VectorMax3<scalar> &p, const VectorMax3<scalar> &e0, const VectorMax3<scalar> &e1, const scalar &dist);

    inline std::array<HEAVISIDE_TYPE, 4> edge_edge_mollifier_type(
        const Vector3<double> &ea0, const Vector3<double> &ea1,
        const Vector3<double> &eb0, const Vector3<double> &eb1, 
        const double &dist);

    template <typename scalar>
    scalar edge_edge_mollifier(
        const Vector3<scalar> &ea0, const Vector3<scalar> &ea1,
        const Vector3<scalar> &eb0, const Vector3<scalar> &eb1, 
        const std::array<HEAVISIDE_TYPE, 4> mtypes,
        const scalar &dist);

    template <typename scalar>
    scalar triangle_mollifier(
        const VectorMax3<scalar> &p, 
        const VectorMax3<scalar> &e0, 
        const VectorMax3<scalar> &e1,
        const VectorMax3<scalar> &e2,
        const scalar &dist);
}

#include "mollifier.tpp"
