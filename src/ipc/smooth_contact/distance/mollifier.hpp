#pragma once
#include "point_face.hpp"
#include "point_edge.hpp"
#include "edge_edge.hpp"

namespace ipc {
template <typename scalar>
scalar point_edge_mollifier(
    const VectorMax3<scalar>& p,
    const VectorMax3<scalar>& e0,
    const VectorMax3<scalar>& e1,
    const scalar& dist_sqr);

inline std::array<HEAVISIDE_TYPE, 4> edge_edge_mollifier_type(
    const Vector3<double>& ea0,
    const Vector3<double>& ea1,
    const Vector3<double>& eb0,
    const Vector3<double>& eb1,
    const double& dist_sqr);

template <typename scalar>
scalar edge_edge_mollifier(
    const Vector3<scalar>& ea0,
    const Vector3<scalar>& ea1,
    const Vector3<scalar>& eb0,
    const Vector3<scalar>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const scalar& dist_sqr);

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
std::pair<double, Vector<double, 13>> edge_edge_mollifier_grad(
    const Vector3<double>& ea0,
    const Vector3<double>& ea1,
    const Vector3<double>& eb0,
    const Vector3<double>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr);

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
std::tuple<double, Vector<double, 13>, Eigen::Matrix<double, 13, 13>>
edge_edge_mollifier_hessian(
    const Vector3<double>& ea0,
    const Vector3<double>& ea1,
    const Vector3<double>& eb0,
    const Vector3<double>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr);

template <typename scalar>
scalar point_face_mollifier(
    const VectorMax3<scalar>& p,
    const VectorMax3<scalar>& e0,
    const VectorMax3<scalar>& e1,
    const VectorMax3<scalar>& e2,
    const scalar& dist_sqr);
} // namespace ipc

#include "mollifier.tpp"
