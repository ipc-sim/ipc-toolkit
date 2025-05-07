#pragma once
#include "point_face.hpp"
#include "point_edge.hpp"
#include "edge_edge.hpp"

namespace ipc {
template <typename scalar, int dim>
scalar point_edge_mollifier(
    const Eigen::Ref<const Vector<scalar, dim>>& p,
    const Eigen::Ref<const Vector<scalar, dim>>& e0,
    const Eigen::Ref<const Vector<scalar, dim>>& e1,
    const scalar& dist_sqr);

std::array<HEAVISIDE_TYPE, 4> edge_edge_mollifier_type(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1,
    const double& dist_sqr);

template <typename scalar>
scalar edge_edge_mollifier(
    const Eigen::Ref<const Vector3<scalar>>& ea0,
    const Eigen::Ref<const Vector3<scalar>>& ea1,
    const Eigen::Ref<const Vector3<scalar>>& eb0,
    const Eigen::Ref<const Vector3<scalar>>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const scalar& dist_sqr);

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
GradType<13> edge_edge_mollifier_gradient(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr);

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
HessianType<13>
edge_edge_mollifier_hessian(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr);

template <typename scalar>
scalar point_face_mollifier(
    const Eigen::Ref<const Vector3<scalar>>& p,
    const Eigen::Ref<const Vector3<scalar>>& e0,
    const Eigen::Ref<const Vector3<scalar>>& e1,
    const Eigen::Ref<const Vector3<scalar>>& e2,
    const scalar& dist_sqr);

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
GradType<13>
point_face_mollifier_gradient(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& e2,
    const double& dist_sqr);

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
HessianType<13>
point_face_mollifier_hessian(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& e2,
    const double& dist_sqr);
} // namespace ipc

#include "mollifier.tpp"
