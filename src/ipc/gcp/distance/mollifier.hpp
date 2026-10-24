#pragma once
#include "edge_edge.hpp"
#include "point_edge.hpp"
#include "point_face.hpp"

namespace ipc {
template <typename scalar, int dim>
scalar point_edge_mollifier(
    Eigen::ConstRef<Eigen::Vector<scalar, dim>> p,
    Eigen::ConstRef<Eigen::Vector<scalar, dim>> e0,
    Eigen::ConstRef<Eigen::Vector<scalar, dim>> e1,
    const scalar& dist_sqr);

std::array<HeavisideType, 4> edge_edge_mollifier_type(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    const double dist_sqr);

template <typename scalar>
scalar edge_edge_mollifier(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1,
    const std::array<HeavisideType, 4>& mtypes,
    const scalar& dist_sqr);

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
GradientType<13> edge_edge_mollifier_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    const std::array<HeavisideType, 4>& mtypes,
    const double dist_sqr);

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
HessianType<13> edge_edge_mollifier_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    const std::array<HeavisideType, 4>& mtypes,
    const double dist_sqr);

template <typename scalar>
scalar point_face_mollifier(
    Eigen::ConstRef<Eigen::Vector3<scalar>> p,
    Eigen::ConstRef<Eigen::Vector3<scalar>> e0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> e1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> e2,
    const scalar& dist_sqr);

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
GradientType<13> point_face_mollifier_gradient(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> e2,
    const double dist_sqr);

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
HessianType<13> point_face_mollifier_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> e2,
    const double dist_sqr);
} // namespace ipc

#include "mollifier.tpp"
