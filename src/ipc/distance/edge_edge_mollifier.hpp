#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the squared norm of the edge-edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The squared norm of the edge-edge cross product.
double edge_edge_cross_squarednorm(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the gradient of the squared norm of the edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The gradient of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
Vector12d edge_edge_cross_squarednorm_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the hessian of the squared norm of the edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The hessian of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
Matrix12d edge_edge_cross_squarednorm_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The mollifier coefficient to premultiply the edge-edge distance.
double edge_edge_mollifier(const double x, const double eps_x);

/// @brief The gradient of the mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The gradient of the mollifier function for edge-edge distance wrt x.
double edge_edge_mollifier_gradient(const double x, const double eps_x);

/// @brief The derivative of the mollifier function for edge-edge distance wrt eps_x.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The derivative of the mollifier function for edge-edge distance wrt eps_x.
double
edge_edge_mollifier_derivative_wrt_eps_x(const double x, const double eps_x);

/// @brief The hessian of the mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The hessian of the mollifier function for edge-edge distance wrt x.
double edge_edge_mollifier_hessian(const double x, const double eps_x);

/// @brief The derivative of the gradient of the mollifier function for edge-edge distance wrt eps_x.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The derivative of the gradient of the mollifier function for edge-edge distance wrt eps_x.
double edge_edge_mollifier_gradient_derivative_wrt_eps_x(
    const double x, const double eps_x);

/// @brief Compute a mollifier for the edge-edge distance.
///
/// This helps smooth the non-smoothness at close to parallel edges.
///
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param eps_x Mollifier activation threshold.
/// @return The mollifier coefficient to premultiply the edge-edge distance.
double edge_edge_mollifier(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    const double eps_x);

/// @brief Compute the gradient of the mollifier for the edge-edge distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param eps_x Mollifier activation threshold.
/// @return The gradient of the mollifier.
Vector12d edge_edge_mollifier_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    const double eps_x);

/// @brief Compute the hessian of the mollifier for the edge-edge distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param eps_x Mollifier activation threshold.
/// @return The hessian of the mollifier.
Matrix12d edge_edge_mollifier_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    const double eps_x);

/// @brief Compute the gradient of the mollifier for the edge-edge distance wrt rest positions.
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The derivative of the mollifier wrt rest positions.
Vector12d edge_edge_mollifier_gradient_wrt_x(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the jacobian of the edge-edge distance mollifier's gradient wrt rest positions.
/// @note This is not the hessian of the mollifier wrt rest positions, but the jacobian wrt rest positions of the mollifier's gradient wrt positions.
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The jacobian of the mollifier's gradient wrt rest positions.
Matrix12d edge_edge_mollifier_gradient_jacobian_wrt_x(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the threshold of the mollifier edge-edge distance.
///
/// This values is computed based on the edges at rest length.
///
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @return Threshold for edge-edge mollification.
double edge_edge_mollifier_threshold(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest);

/// @brief Compute the gradient of the threshold of the mollifier edge-edge distance.
///
/// This values is computed based on the edges at rest length.
///
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @return Gradient of the threshold for edge-edge mollification.
Vector12d edge_edge_mollifier_threshold_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest);

// Symbolically generated derivatives;
namespace autogen {
    void edge_edge_cross_squarednorm_gradient(
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

    void edge_edge_cross_squarednorm_hessian(
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

    void edge_edge_mollifier_threshold_gradient(
        double ea0x,
        double ea0y,
        double ea0z,
        double ea1x,
        double ea1y,
        double ea1z,
        double eb0x,
        double eb0y,
        double eb0z,
        double eb1x,
        double eb1y,
        double eb1z,
        double grad[12],
        double scale = 1e-3);
} // namespace autogen
} // namespace ipc
