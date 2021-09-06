#pragma once

#include <Eigen/Core>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// 2D

bool point_edge_ccd_2D(
    const Eigen::Vector2d& p_t0,
    const Eigen::Vector2d& e0_t0,
    const Eigen::Vector2d& e1_t0,
    const Eigen::Vector2d& p_t1,
    const Eigen::Vector2d& e0_t1,
    const Eigen::Vector2d& e1_t1,
    double& toi,
    double tmax = 1.0,
    double tolerance = 1e-6,
    int max_iterations = 1e7,
    double conservative_rescaling = 0.8);

// 3D

bool point_point_ccd(
    const Eigen::Vector3d& p0_t0,
    const Eigen::Vector3d& p1_t0,
    const Eigen::Vector3d& p0_t1,
    const Eigen::Vector3d& p1_t1,
    double& toi,
    double tmax = 1.0,
    double tolerance = 1e-6,
    int max_iterations = 1e7,
    double conservative_rescaling = 0.8);

bool point_edge_ccd_3D(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& e0_t0,
    const Eigen::Vector3d& e1_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& e0_t1,
    const Eigen::Vector3d& e1_t1,
    double& toi,
    double tmax = 1.0,
    double tolerance = 1e-6,
    int max_iterations = 1e7,
    double conservative_rescaling = 0.8);

bool point_triangle_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double& toi,
    double tmax = 1.0,
    double tolerance = 1e-6,
    int max_iterations = 1e7,
    double conservative_rescaling = 0.8);

bool edge_edge_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double& toi,
    double tmax = 1.0,
    double tolerance = 1e-6,
    int max_iterations = 1e7,
    double conservative_rescaling = 0.8);

// 2D or 3D

bool point_edge_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double& toi,
    double tmax = 1.0,
    double tolerance = 1e-6,
    int max_iterations = 1e7,
    double conservative_rescaling = 0.8);

bool point_static_plane_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& plane_origin,
    const Eigen::Vector3d& plane_normal,
    double& toi,
    double conservative_rescaling = 0.8);

} // namespace ipc
