#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Discrete collision detection

bool point_edge_aabb_cd(
    const VectorMax3d& p,
    const VectorMax3d& e0,
    const VectorMax3d& e1,
    double dist);

bool edge_edge_aabb_cd(
    const VectorMax3d& ea0,
    const VectorMax3d& ea1,
    const VectorMax3d& eb0,
    const VectorMax3d& eb1,
    double dist);

bool point_triangle_aabb_cd(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double dist);

bool edge_triangle_aabb_cd(
    const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double dist);

// Continous collision detection

bool point_edge_aabb_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double dist);

bool edge_edge_aabb_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double dist);

bool point_triangle_aabb_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double dist);

} // namespace ipc
