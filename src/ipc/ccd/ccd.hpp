#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// The default tolerance used with Tight-Inclusion CCD.
static constexpr double DEFAULT_CCD_TOLERANCE = 1e-6;
/// The default maximum number of iterations used with Tight-Inclusion CCD.
static constexpr long DEFAULT_CCD_MAX_ITERATIONS = 10'000'000l;
/// The default conservative rescaling value used to avoid taking steps exactly
/// to impact.
static constexpr double DEFAULT_CCD_CONSERVATIVE_RESCALING = 0.8;

// 2D

bool point_edge_ccd_2D(
    const Eigen::Vector2d& p_t0,
    const Eigen::Vector2d& e0_t0,
    const Eigen::Vector2d& e1_t0,
    const Eigen::Vector2d& p_t1,
    const Eigen::Vector2d& e0_t1,
    const Eigen::Vector2d& e1_t1,
    double& toi,
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

// 3D

bool point_point_ccd(
    const Eigen::Vector3d& p0_t0,
    const Eigen::Vector3d& p1_t0,
    const Eigen::Vector3d& p0_t1,
    const Eigen::Vector3d& p1_t1,
    double& toi,
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

bool point_edge_ccd_3D(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& e0_t0,
    const Eigen::Vector3d& e1_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& e0_t1,
    const Eigen::Vector3d& e1_t1,
    double& toi,
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

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
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

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
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

// 2D or 3D

bool point_edge_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double& toi,
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

} // namespace ipc
