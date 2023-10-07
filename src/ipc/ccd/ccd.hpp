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

/// @brief Computes the time of impact between a point and an edge in 2D using continuous collision detection.
/// @param p_t0 The initial position of the point.
/// @param e0_t0 The initial position of the first endpoint of the edge.
/// @param e1_t0 The initial position of the second endpoint of the edge.
/// @param p_t1 The final position of the point.
/// @param e0_t1 The final position of the first endpoint of the edge.
/// @param e1_t1 The final position of the second endpoint of the edge.
/// @param[out] toi The time of impact between the point and the edge.
/// @param min_distance The minimum distance between the objects.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
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

/// @brief Computes the time of impact between two points in 3D using continuous collision detection.
/// @param p0_t0 The initial position of the first point.
/// @param p1_t0 The initial position of the second point.
/// @param p0_t1 The final position of the first point.
/// @param p1_t1 The final position of the second point.
/// @param[out] toi The time of impact between the two points.
/// @param min_distance The minimum distance between the objects.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
bool point_point_ccd_3D(
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

/// @brief Computes the time of impact between a point and an edge in 3D using continuous collision detection.
/// @param p_t0 The initial position of the point.
/// @param e0_t0 The initial position of the first endpoint of the edge.
/// @param e1_t0 The initial position of the second endpoint of the edge.
/// @param p_t1 The final position of the point.
/// @param e0_t1 The final position of the first endpoint of the edge.
/// @param e1_t1 The final position of the second endpoint of the edge.
/// @param[out] toi The time of impact between the point and the edge.
/// @param min_distance The minimum distance between the objects.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
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

/// @brief Computes the time of impact between a point and a triangle in 3D using continuous collision detection.
/// @param p_t0 The initial position of the point.
/// @param t0_t0 The initial position of the first vertex of the triangle.
/// @param t1_t0 The initial position of the second vertex of the triangle.
/// @param t2_t0 The initial position of the third vertex of the triangle.
/// @param p_t1 The final position of the point.
/// @param t0_t1 The final position of the first vertex of the triangle.
/// @param t1_t1 The final position of the second vertex of the triangle.
/// @param t2_t1 The final position of the third vertex of the triangle.
/// @param[out] toi The time of impact between the point and the triangle.
/// @param min_distance The minimum distance between the objects.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
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

/// @brief Computes the time of impact between two edges in 3D using continuous collision detection.
/// @param ea0_t0 The initial position of the first endpoint of the first edge.
/// @param ea1_t0 The initial position of the second endpoint of the first edge.
/// @param eb0_t0 The initial position of the first endpoint of the second edge.
/// @param eb1_t0 The initial position of the second endpoint of the second edge.
/// @param ea0_t1 The final position of the first endpoint of the first edge.
/// @param ea1_t1 The final position of the second endpoint of the first edge.
/// @param eb0_t1 The final position of the first endpoint of the second edge.
/// @param eb1_t1 The final position of the second endpoint of the second edge.
/// @param[out] toi The time of impact between the two edges.
/// @param min_distance The minimum distance between the objects.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
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

/// @brief Computes the time of impact between two points using continuous collision detection.
/// @param p0_t0 The initial position of the first point.
/// @param p1_t0 The initial position of the second point.
/// @param p0_t1 The final position of the first point.
/// @param p1_t1 The final position of the second point.
/// @param[out] toi The time of impact between the two points.
/// @param min_distance The minimum distance between the points.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
bool point_point_ccd(
    const VectorMax3d& p0_t0,
    const VectorMax3d& p1_t0,
    const VectorMax3d& p0_t1,
    const VectorMax3d& p1_t1,
    double& toi,
    const double min_distance = 0.0,
    const double tmax = 1.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

/// @brief Computes the time of impact between a point and an edge in 2D or 3D using continuous collision detection.
/// @param p_t0 The initial position of the point.
/// @param e0_t0 The initial position of the first endpoint of the edge.
/// @param e1_t0 The initial position of the second endpoint of the edge.
/// @param p_t1 The final position of the point.
/// @param e0_t1 The final position of the first endpoint of the edge.
/// @param e1_t1 The final position of the second endpoint of the edge.
/// @param[out] toi The time of impact between the point and the edge.
/// @param min_distance The minimum distance between the objects.
/// @param tmax The maximum time to check for collisions.
/// @param tolerance The error tolerance for the time of impact.
/// @param max_iterations The maximum number of iterations to perform.
/// @param conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
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
