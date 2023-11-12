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
/// Tolerance for small time of impact which triggers rerunning CCD without a
/// minimum separation.
static constexpr double CCD_SMALL_TOI = 1e-6;

/// @brief Computes the time of impact between two points in 3D using continuous collision detection.
/// @param[in] p0_t0 The initial position of the first point.
/// @param[in] p1_t0 The initial position of the second point.
/// @param[in] p0_t1 The final position of the first point.
/// @param[in] p1_t1 The final position of the second point.
/// @param[out] toi The time of impact between the two points.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[in] tmax The maximum time to check for collisions.
/// @param[in] tolerance The error tolerance for the time of impact.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
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

/// @brief Computes the time of impact between two points using continuous collision detection.
/// @param[in] p0_t0 The initial position of the first point.
/// @param[in] p1_t0 The initial position of the second point.
/// @param[in] p0_t1 The final position of the first point.
/// @param[in] p1_t1 The final position of the second point.
/// @param[out] toi The time of impact between the two points.
/// @param[in] min_distance The minimum distance between the points.
/// @param[in] tmax The maximum time to check for collisions.
/// @param[in] tolerance The error tolerance for the time of impact.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
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

/// @brief Computes the time of impact between a point and an edge in 3D using continuous collision detection.
/// @param[in] p_t0 The initial position of the point.
/// @param[in] e0_t0 The initial position of the first endpoint of the edge.
/// @param[in] e1_t0 The initial position of the second endpoint of the edge.
/// @param[in] p_t1 The final position of the point.
/// @param[in] e0_t1 The final position of the first endpoint of the edge.
/// @param[in] e1_t1 The final position of the second endpoint of the edge.
/// @param[out] toi The time of impact between the point and the edge.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[in] tmax The maximum time to check for collisions.
/// @param[in] tolerance The error tolerance for the time of impact.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
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

/// @brief Computes the time of impact between a point and an edge using continuous collision detection.
/// @param[in] p_t0 The initial position of the point.
/// @param[in] e0_t0 The initial position of the first endpoint of the edge.
/// @param[in] e1_t0 The initial position of the second endpoint of the edge.
/// @param[in] p_t1 The final position of the point.
/// @param[in] e0_t1 The final position of the first endpoint of the edge.
/// @param[in] e1_t1 The final position of the second endpoint of the edge.
/// @param[out] toi The time of impact between the point and the edge.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[in] tmax The maximum time to check for collisions.
/// @param[in] tolerance The error tolerance for the time of impact.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
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

/// @brief Computes the time of impact between two edges in 3D using continuous collision detection.
/// @param[in] ea0_t0 The initial position of the first endpoint of the first edge.
/// @param[in] ea1_t0 The initial position of the second endpoint of the first edge.
/// @param[in] eb0_t0 The initial position of the first endpoint of the second edge.
/// @param[in] eb1_t0 The initial position of the second endpoint of the second edge.
/// @param[in] ea0_t1 The final position of the first endpoint of the first edge.
/// @param[in] ea1_t1 The final position of the second endpoint of the first edge.
/// @param[in] eb0_t1 The final position of the first endpoint of the second edge.
/// @param[in] eb1_t1 The final position of the second endpoint of the second edge.
/// @param[out] toi The time of impact between the two edges.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[in] tmax The maximum time to check for collisions.
/// @param[in] tolerance The error tolerance for the time of impact.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
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

/// @brief Computes the time of impact between a point and a triangle in 3D using continuous collision detection.
/// @param[in] p_t0 The initial position of the point.
/// @param[in] t0_t0 The initial position of the first vertex of the triangle.
/// @param[in] t1_t0 The initial position of the second vertex of the triangle.
/// @param[in] t2_t0 The initial position of the third vertex of the triangle.
/// @param[in] p_t1 The final position of the point.
/// @param[in] t0_t1 The final position of the first vertex of the triangle.
/// @param[in] t1_t1 The final position of the second vertex of the triangle.
/// @param[in] t2_t1 The final position of the third vertex of the triangle.
/// @param[out] toi The time of impact between the point and the triangle.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[in] tmax The maximum time to check for collisions.
/// @param[in] tolerance The error tolerance for the time of impact.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
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

/// @brief Perform the CCD strategy outlined by Li et al. [2020].
/// @param[in] ccd The continuous collision detection function.
/// @param[in] max_iterations The maximum number of iterations to perform.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[in] initial_distance The initial distance between the objects.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
/// @param[out] toi Output time of impact.
/// @return True if a collision was detected, false otherwise.
bool ccd_strategy(
    const std::function<bool(
        long /*max_iterations*/,
        double /*min_distance*/,
        bool /*no_zero_toi*/,
        double& /*toi*/)>& ccd,
    const long max_iterations,
    const double min_distance,
    const double initial_distance,
    const double conservative_rescaling,
    double& toi);

/// @brief Helper function to check if the initial distance is less than the minimum distance.
/// @param[in] initial_distance The initial distance between the objects.
/// @param[in] min_distance The minimum distance between the objects.
/// @param[out] toi Set to 0 if the initial distance is less than the minimum distance.
/// @return True if the initial distance is less than the minimum distance, false otherwise.
bool check_initial_distance(
    const double initial_distance, const double min_distance, double& toi);

} // namespace ipc
