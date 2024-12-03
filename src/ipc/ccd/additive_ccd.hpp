//
// NOTE: These methods are provided for reference comparison with [Li et al.
// 2021] and is not utilized by the high-level functionality. In compairson to
// Tight Inclusion CCD, this CCD method is not provably conservative and so can
// potentially produce false negatives (i.e., miss collisions) due to
// floating-point rounding error. However, it is much faster than Tight
// Inclusion CCD (>100×) and very robust due to the gaps and conservative
// rescaling used.
//

#pragma once

#include <ipc/ccd/narrow_phase_ccd.hpp>

namespace ipc {

/// @brief Additive Continuous Collision Detection (CCD) from [Li et al. 2021].
class AdditiveCCD : public NarrowPhaseCCD {
public:
    /// The default maximum number of iterations used with Tight-Inclusion CCD.
    static constexpr long DEFAULT_MAX_ITERATIONS = 10'000'000l;
    /// Unlimitted number of iterations.
    static constexpr long UNLIMITTED_ITERATIONS = -1;
    /// The default conservative rescaling value used to avoid taking steps
    /// exactly to impact. Value choosen to based on [Li et al. 2021].
    static constexpr double DEFAULT_CONSERVATIVE_RESCALING = 0.9;

    /// @brief Construct a new AdditiveCCD object.
    /// @param conservative_rescaling The conservative rescaling of the time of impact.
    AdditiveCCD(
        const long max_iterations = UNLIMITTED_ITERATIONS,
        const double conservative_rescaling = DEFAULT_CONSERVATIVE_RESCALING);

    /// @brief Computes the time of impact between two points using continuous collision detection.
    /// @param p0_t0 The initial position of the first point.
    /// @param p1_t0 The initial position of the second point.
    /// @param p0_t1 The final position of the first point.
    /// @param p1_t1 The final position of the second point.
    /// @param[out] toi The time of impact between the two points.
    /// @param min_distance The minimum distance between two objects.
    /// @param tmax The maximum time to check for collisions.
    /// @param conservative_rescaling The conservative rescaling of the time of impact.
    /// @return True if a collision was detected, false otherwise.
    bool point_point_ccd(
        const VectorMax3d& p0_t0,
        const VectorMax3d& p1_t0,
        const VectorMax3d& p0_t1,
        const VectorMax3d& p1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const override;

    /// @brief Computes the time of impact between a point and an edge using continuous collision detection.
    /// @param p_t0 The initial position of the point.
    /// @param e0_t0 The initial position of the first endpoint of the edge.
    /// @param e1_t0 The initial position of the second endpoint of the edge.
    /// @param p_t1 The final position of the point.
    /// @param e0_t1 The final position of the first endpoint of the edge.
    /// @param e1_t1 The final position of the second endpoint of the edge.
    /// @param[out] toi The time of impact between the point and the edge.
    /// @param min_distance The minimum distance between two objects.
    /// @param tmax The maximum time to check for collisions.
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
        const double tmax = 1.0) const override;

    /// @brief Computes the time of impact between a point and a triangle using continuous collision detection.
    /// @param p_t0 The initial position of the point.
    /// @param t0_t0 The initial position of the first vertex of the triangle.
    /// @param t1_t0 The initial position of the second vertex of the triangle.
    /// @param t2_t0 The initial position of the third vertex of the triangle.
    /// @param p_t1 The final position of the point.
    /// @param t0_t1 The final position of the first vertex of the triangle.
    /// @param t1_t1 The final position of the second vertex of the triangle.
    /// @param t2_t1 The final position of the third vertex of the triangle.
    /// @param[out] toi The time of impact between the point and the triangle.
    /// @param min_distance The minimum distance between two objects.
    /// @param tmax The maximum time to check for collisions.
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
        const double tmax = 1.0) const override;

    /// @brief Computes the time of impact between two edges using continuous collision detection.
    /// @param ea0_t0 The initial position of the first endpoint of the first edge.
    /// @param ea1_t0 The initial position of the second endpoint of the first edge.
    /// @param eb0_t0 The initial position of the first endpoint of the second edge.
    /// @param eb1_t0 The initial position of the second endpoint of the second edge.
    /// @param ea0_t1 The final position of the first endpoint of the first edge.
    /// @param ea1_t1 The final position of the second endpoint of the first edge.
    /// @param eb0_t1 The final position of the first endpoint of the second edge.
    /// @param eb1_t1 The final position of the second endpoint of the second edge.
    /// @param[out] toi The time of impact between the two edges.
    /// @param min_distance The minimum distance between two objects.
    /// @param tmax The maximum time to check for collisions.
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
        const double tmax = 1.0) const override;

    /// @brief Maximum number of iterations.
    long max_iterations;

    /// @brief The conservative rescaling value used to avoid taking steps exactly to impact.
    double conservative_rescaling;

private:
    /// @brief Computes the time of impact between two objects using additive continuous collision detection.
    /// @param distance_squared A function that computes the squared distance between the two objects at a given time.
    /// @param[out] toi The time of impact between the two objects.
    /// @param min_distance The minimum distance between the objects.
    /// @param tmax The maximum time to check for collisions.
    /// @param conservative_rescaling The amount to rescale the objects by to ensure conservative advancement.
    /// @return True if a collision was detected, false otherwise.
    bool additive_ccd(
        VectorMax12d x,
        const VectorMax12d& dx,
        const std::function<double(const VectorMax12d&)>& distance_squared,
        const double max_disp_mag,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const;
};

} // namespace ipc