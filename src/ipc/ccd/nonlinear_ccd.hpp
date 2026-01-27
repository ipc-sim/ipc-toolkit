#pragma once

#include <ipc/config.hpp>
#include <ipc/ccd/tight_inclusion_ccd.hpp>
#ifdef IPC_TOOLKIT_WITH_FILIB
#include <ipc/math/interval.hpp>
#endif

#include <functional>

namespace ipc {

/// @brief A nonlinear trajectory is a function that maps time to a point in space.
class NonlinearTrajectory {
public:
    virtual ~NonlinearTrajectory() = default;

    /// @brief Compute the point's position at time t
    virtual VectorMax3d operator()(const double t) const = 0;

    /// @brief Compute the maximum distance from the nonlinear trajectory to a linearized trajectory
    /// @param[in] t0 Start time of the trajectory
    /// @param[in] t1 End time of the trajectory
    virtual double
    max_distance_from_linear(const double t0, const double t1) const = 0;
};

#ifdef IPC_TOOLKIT_WITH_FILIB
/// @brief A nonlinear trajectory with an implementation of the max_distance_from_linear function using interval arithmetic.
class IntervalNonlinearTrajectory : virtual public NonlinearTrajectory {
public:
    virtual ~IntervalNonlinearTrajectory() = default;

    using NonlinearTrajectory::operator();

    /// @brief Compute the point's position over a time interval t
    /// @param[in] t The time interval
    /// @return The point's position at time t as an interval
    virtual VectorMax3I operator()(const filib::Interval& t) const = 0;

    /// @brief Compute the maximum distance from the nonlinear trajectory to a linearized trajectory
    /// @note This uses interval arithmetic to compute the maximum distance. If you know a tighter bound on the maximum distance, it is recommended to override this function.
    /// @param[in] t0 Start time of the trajectory
    /// @param[in] t1 End time of the trajectory
    /// @return The maximum distance from the nonlinear trajectory to a linearized trajectory
    virtual double
    max_distance_from_linear(const double t0, const double t1) const;
};
#endif

class NonlinearCCD {
public:
    /// The default tolerance used with Tight-Inclusion CCD.
    static constexpr double DEFAULT_TOLERANCE = 1e-6;
    /// The default maximum number of iterations used with Tight-Inclusion CCD.
    static constexpr long DEFAULT_MAX_ITERATIONS = 10'000'000l;
    /// The default conservative rescaling value used to avoid taking steps
    /// exactly to impact.
    static constexpr double DEFAULT_CONSERVATIVE_RESCALING = 0.8;

    NonlinearCCD(
        const double tolerance = DEFAULT_TOLERANCE,
        const long max_iterations = DEFAULT_MAX_ITERATIONS,
        const double conservative_rescaling = DEFAULT_CONSERVATIVE_RESCALING);

    virtual ~NonlinearCCD() = default;

    /// @brief Perform nonlinear CCD between two points moving along nonlinear trajectories.
    /// @param[in] p0 First point's trajectory
    /// @param[in] p1 Second point's trajectory
    /// @param[out] toi Output time of impact
    /// @param[in] min_distance Minimum separation distance between the two points
    /// @param[in] tmax Maximum time to check for collision
    /// @return True if the two points collide, false otherwise.
    virtual bool point_point_ccd(
        const NonlinearTrajectory& p0,
        const NonlinearTrajectory& p1,
        double& toi,
        const double min_distance = 0,
        const double tmax = 1.0) const;

    /// @brief Perform nonlinear CCD between a point and a linear edge moving along nonlinear trajectories.
    /// @param[in] p Point's trajectory
    /// @param[in] e0 Edge's first endpoint's trajectory
    /// @param[in] e1 Edge's second endpoint's trajectory
    /// @param[out] toi Output time of impact
    /// @param[in] min_distance Minimum separation distance between the point and the edge
    /// @param[in] tmax Maximum time to check for collision
    /// @return True if the point and edge collide, false otherwise.
    virtual bool point_edge_ccd(
        const NonlinearTrajectory& p,
        const NonlinearTrajectory& e0,
        const NonlinearTrajectory& e1,
        double& toi,
        const double min_distance = 0,
        const double tmax = 1.0) const;

    /// @brief Perform nonlinear CCD between two linear edges moving along nonlinear trajectories.
    /// @ingroup ccd
    /// @param[in] ea0 First edge's first endpoint's trajectory
    /// @param[in] ea1 First edge's second endpoint's trajectory
    /// @param[in] eb0 Second edge's first endpoint's trajectory
    /// @param[in] eb1 Second edge's second endpoint's trajectory
    /// @param[out] toi Output time of impact
    /// @param[in] min_distance Minimum separation distance between the two edges
    /// @param[in] tmax Maximum time to check for collision
    /// @return True if the two edges collide, false otherwise.
    virtual bool edge_edge_ccd(
        const NonlinearTrajectory& ea0,
        const NonlinearTrajectory& ea1,
        const NonlinearTrajectory& eb0,
        const NonlinearTrajectory& eb1,
        double& toi,
        const double min_distance = 0,
        const double tmax = 1.0) const;

    /// @brief Perform nonlinear CCD between a point and a linear triangle moving along nonlinear trajectories.
    /// @param[in] p Point's trajectory
    /// @param[in] t0 Triangle's first vertex's trajectory
    /// @param[in] t1 Triangle's second vertex's trajectory
    /// @param[in] t2 Triangle's third vertex's trajectory
    /// @param[out] toi Output time of impact
    /// @param[in] min_distance Minimum separation distance between the two edges
    /// @param[in] tmax Maximum time to check for collision
    /// @return True if the point and triangle collide, false otherwise.
    bool point_triangle_ccd(
        const NonlinearTrajectory& p,
        const NonlinearTrajectory& t0,
        const NonlinearTrajectory& t1,
        const NonlinearTrajectory& t2,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const;

    /// @brief Solver tolerance.
    double tolerance;

    /// @brief Maximum number of iterations.
    long max_iterations;

    /// @brief Conservative rescaling of the time of impact.
    double conservative_rescaling;

public:
    /// @brief Perform conservative piecewise linear CCD of a nonlinear trajectories.
    /// @param[in] distance Return the distance for a given time in [0, 1].
    /// @param[in] max_distance_from_linear Return the maximum distance from the linearized trajectory for a given time interval.
    /// @param[in] linear_ccd Perform linear CCD on a given time interval.
    /// @param[out] toi Output time of impact.
    /// @param[in] tmax Maximum time to check for collision.
    /// @param[in] min_distance Minimum separation distance between the objects.
    /// @param[in] conservative_rescaling Conservative rescaling of the time of impact.
    /// @return True if a collision was detected, false otherwise.
    static bool conservative_piecewise_linear_ccd(
        const std::function<double(const double)>& distance,
        const std::function<double(const double, const double)>&
            max_distance_from_linear,
        const std::function<bool(
            const double /*ti0*/,
            const double /*ti1*/,
            const double /*min_distance*/,
            const bool /*no_zero_toi*/,
            double& /*toi*/)>& linear_ccd,
        double& toi,
        const double min_distance = 0,
        const double tmax = 1.0,
        const double conservative_rescaling = DEFAULT_CONSERVATIVE_RESCALING);
};

} // namespace ipc
