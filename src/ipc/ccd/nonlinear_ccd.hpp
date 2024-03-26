#pragma once

#include <ipc/config.hpp>

#include <ipc/ccd/ccd.hpp>
#ifdef IPC_TOOLKIT_WITH_FILIB
#include <ipc/utils/interval.hpp>
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
    virtual VectorMax3I operator()(const filib::Interval& t) const = 0;

    /// @brief Compute the maximum distance from the nonlinear trajectory to a linearized trajectory
    /// @note This uses interval arithmetic to compute the maximum distance. If you know a tighter bound on the maximum distance, it is recommended to override this function.
    /// @param[in] t0 Start time of the trajectory
    /// @param[in] t1 End time of the trajectory
    virtual double
    max_distance_from_linear(const double t0, const double t1) const;
};
#endif

/// @brief Perform nonlinear CCD between two points moving along nonlinear trajectories.
/// @param[in] p0 First point's trajectory
/// @param[in] p1 Second point's trajectory
/// @param[out] toi Output time of impact
/// @param[in] tmax Maximum time to check for collision
/// @param[in] min_distance Minimum separation distance between the two points
/// @param[in] tolerance Tolerance for the linear CCD algorithm
/// @param[in] max_iterations Maximum number of iterations for the linear CCD algorithm
/// @param[in] conservative_rescaling  Conservative rescaling of the time of impact
/// @return True if the two points collide, false otherwise.
bool point_point_nonlinear_ccd(
    const NonlinearTrajectory& p0,
    const NonlinearTrajectory& p1,
    double& toi,
    const double tmax = 1.0,
    const double min_distance = 0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

/// @brief Perform nonlinear CCD between a point and a linear edge moving along nonlinear trajectories.
/// @param[in] p Point's trajectory
/// @param[in] e0 Edge's first endpoint's trajectory
/// @param[in] e1 Edge's second endpoint's trajectory
/// @param[out] toi Output time of impact
/// @param[in] tmax Maximum time to check for collision
/// @param[in] min_distance Minimum separation distance between the point and the edge
/// @param[in] tolerance Tolerance for the linear CCD algorithm
/// @param[in] max_iterations Maximum number of iterations for the linear CCD algorithm
/// @param[in] conservative_rescaling Conservative rescaling of the time of impact
/// @return True if the point and edge collide, false otherwise.
bool point_edge_nonlinear_ccd(
    const NonlinearTrajectory& p,
    const NonlinearTrajectory& e0,
    const NonlinearTrajectory& e1,
    double& toi,
    const double tmax = 1.0,
    const double min_distance = 0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

/// @brief Perform nonlinear CCD between two linear edges moving along nonlinear trajectories.
/// @param[in] ea0 First edge's first endpoint's trajectory
/// @param[in] ea1 First edge's second endpoint's trajectory
/// @param[in] eb0 Second edge's first endpoint's trajectory
/// @param[in] eb1 Second edge's second endpoint's trajectory
/// @param[out] toi Output time of impact
/// @param[in] tmax Maximum time to check for collision
/// @param[in] min_distance Minimum separation distance between the two edges
/// @param[in] tolerance Tolerance for the linear CCD algorithm
/// @param[in] max_iterations Maximum number of iterations for the linear CCD algorithm
/// @param[in] conservative_rescaling Conservative rescaling of the time of impact
/// @return True if the two edges collide, false otherwise.
bool edge_edge_nonlinear_ccd(
    const NonlinearTrajectory& ea0,
    const NonlinearTrajectory& ea1,
    const NonlinearTrajectory& eb0,
    const NonlinearTrajectory& eb1,
    double& toi,
    const double tmax = 1.0,
    const double min_distance = 0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

/// @brief Perform nonlinear CCD between a point and a linear triangle moving along nonlinear trajectories.
/// @param[in] p Point's trajectory
/// @param[in] t0 Triangle's first vertex's trajectory
/// @param[in] t1 Triangle's second vertex's trajectory
/// @param[in] t2 Triangle's third vertex's trajectory
/// @param[out] toi Output time of impact
/// @param[in] tmax Maximum time to check for collision
/// @param[in] min_distance Minimum separation distance between the two edges
/// @param[in] tolerance Tolerance for the linear CCD algorithm
/// @param[in] max_iterations Maximum number of iterations for the linear CCD algorithm
/// @param[in] conservative_rescaling Conservative rescaling of the time of impact
/// @return True if the point and triangle collide, false otherwise.
bool point_triangle_nonlinear_ccd(
    const NonlinearTrajectory& p,
    const NonlinearTrajectory& t0,
    const NonlinearTrajectory& t1,
    const NonlinearTrajectory& t2,
    double& toi,
    const double tmax = 1.0,
    const double min_distance = 0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

/// @brief Perform conservative piecewise linear CCD of a nonlinear trajectories.
/// @param[in] distance Return the distance for a given time in [0, 1].
/// @param[in] max_distance_from_linear Return the maximum distance from the linearized trajectory for a given time interval.
/// @param[in] linear_ccd Perform linear CCD on a given time interval.
/// @param[out] toi Output time of impact.
/// @param[in] tmax Maximum time to check for collision.
/// @param[in] min_distance Minimum separation distance between the objects.
/// @param[in] conservative_rescaling Conservative rescaling of the time of impact.
/// @return
bool conservative_piecewise_linear_ccd(
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
    const double tmax = 1.0,
    const double min_distance = 0,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

} // namespace ipc
