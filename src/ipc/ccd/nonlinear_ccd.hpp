#pragma once

#include <ipc/config.hpp>
#include <ipc/ccd/narrow_phase_ccd.hpp>
#include <ipc/ccd/tight_inclusion_ccd.hpp>
#ifdef IPC_TOOLKIT_WITH_FILIB
#include <ipc/math/interval.hpp>
#endif

#include <functional>
#include <memory>

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
    /// The default conservative rescaling value used to avoid taking steps
    /// exactly to impact.
    static constexpr double DEFAULT_CONSERVATIVE_RESCALING = 0.8;

    /// @brief Construct a new NonlinearCCD object.
    /// @param conservative_rescaling Conservative rescaling of the time of impact.
    /// @param narrow_phase_ccd The narrow-phase CCD used for the per-piece linear CCD (if null, a default-constructed TightInclusionCCD is used).
    /// @note The constructor sets NarrowPhaseCCD::component_mode on the given object — telling it that it is a component of this piecewise driver — so do not share the instance with standalone (non-component) users.
    NonlinearCCD(
        const double conservative_rescaling = DEFAULT_CONSERVATIVE_RESCALING,
        std::shared_ptr<NarrowPhaseCCD> narrow_phase_ccd = nullptr);

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
    virtual bool point_triangle_ccd(
        const NonlinearTrajectory& p,
        const NonlinearTrajectory& t0,
        const NonlinearTrajectory& t1,
        const NonlinearTrajectory& t2,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const;

    /// @brief Set the narrow-phase CCD used for the per-piece linear CCD.
    /// @note Sets NarrowPhaseCCD::component_mode on the given object —
    /// telling it that it is a component of this piecewise driver — so do not
    /// share the instance with standalone (non-component) users.
    /// @param narrow_phase_ccd The narrow-phase CCD to use (if null, a default-constructed TightInclusionCCD is used).
    void set_narrow_phase_ccd(std::shared_ptr<NarrowPhaseCCD> narrow_phase_ccd);

    /// @brief The narrow-phase CCD used for the per-piece linear CCD (in component mode; never null).
    std::shared_ptr<const NarrowPhaseCCD> narrow_phase_ccd() const
    {
        return m_narrow_phase_ccd;
    }

    /// @brief Conservative rescaling of the time of impact.
    double conservative_rescaling;

private:
    /// @brief The narrow-phase CCD used for the per-piece linear CCD (in component mode; never null).
    /// @note The piecewise driver passes each piece's linearization error as the minimum separation distance, so the narrow phase must support minimum separation (both TightInclusionCCD and AdditiveCCD do).
    std::shared_ptr<const NarrowPhaseCCD> m_narrow_phase_ccd;

public:
    /// @brief Perform conservative piecewise linear CCD of a nonlinear trajectories.
    /// @param[in] distance Return the distance for a given time in [0, 1].
    /// @param[in] max_distance_from_linear Return the maximum distance from the linearized trajectory for a given time interval.
    /// @param[in] linear_ccd Perform linear CCD on a given time interval. Receives the piece interval, the deviation-inflated minimum separation, and a cm_no_zero_toi hint (true on the first piece, whose start is known to be more than min_distance apart).
    /// @param[out] toi Output time of impact.
    /// @param[in] min_distance Minimum separation distance between the objects.
    /// @param[in] tmax Maximum time to check for collision.
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
            const bool /*cm_no_zero_toi*/,
            double& /*toi*/)>& linear_ccd,
        double& toi,
        const double min_distance = 0,
        const double tmax = 1.0,
        const double conservative_rescaling = DEFAULT_CONSERVATIVE_RESCALING);
};

} // namespace ipc
