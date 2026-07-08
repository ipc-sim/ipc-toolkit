#pragma once

#include <ipc/dynamics/rigid/pose.hpp>

#include <deque>
#include <limits>

namespace ipc::demo {

/// @brief Drives a KINEMATIC rigid body to a per-step target pose.
///
/// A kinematic body follows either a scripted sequence of absolute poses (one
/// per step) or, when no script is given, its own prescribed velocity. After
/// an optional maximum drive time the simulator converts the body to STATIC.
///
/// The driver holds no simulation state beyond its script and remaining time;
/// the current pose and velocity are supplied by the caller (extracted from
/// the time integrator), so the driver is agnostic to the DOF layout.
class KinematicDriver {
public:
    KinematicDriver() = default;

    /// @brief A velocity-driven driver (target integrated from the body's
    /// prescribed velocity each step).
    /// @param max_time Drive time before the body converts to STATIC.
    static KinematicDriver velocity_driven(
        const double max_time = std::numeric_limits<double>::infinity())
    {
        return KinematicDriver({}, max_time);
    }

    /// @brief A scripted driver (target is the front pose, advanced per step).
    /// @param poses Absolute target poses, one per step.
    /// @param max_time Drive time before the body converts to STATIC.
    static KinematicDriver scripted(
        std::deque<rigid::Pose> poses,
        const double max_time = std::numeric_limits<double>::infinity())
    {
        return KinematicDriver(std::move(poses), max_time);
    }

    /// @brief The target pose for this step.
    /// @param current The body's current pose (position + rotation
    ///     vector/angle).
    /// @param velocity The body's prescribed velocity (linear in .position,
    ///     angular ω in .rotation); used only when there is no script.
    /// @param dt Time step.
    rigid::Pose target(
        const rigid::Pose& current,
        const rigid::Pose& velocity,
        const double dt) const;

    /// @brief Advance the driver by one step (pop the front scripted pose and
    /// count down the drive time).
    void step(const double dt)
    {
        m_max_time -= dt;
        if (!m_scripted_poses.empty()) {
            m_scripted_poses.pop_front();
        }
    }

    /// @brief Whether the drive time has elapsed (the body should convert to
    /// STATIC before this step).
    bool expired(const double dt) const { return m_max_time < 0.5 * dt; }

    /// @brief Whether this driver follows a scripted pose sequence.
    bool is_scripted() const { return !m_scripted_poses.empty(); }

    /// @brief Remaining drive time.
    double max_time() const { return m_max_time; }

    /// @brief The scripted target poses (empty ⇒ velocity-driven).
    const std::deque<rigid::Pose>& scripted_poses() const
    {
        return m_scripted_poses;
    }

private:
    KinematicDriver(std::deque<rigid::Pose> poses, const double max_time)
        : m_scripted_poses(std::move(poses))
        , m_max_time(max_time)
    {
    }

    /// @brief Scripted absolute target poses (one per step); empty ⇒
    /// velocity-driven.
    std::deque<rigid::Pose> m_scripted_poses;

    /// @brief Remaining drive time before the body converts to STATIC.
    double m_max_time = std::numeric_limits<double>::infinity();
};

} // namespace ipc::demo
