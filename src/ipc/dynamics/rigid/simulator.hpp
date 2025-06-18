#pragma once

#include <ipc/dynamics/rigid/pose.hpp>

#include <functional>
#include <memory>

namespace ipc::rigid {

class RigidBodies;
class ImplicitEuler;
class InertialTerm;

class Simulator {
public:
    /// @brief Create a simulation box
    /// @param corners Corners of the simulation box
    /// @param gravity Gravity vector
    Simulator(
        const std::shared_ptr<RigidBodies> _bodies,
        const std::vector<Pose>& initial_poses,
        const double dt);

    /// @brief Run the simulation
    /// @param t_end End time
    /// @param dt Time step
    /// @param callback
    void run(
        // const double dt,
        const double t_end,
        const std::function<void(void)> callback = []() { });

    /// @brief Step the simulation
    /// @param dt Time step
    void step(
        // double dt
    );

    /// @brief Reset the simulation to time t=0
    void reset();

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    const std::vector<std::vector<Pose>>& poses_history() const
    {
        return m_pose_history;
    }

    const std::shared_ptr<RigidBodies>& bodies() const { return m_bodies; }

    double t() const { return m_t; }

protected:
    /// @brief Bodies in the simulation
    std::shared_ptr<RigidBodies> m_bodies;

    /// @brief Time integrator for the simulation
    std::shared_ptr<ImplicitEuler> m_time_integrator;

    /// @brief Inertial term for the rigid body dynamics
    std::shared_ptr<InertialTerm> m_inertial_term;

    std::vector<std::vector<Pose>> m_pose_history;

    /// @brief t Current simulation time
    double m_t = 0.0;
};

} // namespace ipc::rigid