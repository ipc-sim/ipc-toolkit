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
    /// @param bodies Rigid bodies in the simulation
    /// @param initial_poses Initial poses of the rigid bodies
    /// @param dt Time step
    Simulator(
        const std::shared_ptr<RigidBodies>& bodies,
        const std::vector<Pose>& initial_poses,
        const double dt);

    /// @brief Run the simulation
    /// @param t_end End time
    /// @param callback Callback function to be called at each step
    void run(
        // const double dt,
        const double t_end,
        const std::function<void(void)>& callback = []() { });

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

    const std::list<std::vector<Pose>>& pose_history() const
    {
        return m_pose_history;
    }

    std::shared_ptr<const RigidBodies> bodies() const { return m_bodies; }

    double t() const { return m_t; }
    // double dt() const { return m_time_integrator->dt; }
    // void set_dt(const double dt) { m_time_integrator->set_dt(dt); }

    // VectorMax3d gravity() const { return m_inertial_term->gravity(); }
    // void set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
    // {
    //     m_inertial_term->set_gravity(gravity);
    // }

protected:
    /// @brief Bodies in the simulation
    std::shared_ptr<RigidBodies> m_bodies;

    /// @brief Time integrator for the simulation
    std::shared_ptr<ImplicitEuler> m_time_integrator;

    /// @brief Inertial term for the rigid body dynamics
    std::shared_ptr<InertialTerm> m_inertial_term;

    std::list<std::vector<Pose>> m_pose_history;

    /// @brief t Current simulation time
    double m_t = 0.0;
};

} // namespace ipc::rigid