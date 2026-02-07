#pragma once

#include <ipc/dynamics/rigid/pose.hpp>

#include <functional>
#include <memory>

namespace ipc::rigid {

class BodyForces;
class GroundContact;
class ImplicitEuler;
class InertialTerm;
class RigidBodies;

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
    ///
    /// @param t_end End time
    /// @param callback Callback function to be called at each step.
    ///
    /// The callback function takes a boolean argument indicating whether the
    /// step was successful (i.e., did not fail to converge). The callback is
    /// called after each step, including the final step that reaches t_end.
    ///
    /// If the simulation is already complete (i.e., t >= t_end), the simulation
    /// will not run and the callback will not be called.
    ///
    /// @return True if the simulation ran successfully, false if it was terminated (e.g., due to convergence failure)
    bool run(
        // const double dt,
        const double t_end,
        const std::function<void(bool)>& callback = [](bool) { });

    /// @brief Step the simulation
    /// @param dt Time step
    /// @return True if the step was successful, false if the simulation should be terminated (e.g., due to convergence failure)
    bool step(
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
    double energy(Eigen::ConstRef<Eigen::VectorXd> x);
    Eigen::VectorXd gradient(Eigen::ConstRef<Eigen::VectorXd> x);
    Eigen::MatrixXd hessian(Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Bodies in the simulation
    std::shared_ptr<RigidBodies> m_bodies;

    /// @brief Time integrator for the simulation
    std::shared_ptr<ImplicitEuler> m_time_integrator;

    /// @brief Inertial term for the rigid body dynamics
    std::shared_ptr<InertialTerm> m_inertial_term;

    /// @brief Body forces acting on the rigid bodies
    std::shared_ptr<BodyForces> m_body_forces;

    /// @brief Ground contact handler
    std::shared_ptr<GroundContact> m_ground_contact;

    /// @brief History of poses at each time step
    std::list<std::vector<Pose>> m_pose_history;

    /// @brief t Current simulation time
    double m_t = 0.0;

    bool has_time_remaining(const double dt, const double t_end) const
    {
        return m_t < t_end && std::abs(m_t - t_end) > dt * 1e-3;
    }
};

} // namespace ipc::rigid