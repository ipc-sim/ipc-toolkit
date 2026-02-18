#pragma once

#include <ipc/dynamics/rigid/pose.hpp>

#include <functional>
#include <memory>

namespace ipc {
class AdditiveCCD;
class BarrierPotential;
class BroadPhase;
class Candidates;
class NormalCollisions;
} // namespace ipc

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

    /// @brief Destructor
    ~Simulator();

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
    // Objective, gradient, and hessian for the optimization problem
    // -----------------------------------------------------------------------

    /// @brief Compute the energy of the system at the given state x.
    /// @param x State vector containing the positions and orientations of all rigid bodies.
    /// @return The total incremental potential of the system, including inertial terms, body forces, ground contact, and barrier potential for collisions.
    double energy(Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Compute the gradient of the energy with respect to the state vector x.
    /// @param x State vector containing the positions and orientations of all rigid bodies.
    /// @return The gradient of the total incremental potential with respect to x.
    Eigen::VectorXd gradient(Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Compute the Hessian of the energy with respect to the state vector x.
    /// @param x State vector containing the positions and orientations of all rigid bodies.
    /// @param project_to_psd If true, project the Hessian to be positive semi-definite to improve convergence of the optimization. This is often necessary when using second-order optimization methods, as a non-PSD Hessian can lead to divergence or convergence to saddle points.
    /// @return The Hessian of the total incremental potential with respect to x.
    Eigen::MatrixXd
    hessian(Eigen::ConstRef<Eigen::VectorXd> x, bool project_to_psd = true);

    // -----------------------------------------------------------------------
    // Collision handling
    // -----------------------------------------------------------------------

    /// @brief Update the simulation state before stepping, including updating the collision sets based on the current state.
    void initialize_step();

    /// @brief Update the collision sets based on the current state x.
    /// @param x State vector containing the positions and orientations of all rigid bodies.
    /// @param update_candidates If true, also update the candidate collisions based on the current state. If false, only update the normal collision set based on the current candidate collisions.
    void update_collisions(
        Eigen::ConstRef<Eigen::VectorXd> x, bool update_candidates = true);

    /// @brief Update the candidate collisions based on the current state at time t0 and the predicted state at time t1.
    /// @param x_t0 State vector at time t0 containing the positions and orientations of all rigid bodies.
    /// @param x_t1 State vector at time t1 containing the positions and orientations of all rigid bodies.
    void update_candidates(
        Eigen::ConstRef<Eigen::VectorXd> x_t0,
        Eigen::ConstRef<Eigen::VectorXd> x_t1);

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    const std::list<std::vector<Pose>>& pose_history() const
    {
        return m_pose_history;
    }

    const RigidBodies& bodies() const { return *m_bodies; }

    double t() const { return m_t; }
    // double dt() const { return m_time_integrator->dt; }
    // void set_dt(const double dt) { m_time_integrator->set_dt(dt); }

    // VectorMax3d gravity() const { return m_inertial_term->gravity(); }
    // void set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
    // {
    //     m_inertial_term->set_gravity(gravity);
    // }

    /// @brief Get the barrier potential for collision handling
    const BarrierPotential& barrier_potential() const
    {
        return *m_barrier_potential;
    }

    /// @brief Get the candidate collisions for the current time step
    const Candidates& candidates() const { return *m_candidates; }

    /// @brief Get the normal collision set for the current time step.
    const NormalCollisions& normal_collisions() const
    {
        return *m_normal_collisions;
    }

protected:
    /// @brief Check if there is time remaining in the simulation to take another step.
    ///
    /// This function accounts for floating point precision issues by checking
    /// if m_t is less than t_end and that the difference between m_t and t_end
    /// is greater than a small fraction of dt.
    ///
    /// @param dt Time step size
    /// @param t_end End time of the simulation
    /// @return True if there is time remaining to take another step, false otherwise
    bool has_time_remaining(const double dt, const double t_end) const
    {
        return m_t < t_end && std::abs(m_t - t_end) > dt * 1e-3;
    }

    /// @brief Bodies in the simulation
    std::shared_ptr<RigidBodies> m_bodies;

    /// @brief Time integrator for the simulation
    std::shared_ptr<ImplicitEuler> m_time_integrator;

    /// @brief Inertial term for the rigid body dynamics
    std::unique_ptr<InertialTerm> m_inertial_term;

    /// @brief Body forces acting on the rigid bodies
    std::unique_ptr<BodyForces> m_body_forces;

    /// @brief History of poses at each time step
    std::list<std::vector<Pose>> m_pose_history;

    /// @brief Barrier potential for collision handling
    std::unique_ptr<BarrierPotential> m_barrier_potential;

    /// @brief Candidate collisions for the current time step
    std::unique_ptr<Candidates> m_candidates;

    /// @brief Normal collision set
    std::unique_ptr<NormalCollisions> m_normal_collisions;

    /// @brief t Current simulation time
    double m_t = 0.0;

    /// @brief Broad phase collision handler for efficiently finding potential collisions.
    std::unique_ptr<BroadPhase> m_broad_phase;

    /// @brief Additive CCD handler for continuous collision detection.
    std::unique_ptr<AdditiveCCD> m_additive_ccd;
};

} // namespace ipc::rigid