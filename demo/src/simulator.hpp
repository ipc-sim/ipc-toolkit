#pragma once

#include "kinematic_driver.hpp"
#include "kinematics.hpp"

#include <ipc/dynamics/affine/body_forces.hpp>
#include <ipc/dynamics/affine/inertial_term.hpp>
#include <ipc/dynamics/affine/orthogonality_potential.hpp>
#include <ipc/dynamics/rigid/body_forces.hpp>
#include <ipc/dynamics/rigid/inertial_term.hpp>
#include <ipc/dynamics/rigid/pose.hpp>

#include <nlohmann/json.hpp>
#include <polysolve/nonlinear/Problem.hpp>
#include <spdlog/common.h>

#include <cstdint>
#include <functional>
#include <list>
#include <memory>
#include <optional>

namespace spdlog {
class logger;
} // namespace spdlog

namespace ipc {
class BarrierPotential;
class BroadPhase;
class Candidates;
class FrictionPotential;
class NormalCollisions;
class TangentialCollisions;
} // namespace ipc

namespace ipc::rigid {
class RigidBodies;
} // namespace ipc::rigid

namespace ipc::affine {
class JointConstraints;
} // namespace ipc::affine

#include <ipc/dynamics/affine/augmented_lagrangian.hpp>
#include <ipc/dynamics/rigid/augmented_lagrangian.hpp>
#include <ipc/dynamics/rigid/rigid_body.hpp>

namespace ipc::dynamics {
class ImplicitTimeIntegrator;
} // namespace ipc::dynamics

namespace polysolve::nonlinear {
class Solver;
} // namespace polysolve::nonlinear

namespace ipc::demo {

/// @brief Body dynamics simulator with a rigid/affine toggle.
///
/// Minimizes an incremental potential (inertia + body forces
/// [+ orthogonality in affine mode] + barrier) over the DOFs of the selected
/// dynamics model with polysolve's Newton solver:
/// - RIGID [Ferguson et al. 2021]: 6-DOF [p; θ] per body with
///   curved-trajectory (nonlinear) CCD.
/// - AFFINE [Lan et al. 2022]: 12-DOF [p; vec(A)] per body with a stiff
///   orthogonality potential and linear CCD (exact for affine trajectories);
///   supports linear equality joint constraints [Chen et al. 2022].
///
/// The simulator is itself the polysolve::nonlinear::Problem it solves: the
/// Problem interface (value/gradient/hessian, CCD max_step_size, collision
/// state hooks) is implemented below. When joint constraints are present the
/// Problem interface operates in the reduced coordinates z (x = Uz with the
/// first m entries pinned [Chen et al. 2022, §4.1]); otherwise its
/// coordinates are the DOFs x directly.
class Simulator : public polysolve::nonlinear::Problem {
public:
    /// @brief The dynamics model used to simulate the bodies.
    enum class BodyDynamics : std::uint8_t {
        /// @brief Rigid body dynamics [Ferguson et al. 2021]: 6-DOF per body
        /// ([position; rotation vector]) with curved-trajectory (nonlinear)
        /// CCD.
        RIGID,
        /// @brief Affine body dynamics [Lan et al. 2022]: 12-DOF per body
        /// ([position; vec(A) column-major]) with a stiff orthogonality
        /// potential and linear CCD (exact for affine per-iterate
        /// trajectories).
        AFFINE,
    };

    /// @brief User-configurable simulation parameters.
    struct Settings {
        /// @brief The dynamics model used for all bodies.
        BodyDynamics body_dynamics = BodyDynamics::RIGID;
        /// @brief Order of the BDF time integrator (1 ≤ n ≤ 6; 1 is implicit
        /// Euler). The effective order ramps up from 1 over the first n steps
        /// as the required history accumulates.
        int bdf_order = 1;
        /// @brief Barrier activation distance.
        double dhat = 1e-3;
        /// @brief Gravitational acceleration.
        VectorMax3d gravity = Eigen::Vector3d(0, -9.81, 0);
        /// @brief Scale used to premultiply the minimum barrier stiffness.
        double min_barrier_stiffness_scale = 1e11;
        /// @brief Adaptive barrier stiffness updates if the minimum distance
        /// is less than this fraction of the bounding box diagonal.
        double dhat_epsilon_scale = 1e-9;
        /// @brief Coefficient of friction (0 disables friction)
        /// [Ferguson et al. 2021].
        /// @note Friction forces are proportional to the (lagged) barrier
        /// normal forces, so their accuracy follows the contact-force
        /// accuracy of the solve: for quantitative friction, iterate the
        /// lagging to convergence (friction_iterations < 0) and tighten
        /// velocity_conv_tol (e.g., 1e-3) so the solver cannot stop while
        /// the contact forces are unbalanced.
        double friction_coefficient = 0;
        /// @brief Smooth friction mollifier speed bound εᵥ (m/s): sliding
        /// speeds below this are treated as static [Ferguson et al. 2021].
        double static_friction_speed_bound = 1e-3;
        /// @brief Friction lagging iterations: the friction operators (tangent
        /// bases and normal forces) are frozen during each solve and rebuilt
        /// between solves. 0 disables friction; > 0 caps the number of solves
        /// per step (1 = a single lagged solve); < 0 iterates until the
        /// momentum balance ‖∇E + κ∇B + ∇D‖ ≤ 1e-2·bbox_diagonal (up to
        /// max_outer_iterations).
        int friction_iterations = 1;
        /// @brief Hard cap on the outer (friction-lagging / augmented
        /// Lagrangian) solves per step.
        int max_outer_iterations = 50;
        /// @brief Initial penalty of the kinematic-body augmented Lagrangian
        /// (reset each step) [Ferguson et al. 2021].
        double al_initial_penalty = 1e3;
        /// @brief Maximum AL penalty (stop doubling).
        double al_max_penalty = 1e8;
        /// @brief AL progress η at which a channel is satisfied (DOFs freeze).
        double al_satisfied_progress = 0.999;
        /// @brief AL progress below which the penalty doubles (otherwise the
        /// multipliers update).
        double al_stall_progress = 0.99;
        /// @brief Stiffness of the orthogonality potential (affine mode).
        /// @note The potential also scales by each body's volume and enters
        /// the incremental potential as Δt²V⊥ [Lan et al. 2022, Eq. 9].
        /// Lan et al. recommend κ > 100 GPa to keep deformation negligible.
        double orthogonality_stiffness = 1e11;
        /// @brief Whether to use area weighting for the collision set.
        bool use_area_weighting = true;
        /// @brief If true, stop the simulation when a step fails to converge.
        bool abort_on_convergence_failure = true;
        /// @brief Velocity convergence tolerance: converged when the proposed
        /// Newton step moves every vertex slower than this × the bounding box
        /// diagonal (per unit time) [Ferguson et al. 2021]. Measured in world
        /// space through the step_norm() override and enforced through
        /// polysolve's Δx criterion (solver_params["x_delta_tol"], set
        /// automatically unless already present). Set ≤ 0 to disable.
        double velocity_conv_tol = 1e-2;
        /// @brief polysolve nonlinear solver parameters.
        /// Validated by polysolve; entries not set here get the defaults of
        /// https://github.com/polyfem/polysolve/blob/main/nonlinear-solver-spec.json
        // NOTE: The strategy ladder starts PSD-projected (like the original
        //       Rigid IPC solver) instead of polysolve's default plain Newton,
        //       whose indefinite-Hessian direction is discarded at every
        //       contact step; gradient descent is excluded because its
        //       barrier-scaled directions (‖Δx‖ ~ 1e14) make the CCD-filtered
        //       line search pathologically slow and never produce a usable
        //       step. "Backtracking" accepts any energy decrease, matching the
        //       original Rigid IPC line search; Armijo's sufficient-decrease
        //       condition is unattainable along CCD-capped steps into the stiff
        //       barrier (the predicted decrease c·α·|gᵀd| dwarfs the achievable
        //       one), which stalls the solver at first contact.
        // NOTE: residual_tolerance gates directions on the ABSOLUTE linear
        //       solve residual ‖HΔx + g‖, which is meaningless at our scales
        //       (κ ~ 1e11 puts ‖g‖ at ~1e13); it is effectively disabled,
        //       leaving the descent-direction check and the line search as the
        //       safety net (exactly the original Rigid IPC solver's guards).
        nlohmann::json solver_params = {
            { "solver",
              { { { "type", "ProjectedNewton" },
                  { "residual_tolerance", 1e100 } },
                { { "type", "RegularizedProjectedNewton" },
                  { "residual_tolerance", 1e100 } } } },
            { "allow_non_grad_convergence", true },
            { "max_iterations", 100 },
            { "grad_norm_tol", 1e-6 },
            { "line_search",
              { { "method", "Backtracking" }, { "min_step_size", 1e-12 } } },
        };
        /// @brief polysolve linear solver parameters.
        /// @note SimplicialLDLT (a direct solver, as in the original Rigid IPC
        ///       solver) is selected explicitly: polysolve's fallback default
        ///       is the iterative Eigen::BiCGSTAB, which does not converge on
        ///       the severely ill-conditioned barrier systems (entries spanning
        ///       ~1e2 to ~1e14), causing spurious direction failures.
        nlohmann::json linear_solver_params = {
            { "solver", "Eigen::SimplicialLDLT" },
        };
        /// @brief Verbosity of polysolve's own logger (separate from the
        /// toolkit's logger).
        /// @note Silenced by default: polysolve logs a per-solve "Finished"
        /// report at *error* level whenever the stop was not the gradient
        /// criterion — i.e., on every velocity-criterion stop — and
        /// per-iteration [timing] lines at debug. Lower this (e.g., to warn
        /// or debug) to see them; solve failures are always reported through
        /// the toolkit's logger regardless.
        spdlog::level::level_enum solver_log_level = spdlog::level::info;
    };

    /// @brief Create a simulator.
    /// @param bodies Bodies in the simulation (rest geometry and mass properties)
    /// @param initial_poses Initial (rigid) poses of the bodies
    /// @param dt Time step
    Simulator(
        const std::shared_ptr<rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses,
        const double dt);

    /// @copydoc Simulator(const std::shared_ptr<rigid::RigidBodies>&, const std::vector<rigid::Pose>&, const double)
    /// @param settings Simulation settings
    Simulator(
        const std::shared_ptr<rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses,
        const double dt,
        const Settings& settings);

    /// @brief Create a simulator with initial velocities.
    /// @param initial_velocities Initial velocities of the bodies: position
    ///     is the linear velocity and rotation is the world-frame angular
    ///     velocity ω (Q̇ = [ω]× Q).
    Simulator(
        const std::shared_ptr<rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses,
        const std::vector<rigid::Pose>& initial_velocities,
        const double dt,
        const Settings& settings);

    /// @brief Create a simulator with joint constraints.
    /// @param joints Linear equality joint constraints (finalized by this
    ///     constructor if not already).
    /// @throws std::invalid_argument unless settings.body_dynamics == AFFINE
    ///     (material-point constraints are nonlinear in the rigid rotation
    ///     vectors).
    Simulator(
        const std::shared_ptr<rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses,
        const std::shared_ptr<affine::JointConstraints>& joints,
        const double dt);

    /// @copydoc Simulator(const std::shared_ptr<rigid::RigidBodies>&, const std::vector<rigid::Pose>&, const std::shared_ptr<affine::JointConstraints>&, const double)
    Simulator(
        const std::shared_ptr<rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses,
        const std::shared_ptr<affine::JointConstraints>& joints,
        const double dt,
        const Settings& settings);

    /// @brief Create a simulator (master constructor).
    Simulator(
        const std::shared_ptr<rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses,
        const std::vector<rigid::Pose>& initial_velocities,
        const std::shared_ptr<affine::JointConstraints>& joints,
        const double dt,
        const Settings& settings);

    ~Simulator();

    /// @brief Run the simulation.
    ///
    /// @param t_end End time
    /// @param callback Callback function to be called at each step.
    ///
    /// The callback function takes a boolean argument indicating whether the
    /// step was successful (i.e., did not fail to converge). The callback is
    /// called after each step, including the final step that reaches t_end.
    ///
    /// If the simulation is already complete (i.e., t >= t_end), the
    /// simulation will not run and the callback will not be called.
    ///
    /// @return True if the simulation ran successfully, false if it was terminated (e.g., due to convergence failure)
    bool run(
        const double t_end,
        const std::function<void(bool)>& callback = [](bool) { });

    /// @brief Step the simulation.
    /// @return True if the step was successful, false if the simulation should be terminated (e.g., due to convergence failure)
    bool step();

    /// @brief Reset the simulation to time t=0.
    void reset();

    // -----------------------------------------------------------------------
    // polysolve::nonlinear::Problem implementation
    //
    // NOTE: When joint constraints are present these operate in the reduced
    // coordinates z (x = Uz); otherwise in the DOFs x directly.
    // -----------------------------------------------------------------------

    /// @brief Compute the incremental potential at v.
    double value(const TVector& v) override;

    /// @brief Compute the gradient of the incremental potential at v.
    void gradient(const TVector& v, TVector& grad) override;

    /// @brief Compute the (sparse) Hessian of the incremental potential at v.
    /// @note Projected to PSD per term unless set_project_to_psd(false).
    void hessian(const TVector& v, THessian& hess) override;

    /// @brief Compute a collision-free maximum step size from v0 to v1 ∈ [0, 1].
    /// Also updates the collision candidates for the swept interval (CCD).
    double max_step_size(const TVector& v0, const TVector& v1) override;

    /// @brief Rebuild the collision set at the line-search start (the
    /// candidates may have changed in max_step_size).
    void line_search_begin(const TVector& v0, const TVector& v1) override;

    /// @brief Rebuild the collision set at the new iterate.
    void solution_changed(const TVector& v) override;

    /// @brief Adaptively update the barrier stiffness after an accepted iterate.
    void post_step(const polysolve::nonlinear::PostStepData& data) override;

    /// @brief Set whether hessian() projects each term to PSD.
    void set_project_to_psd(bool val) override { m_project_to_psd = val; }

    /// @brief Norm of a proposed step: the maximum world-space vertex
    /// displacement it produces (from the current iterate).
    ///
    /// polysolve's Δx stopping criterion uses this norm, so with
    /// Settings::velocity_conv_tol this reproduces the velocity convergence
    /// criterion of Rigid IPC [Ferguson et al. 2021] in a space where a
    /// single scalar tolerance is meaningful (raw DOF-space Δx mixes
    /// translations with rotation-vector radians).
    double step_norm(
        const TVector& delta_v,
        const polysolve::nonlinear::NormType norm_type) const override;

    // -----------------------------------------------------------------------
    // Collision handling
    // -----------------------------------------------------------------------

    /// @brief Update the simulation state before stepping (refreshes the energy terms).
    void initialize_step();

    /// @brief Update the collision sets based on the current state x.
    /// @param x DOF vector.
    /// @param update_candidates If true, also update the candidate collisions based on the current state. If false, only update the normal collision set based on the current candidate collisions.
    void update_collisions(
        Eigen::ConstRef<Eigen::VectorXd> x, bool update_candidates = true);

    /// @brief Capture the vertex-space velocity history for the friction term
    /// (Σᵢ wᵢ V(xᵗ⁻ⁱ) from the time integrator's velocity formula). Called at
    /// the start of every step; kept fixed through all lagging iterations.
    void initialize_friction_step();

    /// @brief Rebuild the tangential (friction) collision set at the lagged
    /// state x (the normal collision set must be up to date at x). The
    /// tangent bases and normal force magnitudes are baked in at x — the
    /// friction lagging [Ferguson et al. 2021].
    void update_friction_collisions(Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Whether any body is KINEMATIC.
    bool has_kinematic_bodies() const;

    /// @brief Attach a kinematic driver to a KINEMATIC body.
    ///
    /// KINEMATIC bodies default to a velocity-driven driver at construction;
    /// call this to script their motion or set a finite drive time. Call after
    /// construction (the body must already be KINEMATIC). The driver is also
    /// captured as the reset() state.
    /// @param body Index of the body (must be KINEMATIC).
    /// @param driver The kinematic driver.
    void
    set_kinematic_driver(const size_t body, const KinematicDriver& driver);

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// @brief History of the (affine) poses at each time step.
    /// @note In rigid mode the rotation part is exactly R(θ).
    const std::list<std::vector<affine::Pose>>& pose_history() const
    {
        return m_pose_history;
    }

    /// @brief History of the poses converted to rigid poses (log map).
    /// @note Exact in rigid mode; in affine mode the rotation part of the
    /// pose must be (numerically) a rotation matrix.
    std::list<std::vector<rigid::Pose>> rigid_pose_history() const;

    /// @brief The current (most recent) affine poses.
    /// @note O(num_bodies); prefer this over pose_history().back() for
    /// per-frame access, since the full history is copied on each access.
    const std::vector<affine::Pose>& poses() const
    {
        return m_pose_history.back();
    }

    /// @brief The current (most recent) poses as rigid poses (log map).
    /// @note O(num_bodies); the per-frame counterpart to rigid_pose_history().
    std::vector<rigid::Pose> rigid_poses() const;

    const rigid::RigidBodies& bodies() const { return *m_bodies; }

    double t() const { return m_t; }

    const Settings& settings() const { return m_settings; }

    VectorMax3d gravity() const { return m_settings.gravity; }
    void set_gravity(Eigen::ConstRef<VectorMax3d> gravity);

    /// @brief Get the barrier potential for collision handling.
    const BarrierPotential& barrier_potential() const
    {
        return *m_barrier_potential;
    }

    /// @brief Get the friction potential (dissipative, velocity-based).
    const FrictionPotential& friction_potential() const
    {
        return *m_friction_potential;
    }

    /// @brief Get the tangential (friction) collision set.
    const TangentialCollisions& tangential_collisions() const
    {
        return *m_tangential_collisions;
    }

    /// @brief The momentum balance ‖∇E + κ∇B + ∇D‖ after the last friction
    /// lagging iteration (-1 before any friction solve).
    double last_momentum_balance() const { return m_last_momentum_balance; }

    /// @brief Get the candidate collisions for the current time step.
    const Candidates& candidates() const { return m_kinematics->candidates(); }

    /// @brief Get the normal collision set for the current time step.
    const NormalCollisions& normal_collisions() const
    {
        return *m_normal_collisions;
    }

    /// @brief Get the DOF model (DOF layout, world map, and CCD strategy).
    const Kinematics& kinematics() const { return *m_kinematics; }

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

    /// @brief Build the time-integrator state from poses and velocities.
    void initialize_time_integrator(
        const std::vector<rigid::Pose>& poses,
        const std::vector<rigid::Pose>& velocities,
        const double dt);

    /// @brief (Re)build the energy terms of the incremental potential.
    /// Called on construction and reset because the terms hold a reference to
    /// the time integrator.
    void initialize_terms();

    /// @brief Refresh the terms at the start of a time step (predicted
    /// positions and body-force wrenches).
    void update_terms();

    /// @brief Compute the initial (adaptive) barrier stiffness for this step.
    void initialize_barrier_stiffness(Eigen::ConstRef<Eigen::VectorXd> x);

    // ---- Joint reduction [Chen et al. 2022, §4.1] --------------------------

    /// @brief Whether joint constraints are active.
    bool has_joints() const;

    /// @brief Map the solver coordinates v to the full DOFs x (x = Uv when
    /// joints are present; identity otherwise).
    Eigen::VectorXd to_full(Eigen::ConstRef<Eigen::VectorXd> v) const;

    // ---- Full-space (DOF) energy evaluation ---------------------------------

    /// @brief Whether the affine terms are engaged (rigid otherwise).
    bool is_affine() const { return m_affine_inertial.has_value(); }

    /// @brief Total incremental potential at the full DOFs x.
    double energy(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Gradient of every term except the barrier, as the terms enter
    /// the objective (i.e., scaled by (βΔt)² where applicable; also used to
    /// initialize the adaptive barrier stiffness).
    Eigen::VectorXd
    non_barrier_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Barrier energy gradient chained to the DOFs: Jᵀ ∇B.
    /// @note Not scaled by the (βΔt)² acceleration scaling; callers apply it.
    Eigen::VectorXd barrier_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Whether the friction term is active.
    bool friction_enabled() const
    {
        return m_settings.friction_coefficient > 0
            && m_settings.friction_iterations != 0;
    }

    /// @brief Per-vertex velocities at x from the time integrator's velocity
    /// formula: (V(x) − Σᵢ wᵢ V(xᵗ⁻ⁱ)) / (βΔt).
    Eigen::MatrixXd
    friction_velocities(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Friction term gradient as it enters the objective:
    /// (βΔt)² Jᵀ ∇ᵥD (zero when friction is disabled or contact-free).
    Eigen::VectorXd
    friction_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const;

    // ---- Kinematic/static bodies -------------------------------------------

    /// @brief Per-body target poses for this step (KINEMATIC bodies: script
    /// front or velocity-integrated; others: current pose, unused).
    std::vector<rigid::Pose> kinematic_targets() const;

    /// @brief Reset the augmented Lagrangian for a new step at state x.
    void al_init(Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Whether the augmented Lagrangian still has unsatisfied channels.
    bool al_active() const;

    /// @brief Apply the AL update policy at the (energy-converged) state x.
    void al_update(Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief AL energy (unscaled by the acceleration scaling).
    double al_energy(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief AL gradient as it enters the objective: (βΔt)² ∇E_AL.
    Eigen::VectorXd al_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Rebuild m_pinned_dofs (STATIC bodies, is_dof_fixed flags, and
    /// AL-satisfied channels of KINEMATIC bodies) in solver coordinates.
    void rebuild_dof_mask();

    /// @brief Advance kinematic scripts and convert expired bodies to STATIC
    /// (called at the end of every step).
    void step_kinematic_bodies();

    /// @brief Total Hessian at the full DOFs x (PSD projected per term if
    /// m_project_to_psd).
    Eigen::SparseMatrix<double>
    full_hessian(Eigen::ConstRef<Eigen::VectorXd> x) const;

    // -----------------------------------------------------------------------

    /// @brief Simulation settings.
    Settings m_settings;

    /// @brief Bodies in the simulation.
    std::shared_ptr<rigid::RigidBodies> m_bodies;

    /// @brief DOF model strategy (rigid or affine parameterization).
    std::shared_ptr<Kinematics> m_kinematics;

    /// @brief Time integrator ([p; vec(Q)] per-body state for both models).
    std::shared_ptr<dynamics::ImplicitTimeIntegrator> m_time_integrator;

    /// @brief Optional linear equality joint constraints (affine mode only).
    std::shared_ptr<affine::JointConstraints> m_joints;

    /// @brief History of poses at each time step.
    std::list<std::vector<affine::Pose>> m_pose_history;

    /// @brief Initial poses (used by reset()).
    std::vector<rigid::Pose> m_initial_poses;

    /// @brief Initial velocities (used by reset()).
    std::vector<rigid::Pose> m_initial_velocities;

    /// @brief Barrier potential for collision handling.
    std::shared_ptr<BarrierPotential> m_barrier_potential;

    /// @brief Normal collision set.
    std::shared_ptr<NormalCollisions> m_normal_collisions;

    /// @brief Friction (dissipative) potential; evaluated on velocities.
    std::shared_ptr<FrictionPotential> m_friction_potential;

    /// @brief Tangential (friction) collision set (lagged: rebuilt between
    /// solves, frozen during each solve).
    std::shared_ptr<TangentialCollisions> m_tangential_collisions;

    /// @brief Vertex-space velocity history Σᵢ wᵢ V(xᵗ⁻ⁱ) — the reference for
    /// the friction velocities (captured at step start, fixed for the step).
    Eigen::MatrixXd m_friction_vertex_history;

    /// @brief Momentum balance after the last friction lagging iteration.
    double m_last_momentum_balance = -1.0;

    // ---- Kinematic/static bodies -------------------------------------------

    /// @brief Augmented Lagrangian for KINEMATIC bodies (mode-dependent).
    std::optional<rigid::AugmentedLagrangian> m_rigid_al;
    std::optional<affine::AugmentedLagrangian> m_affine_al;

    /// @brief Pinned solver-coordinate indices (zero search direction):
    /// STATIC bodies, fixed DOFs, and AL-satisfied kinematic channels.
    std::vector<int> m_pinned_dofs;

    /// @brief Per-body kinematic drivers (empty entry ⇒ not kinematic).
    /// KINEMATIC bodies get a default velocity-driven driver at construction.
    std::vector<std::optional<KinematicDriver>> m_kinematic_drivers;

    /// @brief Body state at construction (restored by reset()).
    std::vector<rigid::RigidBody::Type> m_initial_body_types;
    std::vector<VectorMax6b> m_initial_is_dof_fixed;
    std::vector<std::optional<KinematicDriver>> m_initial_kinematic_drivers;

    /// @brief Broad phase collision handler for efficiently finding potential collisions.
    std::unique_ptr<BroadPhase> m_broad_phase;

    /// @brief Dedicated logger for polysolve (see Settings::solver_log_level).
    /// @note Declared before m_solver: the solver holds a reference to it, so
    /// it must be destroyed after the solver.
    std::shared_ptr<spdlog::logger> m_solver_logger;

    /// @brief polysolve Newton solver used to minimize the incremental potential.
    std::unique_ptr<polysolve::nonlinear::Solver> m_solver;

    // ---- Energy terms (mode-dependent; rebuilt by initialize_terms) --------

    std::optional<rigid::InertialTerm> m_rigid_inertial;
    std::optional<rigid::BodyForces> m_rigid_body_forces;

    std::optional<affine::InertialTerm> m_affine_inertial;
    std::optional<affine::BodyForces> m_affine_body_forces;
    std::optional<affine::OrthogonalityPotential> m_orthogonality;

    /// @brief Whether hessian() projects each term to PSD (set by polysolve
    /// per descent strategy).
    bool m_project_to_psd = true;

    /// @brief The solver's current iterate (tracked through
    /// solution_changed); the base point for step_norm().
    Eigen::VectorXd m_solver_v;

    /// @brief Whether the solver has moved off the starting iterate this
    /// solve; until then step_norm() disables the velocity convergence
    /// criterion so at least one Newton step is always taken (the original
    /// Rigid IPC solver's `iter > 0` guard).
    bool m_solver_stepped = false;

    // -----------------------------------------------------------------------

    /// @brief Current simulation time.
    double m_t = 0.0;

    /// @brief World bounding box diagonal (updated each step).
    double m_bbox_diagonal = 1.0;

    /// @brief Maximum barrier stiffness (from initial_barrier_stiffness).
    double m_max_barrier_stiffness = 0.0;

    /// @brief Minimum distance at the previous Newton iterate (adaptive κ).
    double m_prev_min_distance = -1.0;
};

} // namespace ipc::demo
