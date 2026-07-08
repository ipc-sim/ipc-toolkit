#include "simulator.hpp"

#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/dynamics/affine/joints.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/time_integration/bdf.hpp>
#include <ipc/geometry/normal.hpp> // cross_product_matrix
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/logger.hpp>

#include <polysolve/nonlinear/PostStepData.hpp>
#include <polysolve/nonlinear/Solver.hpp>
#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <cassert>
#include <exception>
#include <limits>
#include <stdexcept>

namespace ipc::demo {

using namespace ipc::dynamics;

namespace {
    // TODO: Why does this not work with 0.5?
    inline constexpr double INFLATION_RADIUS_MULTIPLIER = 1.0;

    std::vector<rigid::Pose>
    zero_velocities(const size_t num_bodies, const int dim)
    {
        return std::vector<rigid::Pose>(num_bodies, rigid::Pose::Identity(dim));
    }

    /// 2D skew "matrix" generator: ω S with S = [[0, -1], [1, 0]].
    inline Eigen::Matrix2d skew_2d(const double omega)
    {
        Eigen::Matrix2d S; // NOLINT
        S << 0, -omega, omega, 0;
        return S;
    }
} // namespace

Simulator::Simulator(
    const std::shared_ptr<rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    const double dt)
    : Simulator(bodies, initial_poses, dt, Settings())
{
}

Simulator::Simulator(
    const std::shared_ptr<rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    const double dt,
    const Settings& settings)
    : Simulator(
          bodies,
          initial_poses,
          zero_velocities(bodies->num_bodies(), bodies->dim()),
          /*joints=*/nullptr,
          dt,
          settings)
{
}

Simulator::Simulator(
    const std::shared_ptr<rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    const std::vector<rigid::Pose>& initial_velocities,
    const double dt,
    const Settings& settings)
    : Simulator(
          bodies,
          initial_poses,
          initial_velocities,
          /*joints=*/nullptr,
          dt,
          settings)
{
}

Simulator::Simulator(
    const std::shared_ptr<rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    const std::shared_ptr<affine::JointConstraints>& joints,
    const double dt)
    : Simulator(bodies, initial_poses, joints, dt, Settings())
{
}

Simulator::Simulator(
    const std::shared_ptr<rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    const std::shared_ptr<affine::JointConstraints>& joints,
    const double dt,
    const Settings& settings)
    : Simulator(
          bodies,
          initial_poses,
          zero_velocities(bodies->num_bodies(), bodies->dim()),
          joints,
          dt,
          settings)
{
}

Simulator::Simulator(
    const std::shared_ptr<rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    const std::vector<rigid::Pose>& initial_velocities,
    const std::shared_ptr<affine::JointConstraints>& joints,
    const double dt,
    const Settings& settings)
    : m_settings(settings)
    , m_bodies(bodies)
    , m_joints(joints)
    , m_initial_poses(initial_poses)
    , m_initial_velocities(initial_velocities)
{
    assert(m_bodies != nullptr);
    assert(initial_poses.size() == m_bodies->num_bodies());
    assert(initial_velocities.size() == m_bodies->num_bodies());

    if (m_settings.gravity.size() != m_bodies->dim()) {
        // The default gravity is 3D; truncate for 2D scenes.
        logger().debug(
            "truncating gravity to the scene dimension ({})", m_bodies->dim());
        m_settings.gravity =
            VectorMax3d(m_settings.gravity.head(m_bodies->dim()));
    }

    if (m_joints != nullptr && !m_joints->empty()
        && m_settings.body_dynamics != BodyDynamics::AFFINE) {
        throw std::invalid_argument(
            "Joint constraints require affine body dynamics "
            "(material-point constraints are nonlinear in the rigid "
            "rotation vectors)");
    }

    switch (m_settings.body_dynamics) {
    case BodyDynamics::RIGID:
        m_kinematics = Kinematics::create_rigid(m_bodies);
        break;
    case BodyDynamics::AFFINE:
        m_kinematics = Kinematics::create_affine(m_bodies);
        break;
    }

    m_pose_history.push_back(
        m_kinematics->poses(m_kinematics->dof(initial_poses)));

    m_barrier_potential = std::make_shared<BarrierPotential>(
        m_settings.dhat, /*stiffness=*/1.0, /*use_physical_barrier=*/true);
    m_normal_collisions = std::make_shared<NormalCollisions>();
    m_normal_collisions->set_use_area_weighting(m_settings.use_area_weighting);
    m_normal_collisions->set_collision_set_type(
        NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX);
    m_broad_phase = std::make_unique<LBVH>();

    if (friction_enabled() && m_settings.static_friction_speed_bound <= 0) {
        throw std::invalid_argument(
            "static_friction_speed_bound must be positive when friction is "
            "enabled");
    }
    // εᵥ is a true speed bound (the friction potential is evaluated on
    // velocities); keep the potential constructible when friction is off.
    m_friction_potential = std::make_shared<FrictionPotential>(
        std::max(m_settings.static_friction_speed_bound, 1e-16));
    m_tangential_collisions = std::make_shared<TangentialCollisions>();

    initialize_time_integrator(initial_poses, initial_velocities, dt);
    initialize_terms();

    if (m_joints != nullptr && !m_joints->empty()) {
        m_joints->finalize();
    }

    // Snapshot the body state so reset() can restore it (kinematic bodies
    // mutate to STATIC on expiry). KINEMATIC bodies get a default
    // velocity-driven driver; set_kinematic_driver() replaces it to script a
    // body or set a finite drive time.
    m_initial_body_types.resize(m_bodies->num_bodies());
    m_initial_is_dof_fixed.resize(m_bodies->num_bodies());
    m_kinematic_drivers.assign(m_bodies->num_bodies(), std::nullopt);
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        const rigid::RigidBody& body = (*m_bodies)[i];
        m_initial_body_types[i] = body.type();
        m_initial_is_dof_fixed[i] = body.is_dof_fixed();
        if (body.type() == rigid::RigidBody::Type::KINEMATIC) {
            m_kinematic_drivers[i] = KinematicDriver::velocity_driven();
        }

        if (is_affine()) {
            // Partial rotational fixing has no per-entry vec(A) analogue.
            const auto rot_fixed = body.is_dof_fixed().tail(
                body.is_dof_fixed().size() - m_bodies->dim());
            if (rot_fixed.any() && !rot_fixed.all()) {
                throw std::invalid_argument(
                    "Partial rotational DOF fixing is not supported in "
                    "affine mode (fix all rotational DOFs or none)");
            }
        }

        // A joint anchored to a prescribed body is over-constrained: the
        // joint RHS is fixed at the initial configuration while the script
        // or pin moves/holds the body.
        if (has_joints() && (!body.is_dynamic() || body.is_dof_fixed().any())
            && m_joints->is_body_constrained(i)) {
            throw std::invalid_argument(
                "Joints may not reference STATIC/KINEMATIC or fixed-DOF "
                "bodies (use add_fixed_point/add_fixed_body instead)");
        }
    }
    m_initial_kinematic_drivers = m_kinematic_drivers;

    rebuild_dof_mask();

    // Create the polysolve Newton solver (the tolerances are rescaled by the
    // characteristic length; use the initial world bounding box diagonal).
    const Eigen::MatrixXd V0 =
        m_kinematics->world_vertices(m_kinematics->dof(initial_poses));
    m_bbox_diagonal =
        (V0.colwise().maxCoeff() - V0.colwise().minCoeff()).norm();
    assert(m_bbox_diagonal > 0);

    // The velocity convergence criterion [Ferguson et al. 2021]: converged
    // when the proposed step moves every vertex slower than the tolerance.
    // The step_norm() override measures steps as world-space vertex
    // displacements, so polysolve's Δx criterion compares against
    // tol × bbox × dt.
    nlohmann::json solver_params = m_settings.solver_params;
    if (m_settings.velocity_conv_tol > 0
        && !solver_params.contains("x_delta_tol")) {
        solver_params["x_delta_tol"] =
            m_settings.velocity_conv_tol * m_bbox_diagonal * dt;
    }

    // polysolve logs through a dedicated logger so its chatter (a per-solve
    // "Finished" report at error level whenever the stop was not the
    // gradient criterion, per-iteration [timing] lines at debug, ...) can be
    // silenced independently of the toolkit's logger. Created without
    // registering it in spdlog's global registry (see ipc/utils/logger.cpp).
    m_solver_logger = std::make_shared<spdlog::logger>(
        "polysolve", std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    m_solver_logger->set_level(m_settings.solver_log_level);

    m_solver = polysolve::nonlinear::Solver::create(
        solver_params, m_settings.linear_solver_params,
        /*characteristic_length=*/m_bbox_diagonal, *m_solver_logger);
}

Simulator::~Simulator() = default;

void Simulator::initialize_time_integrator(
    const std::vector<rigid::Pose>& poses,
    const std::vector<rigid::Pose>& velocities,
    const double dt)
{
    const int dim = m_bodies->dim();

    // The integrator state is affine-shaped [p (dim); vec(Q) (dim²)] per body
    // in both models and dimensions (the rigid optimization DOFs [p; θ] are
    // mapped to/from this state by the parameterization's to-affine map).
    const int pos_ndof = dim;
    const int rot_ndof = dim * dim;
    const int state_ndof = pos_ndof + rot_ndof;

    Eigen::VectorXd x0 =
        Eigen::VectorXd::Zero(state_ndof * m_bodies->num_bodies());
    Eigen::VectorXd v0 =
        Eigen::VectorXd::Zero(state_ndof * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x0.segment(state_ndof * i, dim) = poses[i].position;

        // Zero the velocity components of fixed DOFs.
        rigid::Pose velocity = velocities[i];
        const auto& fixed = (*m_bodies)[i].is_dof_fixed();
        for (int k = 0; k < dim; ++k) {
            if (fixed(k)) {
                velocity.position(k) = 0;
            }
        }
        for (int k = 0; k < velocity.rotation.size(); ++k) {
            if (fixed(dim + k)) {
                velocity.rotation(k) = 0;
            }
        }

        v0.segment(state_ndof * i, dim) = velocity.position;

        if (rot_ndof == 4) {
            // 2D (both models): Q̇ = ω S Q with S the 2D skew generator.
            const Eigen::Matrix2d Q_t0 = poses[i].rotation_matrix();
            x0.segment<4>(state_ndof * i + dim) = Q_t0.reshaped();
            v0.segment<4>(state_ndof * i + dim) =
                (skew_2d(velocity.rotation(0)) * Q_t0).reshaped();
        } else {
            // Q̇ = [ω]× Q with ω the world-frame angular velocity
            const Eigen::Matrix3d Q_t0 = poses[i].rotation_matrix();
            x0.segment<9>(state_ndof * i + 3) = Q_t0.reshaped();
            v0.segment<9>(state_ndof * i + 3) =
                (cross_product_matrix(velocity.rotation) * Q_t0).reshaped();
        }
    }
    const Eigen::VectorXd a0 =
        Eigen::VectorXd::Zero(state_ndof * m_bodies->num_bodies());

    auto bdf = std::make_shared<BDF>(m_settings.bdf_order);
    bdf->set_dt(dt);
    bdf->init(x0, v0, a0, m_bodies->num_bodies(), pos_ndof, rot_ndof);
    m_time_integrator = bdf;
}

// ----------------------------------------------------------------------------
// Energy terms
// ----------------------------------------------------------------------------

void Simulator::initialize_terms()
{
    affine::AugmentedLagrangian::Params al_params;
    al_params.initial_penalty = m_settings.al_initial_penalty;
    al_params.max_penalty = m_settings.al_max_penalty;
    al_params.satisfied_progress = m_settings.al_satisfied_progress;
    al_params.stall_progress = m_settings.al_stall_progress;

    // The orthogonality penalty only applies to the affine (identity) map; the
    // BodyPotentials ignores it for the rigid map.
    m_body_potentials.emplace(
        *m_bodies, m_time_integrator, m_kinematics->to_affine(),
        m_settings.orthogonality_stiffness, al_params);
    m_body_potentials->set_gravity(m_settings.gravity);
}

void Simulator::update_terms() { m_body_potentials->update(*m_bodies); }

void Simulator::set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
{
    m_settings.gravity = gravity;
    if (m_body_potentials) {
        m_body_potentials->set_gravity(gravity);
    }
}

// ----------------------------------------------------------------------------
// Full-space (DOF) energy evaluation
// ----------------------------------------------------------------------------

double Simulator::energy(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    // Every term except the inertia is a pure potential (energy units) and is
    // scaled by (βΔt)² so the incremental potential is uniformly in kg·m²
    // [Lan et al. 2022, Eq. 9]; β changes as the BDF order ramps up, so query
    // the integrator on every evaluation.
    const double s = m_time_integrator->acceleration_scaling();

    // Inertia + s·(body forces + orthogonality + augmented Lagrangian), pulled
    // back to the DOFs through one to-affine map.
    double energy = m_body_potentials->energy(*m_bodies, x, s);

    energy +=
        s
        * (*m_barrier_potential)(
            *m_normal_collisions, *m_bodies, m_kinematics->world_vertices(x));

    if (friction_enabled() && !m_tangential_collisions->empty()) {
        // The friction term (βΔt)²·(βΔt)·D(v(x)) contributes the friction
        // force −μλ to the equations of motion: its x-gradient is
        // (βΔt)²·μλf₁t̂, matching the (βΔt)²-scaled force terms.
        energy +=
            s * m_time_integrator->velocity_scaling()
            * (*m_friction_potential)(
                *m_tangential_collisions, *m_bodies, friction_velocities(x));
    }

    return energy;
}

Eigen::VectorXd
Simulator::non_barrier_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    // The physical body forces only (augmented Lagrangian excluded) so the
    // adaptive barrier stiffness is seeded from the physical forces, as the
    // original Rigid IPC solver does.
    return m_body_potentials->gradient(
        *m_bodies, x, m_time_integrator->acceleration_scaling(),
        /*include_al=*/false);
}

Eigen::VectorXd
Simulator::barrier_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    return m_kinematics->map_vertex_gradient(
        x,
        m_barrier_potential->gradient(
            *m_normal_collisions, *m_bodies, m_kinematics->world_vertices(x)));
}

Eigen::MatrixXd
Simulator::friction_velocities(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    // v(x) = (V(x) − Σᵢ wᵢ V(xᵗ⁻ⁱ)) / (βΔt): the time integrator's velocity
    // formula applied to the world-vertex trajectories.
    assert(m_friction_vertex_history.size() > 0);
    return (m_kinematics->world_vertices(x) - m_friction_vertex_history)
        / m_time_integrator->velocity_scaling();
}

Eigen::VectorXd
Simulator::friction_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    if (!friction_enabled() || m_tangential_collisions->empty()) {
        return Eigen::VectorXd::Zero(x.size());
    }
    // ∇ₓ[(βΔt)³ D(v(x))] = (βΔt)² Jᵀ ∇ᵥD since ∂v/∂V = 1/(βΔt).
    return m_time_integrator->acceleration_scaling()
        * m_kinematics->map_vertex_gradient(
            x,
            m_friction_potential->gradient(
                *m_tangential_collisions, *m_bodies, friction_velocities(x)));
}

Eigen::SparseMatrix<double>
Simulator::full_hessian(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const double s = m_time_integrator->acceleration_scaling();

    // Inertia + s·(body forces + orthogonality + augmented Lagrangian), pulled
    // back to the DOFs through one to-affine map. The rotation block is PSD-
    // projected once (ABS for rigid, CLAMP for affine, matching the historical
    // per-term choices).
    const PSDProjectionMethod body_psd = m_project_to_psd
        ? (is_affine() ? PSDProjectionMethod::CLAMP : PSDProjectionMethod::ABS)
        : PSDProjectionMethod::NONE;
    Eigen::SparseMatrix<double> hess = m_body_potentials->hessian(
        *m_bodies, x, s, body_psd, /*include_al=*/true);

    // Barrier: chain the vertex-space Hessian to the DOFs (Jᵀ H J plus the
    // second-order term for the rigid parameterization, which needs the
    // vertex-space gradient).
    const Eigen::MatrixXd V = m_kinematics->world_vertices(x);
    const Eigen::VectorXd barrier_grad =
        m_barrier_potential->gradient(*m_normal_collisions, *m_bodies, V);
    const Eigen::SparseMatrix<double> barrier_hess =
        m_barrier_potential->hessian(
            *m_normal_collisions, *m_bodies, V,
            m_project_to_psd ? PSDProjectionMethod::CLAMP
                             : PSDProjectionMethod::NONE);

    hess = hess
        + s * m_kinematics->map_vertex_hessian(x, barrier_grad, barrier_hess);

    if (friction_enabled() && !m_tangential_collisions->empty()) {
        // ∇²ₓ[(βΔt)³ D(v(x))] = (βΔt) Jᵀ ∇²ᵥD J + (βΔt)² Σ ∇ᵥD·d²V/dx²
        // (∂v/∂V = 1/(βΔt) is constant); pass the pre-scaled vertex-space
        // gradient/Hessian so map_vertex_hessian applies both pieces.
        const Eigen::MatrixXd velocities = friction_velocities(x);
        const Eigen::VectorXd friction_grad = s
            * m_friction_potential->gradient(
                *m_tangential_collisions, *m_bodies, velocities);
        const Eigen::SparseMatrix<double> friction_hess =
            (s / m_time_integrator->velocity_scaling())
            * m_friction_potential->hessian(
                *m_tangential_collisions, *m_bodies, velocities,
                m_project_to_psd ? PSDProjectionMethod::CLAMP
                                 : PSDProjectionMethod::NONE);
        hess = hess
            + m_kinematics->map_vertex_hessian(x, friction_grad, friction_hess);
    }

    return hess;
}

// ----------------------------------------------------------------------------
// Joint reduction [Chen et al. 2022, §4.1]
// ----------------------------------------------------------------------------

bool Simulator::has_joints() const
{
    return m_joints != nullptr && !m_joints->empty();
}

Eigen::VectorXd Simulator::to_full(Eigen::ConstRef<Eigen::VectorXd> v) const
{
    return has_joints() ? m_joints->to_full(v) : Eigen::VectorXd(v);
}

// ----------------------------------------------------------------------------
// polysolve::nonlinear::Problem implementation
// ----------------------------------------------------------------------------

double Simulator::value(const TVector& v) { return energy(to_full(v)); }

void Simulator::gradient(const TVector& v, TVector& grad)
{
    const Eigen::VectorXd x = to_full(v);
    const double s = m_time_integrator->acceleration_scaling();
    // Body potentials incl. the augmented Lagrangian (non_barrier_gradient
    // excludes it), plus the contact terms via the vertex-space chain rule.
    grad = m_body_potentials->gradient(*m_bodies, x, s, /*include_al=*/true)
        + s * barrier_gradient(x) + friction_gradient(x);
    if (has_joints()) {
        // The pinned entries never move: zero their gradient so the search
        // direction leaves them fixed.
        grad = m_joints->U().transpose() * grad;
        grad.head(m_joints->num_constraints()).setZero();
    }
    // Same treatment for the masked DOFs (STATIC bodies, fixed DOFs, and
    // AL-satisfied kinematic channels).
    for (const int i : m_pinned_dofs) {
        grad(i) = 0;
    }
}

void Simulator::hessian(const TVector& v, THessian& hess)
{
    hess = full_hessian(to_full(v));
    if (has_joints()) {
        // Reduced Hessian UᵀHU with the pinned rows/columns replaced by
        // identity so the direction solve never moves them.
        const Eigen::MatrixXd& U = m_joints->U();
        Eigen::MatrixXd H = U.transpose() * (hess * U);
        const int m = int(m_joints->num_constraints());
        H.topRows(m).setZero();
        H.leftCols(m).setZero();
        H.diagonal().head(m).setOnes();
        for (const int i : m_pinned_dofs) {
            H.row(i).setZero();
            H.col(i).setZero();
            H(i, i) = 1.0;
        }
        hess = H.sparseView();
    } else if (!m_pinned_dofs.empty()) {
        // Zero the pinned rows/columns and put ones on their diagonal so the
        // direction solve never moves them (S H S + (I − S) with a diagonal
        // selector S, staying sparse).
        Eigen::SparseMatrix<double> S(hess.rows(), hess.cols()); // NOLINT
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(hess.rows());
        std::vector<bool> pinned(hess.rows(), false);
        for (const int i : m_pinned_dofs) {
            pinned[i] = true;
        }
        for (int i = 0; i < hess.rows(); ++i) {
            if (!pinned[i]) {
                triplets.emplace_back(i, i, 1.0);
            }
        }
        S.setFromTriplets(triplets.begin(), triplets.end());
        hess = S * hess * S;
        for (const int i : m_pinned_dofs) {
            hess.coeffRef(i, i) = 1.0;
        }
    }
}

double Simulator::max_step_size(const TVector& v0, const TVector& v1)
{
    // CCD filter: also updates the collision candidates for the swept
    // interval.
    return m_kinematics->max_step_size(
        to_full(v0), to_full(v1),
        INFLATION_RADIUS_MULTIPLIER * m_barrier_potential->dhat(),
        m_broad_phase.get());
}

void Simulator::line_search_begin(const TVector& v0, const TVector& /*v1*/)
{
    // The candidate set may have changed in max_step_size; rebuild the
    // collision set at the line-search start so the energy comparisons made
    // by the line search are consistent.
    update_collisions(to_full(v0), /*update_candidates=*/false);
}

void Simulator::solution_changed(const TVector& v)
{
    // The solver has moved off this solve's starting iterate (the first
    // call, from minimize()'s initialization, passes the start point
    // itself; later calls come from line-search trials and accepted steps).
    if (m_solver_v.size() != v.size()
        || (v.array() != m_solver_v.array()).any()) {
        m_solver_stepped = true;
    }

    m_solver_v = v;
    update_collisions(to_full(v), /*update_candidates=*/false);
}

double Simulator::step_norm(
    const TVector& delta_v,
    const polysolve::nonlinear::NormType /*norm_type*/) const
{
    // Guard the velocity convergence criterion until at least one step has
    // been taken this solve (the original Rigid IPC solver's `iter > 0`
    // check): polysolve tests its Δx criterion on the *proposed* direction
    // before the first step, so a quasi-static scene whose per-step motion
    // (~g·dt²) is below the tolerance would otherwise freeze forever — the
    // step is never taken, and the integrator then re-derives zero velocity.
    if (!m_solver_stepped) {
        return std::numeric_limits<double>::infinity();
    }
    if (m_solver_v.size() != delta_v.size()) {
        return delta_v.norm(); // fallback: no base point available
    }
    // Maximum world-space vertex displacement of the proposed step.
    return (m_kinematics->world_vertices(to_full(m_solver_v + delta_v))
            - m_kinematics->world_vertices(to_full(m_solver_v)))
        .lpNorm<Eigen::Infinity>();
}

void Simulator::post_step(const polysolve::nonlinear::PostStepData& data)
{
    const Eigen::VectorXd x = to_full(data.x);

    const double min_distance = m_normal_collisions->compute_minimum_distance(
        *m_bodies, m_kinematics->world_vertices(x));

    const double kappa = update_barrier_stiffness(
        m_prev_min_distance, min_distance, m_max_barrier_stiffness,
        m_barrier_potential->stiffness(), m_bbox_diagonal,
        m_settings.dhat_epsilon_scale);

    if (kappa != m_barrier_potential->stiffness()) {
        logger().debug(
            "updated barrier stiffness: κ={:g} min_distance={:g}", kappa,
            min_distance);
        m_barrier_potential->set_stiffness(kappa);
    }

    m_prev_min_distance = min_distance;
}

// ----------------------------------------------------------------------------
// Stepping
// ----------------------------------------------------------------------------

bool Simulator::run(
    const double t_end, const std::function<void(bool)>& callback)
{
    if (!has_time_remaining(m_time_integrator->dt(), t_end)) {
        logger().warn(
            "simulation already complete: t={:g} t_end={:g}", m_t, t_end);
        return false; // Simulation already complete
    }

    bool step_succeeded = true;
    while (step_succeeded
           && has_time_remaining(m_time_integrator->dt(), t_end)) {
        step_succeeded = step();
        callback(step_succeeded);
    }
    return step_succeeded;
}

void Simulator::initialize_step() { update_terms(); }

void Simulator::initialize_barrier_stiffness(Eigen::ConstRef<Eigen::VectorXd> x)
{
    const Eigen::MatrixXd V = m_kinematics->world_vertices(x);
    m_bbox_diagonal = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
    assert(m_bbox_diagonal > 0);

    double average_mass = 0;
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        average_mass += (*m_bodies)[i].mass();
    }
    average_mass /= m_bodies->num_bodies();

    // Sum the gradients of every term except the barrier (as they enter the
    // objective, i.e., scaled by (βΔt)² where applicable).
    const Eigen::VectorXd grad_energy = non_barrier_gradient(x);

    // initial_barrier_stiffness works in the units of the *raw* barrier, but
    // the barrier enters the objective scaled by (a) the potential's physical
    // scaling dhat/units(dhat²) (if enabled) and (b) the (βΔt)² acceleration
    // scaling applied by the simulator. Remove the stiffness and physical
    // scaling from the gradient (both are linear), and convert the resulting
    // stiffness back to the potential's units accounting for both scales.
    const double dhat = m_barrier_potential->dhat();
    const double physical_scale = m_barrier_potential->use_physical_barrier()
        ? (dhat / m_barrier_potential->barrier().units(dhat * dhat))
        : 1.0;
    const double barrier_scale =
        m_time_integrator->acceleration_scaling() * physical_scale;

    const Eigen::VectorXd grad_barrier = barrier_gradient(x)
        / (m_barrier_potential->stiffness() * physical_scale);

    const double kappa = initial_barrier_stiffness(
        m_bbox_diagonal, m_barrier_potential->barrier(), dhat, average_mass,
        grad_energy, grad_barrier, m_max_barrier_stiffness,
        m_settings.min_barrier_stiffness_scale);
    m_barrier_potential->set_stiffness(kappa / barrier_scale);
    m_max_barrier_stiffness /= barrier_scale;

    m_prev_min_distance =
        m_normal_collisions->compute_minimum_distance(*m_bodies, V);

    logger().debug(
        "initial barrier stiffness: κ={:g} κ_max={:g} min_distance={:g}",
        m_barrier_potential->stiffness(), m_max_barrier_stiffness,
        m_prev_min_distance);
}

bool Simulator::step()
{
    // Convert expired kinematic bodies to STATIC before the terms refresh
    // (the affine mass matrix depends on the types).
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        if (m_kinematic_drivers[i].has_value()
            && m_kinematic_drivers[i]->expired(m_time_integrator->dt())) {
            logger().debug(
                "kinematic body {} expired; converting to static", i);
            (*m_bodies)[i].convert_to_static();
            m_kinematic_drivers[i].reset();
        }
    }

    initialize_step();

    // For the rigid model this maps the rotation matrices back to canonical
    // rotation vectors (‖θ‖ ≤ π), keeping the optimization variable bounded.
    Eigen::VectorXd x =
        m_kinematics->from_integrator_state(m_time_integrator->x_prev(0));

    update_collisions(x, /*update_candidates=*/true);
    initialize_barrier_stiffness(x);

    // Reset the augmented Lagrangian toward this step's kinematic targets and
    // rebuild the DOF mask (types may have changed above).
    al_init(x);
    rebuild_dof_mask();

    // Friction lagging convergence tolerance [Ferguson et al. 2021]; fixed
    // for the whole step (initialize_barrier_stiffness refreshes the bbox).
    const double eps_d = 1e-2 * m_bbox_diagonal;
    m_last_momentum_balance = -1;
    initialize_friction_step();
    update_friction_collisions(x);

    // Solve in the reduced coordinates z = Vx with the constrained entries
    // pinned when joints are present (this also projects x exactly onto the
    // constraints).
    Eigen::VectorXd v = x;
    if (has_joints()) {
        v = m_joints->to_reduced(x);
        v.head(m_joints->num_constraints()) = m_joints->rhs();
    }

    // The unified outer loop: each iteration is a full Newton solve followed
    // by (a) an augmented Lagrangian update while kinematic targets are
    // unsatisfied and (b) a friction re-lag (rebuild the collision and
    // friction sets at the new state) until the momentum balance
    // ‖∇E + κ∇B + ∇D‖ ≤ eps_d [Ferguson et al. 2021].
    bool converged = true;
    double prev_momentum_balance = std::numeric_limits<double>::infinity();
    for (int outer_iter = 0;; ++outer_iter) {
        m_solver_v = v; // base point for step_norm() until solution_changed
        m_solver_stepped = false; // arm the velocity criterion after one step

        try {
            m_solver->minimize(*this, v);
        } catch (const std::exception& e) {
            // polysolve throws on convergence failures (line search failure,
            // iteration limit, non-finite energy, ...).
            logger().warn("solve failed: {}", e.what());
            converged = false;
        }
        x = to_full(v);

        if (!converged) {
            break;
        }

        bool needs_more = false;

        // Augmented Lagrangian channel: update the penalties/multipliers at
        // the energy-converged state; satisfied channels freeze their DOFs.
        if (al_active()) {
            al_update(x);
            rebuild_dof_mask();
            needs_more = al_active();
        }

        // Friction channel: re-lag and measure how far the solution is from
        // a fixed point.
        if (friction_enabled()) {
            update_collisions(x, /*update_candidates=*/true);
            update_friction_collisions(x);

            Eigen::VectorXd grad;
            gradient(v, grad);
            m_last_momentum_balance = grad.norm();
            logger().debug(
                "friction lagging: iteration={:d} momentum_balance={:g} "
                "eps_d={:g}",
                outer_iter + 1, m_last_momentum_balance, eps_d);

            // Stop on convergence, on the iteration cap (normal termination,
            // as in the original Rigid IPC), or when the momentum balance
            // stalls: the inner solver's own stopping criteria bound the
            // reachable gradient residual (e.g., a resting contact held by
            // static friction stops on the velocity criterion with an
            // O(force) residual), so require ≥10% improvement per lagging
            // iteration to keep going.
            const bool improving =
                m_last_momentum_balance < 0.9 * prev_momentum_balance;
            prev_momentum_balance = m_last_momentum_balance;

            const bool friction_needs_more = m_last_momentum_balance > eps_d
                && improving
                && (m_settings.friction_iterations < 0
                    || outer_iter + 1 < m_settings.friction_iterations);
            if (!friction_needs_more && m_last_momentum_balance > eps_d) {
                logger().debug(
                    "friction lagging ended above tolerance: "
                    "momentum_balance={:g} eps_d={:g} iterations={:d}",
                    m_last_momentum_balance, eps_d, outer_iter + 1);
            }
            needs_more = needs_more || friction_needs_more;
        } else if (needs_more) {
            // The AL wants another solve: refresh the collision set at the
            // new state (the friction branch above does this otherwise).
            update_collisions(x, /*update_candidates=*/true);
        }

        if (needs_more && outer_iter + 1 >= m_settings.max_outer_iterations) {
            if (al_active()) {
                // A kinematic body genuinely failed to reach its target.
                logger().warn(
                    "augmented Lagrangian did not satisfy the kinematic "
                    "targets in {} outer iterations",
                    m_settings.max_outer_iterations);
                converged = false;
            }
            needs_more = false;
        }
        if (!needs_more) {
            break;
        }

        // Re-seed the barrier stiffness for the next solve (the original
        // Rigid IPC solver re-initializes κ at the start of every solve).
        initialize_barrier_stiffness(x);
    }

    if (!converged) {
        logger().warn("step failed to converge: t={:g}", m_t);
        if (m_settings.abort_on_convergence_failure) {
            return false;
        }
    }

    m_pose_history.push_back(m_kinematics->poses(x));

    m_time_integrator->update(m_kinematics->to_integrator_state(x));

    m_t += m_time_integrator->dt();

    // Advance the kinematic scripts (pop scripted poses, count down the
    // drive time).
    step_kinematic_bodies();

    return true;
}

void Simulator::update_collisions(
    Eigen::ConstRef<Eigen::VectorXd> x, bool update_candidates)
{
    const Eigen::MatrixXd V = m_kinematics->world_vertices(x);

    if (update_candidates) {
        m_kinematics->update_candidates(
            x, INFLATION_RADIUS_MULTIPLIER * m_barrier_potential->dhat(),
            m_broad_phase.get());
    }

    m_normal_collisions->build(
        m_kinematics->candidates(), *m_bodies, V, m_barrier_potential->dhat());
}

void Simulator::initialize_friction_step()
{
    if (!friction_enabled()) {
        return;
    }

    // Σᵢ wᵢ V(xᵗ⁻ⁱ): the history part of the integrator's velocity formula
    // applied to the world-vertex trajectories.
    const std::vector<double> weights =
        m_time_integrator->velocity_history_weights();
    assert(!weights.empty());
    m_friction_vertex_history = weights[0]
        * m_kinematics->world_vertices(
            m_kinematics->from_integrator_state(m_time_integrator->x_prev(0)));
    for (size_t i = 1; i < weights.size(); ++i) {
        m_friction_vertex_history += weights[i]
            * m_kinematics->world_vertices(m_kinematics->from_integrator_state(
                m_time_integrator->x_prev(i)));
    }
}

void Simulator::update_friction_collisions(Eigen::ConstRef<Eigen::VectorXd> x)
{
    if (!friction_enabled()) {
        return;
    }

    // The normal force magnitudes (κ through the barrier potential) and the
    // tangent bases are baked in at x and frozen until the next rebuild.
    m_tangential_collisions->build(
        *m_bodies, m_kinematics->world_vertices(x), *m_normal_collisions,
        *m_barrier_potential, m_settings.friction_coefficient);
}

// ----------------------------------------------------------------------------
// Kinematic/static bodies
// ----------------------------------------------------------------------------

bool Simulator::has_kinematic_bodies() const
{
    return m_bodies->count_type(rigid::RigidBody::Type::KINEMATIC) > 0;
}

std::vector<rigid::Pose> Simulator::kinematic_targets() const
{
    // Start from the current poses (non-KINEMATIC entries are unused).
    std::vector<rigid::Pose> targets(m_bodies->num_bodies());

    const int dim = m_bodies->dim();
    const double h = m_time_integrator->dt();
    const Eigen::VectorXd& v_prev = m_time_integrator->v_prev(0);
    const int state_ndof =
        int(m_time_integrator->pos_ndof() + m_time_integrator->rot_ndof());

    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        const affine::Pose pose = m_time_integrator->pose(i);

        // The current pose (position + rotation vector/angle).
        rigid::Pose current;
        current.position = pose.position;
        current.rotation = pose.rotation_vector();

        targets[i] = current;

        if (!m_kinematic_drivers[i].has_value()) {
            continue;
        }

        // The prescribed velocity: linear in .position, angular ω in
        // .rotation. The driver applies the script-or-integrate policy; the
        // ω extraction from the integrator state layout lives here.
        rigid::Pose velocity;
        velocity.position = v_prev.segment(state_ndof * i, dim);
        if (dim == 2) {
            // 2D (both models): ω from Q̇ = ω S Q ⇒ ω = (Q̇ Qᵀ)(1, 0).
            const Eigen::Matrix2d Q = pose.rotation;
            const Eigen::Matrix2d Qdot =
                v_prev.segment<4>(state_ndof * i + 2).reshaped(2, 2);
            velocity.rotation =
                VectorMax3d::Constant(1, (Qdot * Q.transpose())(1, 0));
        } else {
            // 3D: [ω]× = Q̇ Qᵀ (projected to antisymmetric for robustness).
            const Eigen::Matrix3d Q = pose.rotation;
            const Eigen::Matrix3d Qdot =
                v_prev.segment<9>(state_ndof * i + 3).reshaped(3, 3);
            const Eigen::Matrix3d omega_hat = 0.5
                * (Qdot * Q.transpose() - (Qdot * Q.transpose()).transpose());
            velocity.rotation = Eigen::Vector3d(
                omega_hat(2, 1), omega_hat(0, 2), omega_hat(1, 0));
        }

        targets[i] = m_kinematic_drivers[i]->target(current, velocity, h);
    }
    return targets;
}

void Simulator::al_init(Eigen::ConstRef<Eigen::VectorXd> x)
{
    const std::vector<rigid::Pose> targets = kinematic_targets();
    m_body_potentials->init_augmented_lagrangian(*m_bodies, x, targets);
}

bool Simulator::al_active() const
{
    return m_body_potentials->augmented_lagrangian_active();
}

void Simulator::al_update(Eigen::ConstRef<Eigen::VectorXd> x)
{
    m_body_potentials->update_augmented_lagrangian(*m_bodies, x);
}

void Simulator::rebuild_dof_mask()
{
    m_pinned_dofs.clear();

    const int dim = m_bodies->dim();
    const int body_ndof = is_affine() ? (dim + dim * dim) : (dim == 2 ? 3 : 6);
    const bool linear_satisfied =
        m_body_potentials->augmented_lagrangian_linear_satisfied();
    const bool angular_satisfied =
        m_body_potentials->augmented_lagrangian_angular_satisfied();

    std::vector<bool> pinned(body_ndof * m_bodies->num_bodies(), false);
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        const rigid::RigidBody& body = (*m_bodies)[i];
        const auto& fixed = body.is_dof_fixed();
        const int pos_ndof = dim, rot_flags = fixed.size() - pos_ndof;
        const bool kinematic = body.type() == rigid::RigidBody::Type::KINEMATIC;
        const bool is_static = body.type() == rigid::RigidBody::Type::STATIC;

        // Positional DOFs (world axes; identical layout in both modes).
        for (int k = 0; k < pos_ndof; ++k) {
            if (is_static || fixed(k) || (kinematic && linear_satisfied)) {
                pinned[body_ndof * i + k] = true;
            }
        }

        // Rotational DOFs.
        const bool pin_rotation = is_static || (kinematic && angular_satisfied)
            || fixed.tail(rot_flags).all();
        if (is_affine()) {
            // Whole-rotation only (partial rejected at construction).
            if (pin_rotation) {
                for (int k = dim; k < body_ndof; ++k) {
                    pinned[body_ndof * i + k] = true;
                }
            }
        } else {
            for (int k = 0; k < rot_flags; ++k) {
                if (pin_rotation || fixed(pos_ndof + k)) {
                    pinned[body_ndof * i + pos_ndof + k] = true;
                }
            }
        }
    }

    // Translate to solver coordinates (identity without joints; through the
    // reduction's identity-on-unjointed-DOF property with joints).
    for (size_t j = 0; j < pinned.size(); ++j) {
        if (!pinned[j]) {
            continue;
        }
        m_pinned_dofs.push_back(
            has_joints() ? m_joints->free_reduced_index(int(j)) : int(j));
    }
}

void Simulator::step_kinematic_bodies()
{
    const double dt = m_time_integrator->dt();
    for (auto& driver : m_kinematic_drivers) {
        if (driver.has_value()) {
            driver->step(dt);
        }
    }
}

void Simulator::set_kinematic_driver(
    const size_t body, const KinematicDriver& driver)
{
    if ((*m_bodies)[body].type() != rigid::RigidBody::Type::KINEMATIC) {
        throw std::invalid_argument(
            "set_kinematic_driver requires the body to be KINEMATIC");
    }
    // Also update the reset() snapshot so reset() restores this driver.
    m_kinematic_drivers[body] = driver;
    m_initial_kinematic_drivers[body] = driver;
}

// ----------------------------------------------------------------------------
// Pose history
// ----------------------------------------------------------------------------

namespace {
    std::vector<rigid::Pose>
    to_rigid_poses(const std::vector<affine::Pose>& affine_poses)
    {
        std::vector<rigid::Pose> rigid_poses(affine_poses.size());
        for (size_t i = 0; i < affine_poses.size(); ++i) {
            rigid_poses[i].position = affine_poses[i].position;
            rigid_poses[i].rotation = affine_poses[i].rotation_vector();
        }
        return rigid_poses;
    }
} // namespace

std::vector<rigid::Pose> Simulator::rigid_poses() const
{
    return to_rigid_poses(m_pose_history.back());
}

std::list<std::vector<rigid::Pose>> Simulator::rigid_pose_history() const
{
    std::list<std::vector<rigid::Pose>> history;
    for (const auto& poses : m_pose_history) {
        history.push_back(to_rigid_poses(poses));
    }
    return history;
}

void Simulator::reset()
{
    m_t = 0.0;
    // Reset pose history to only contain the initial poses
    m_pose_history.resize(1);

    // Rebuild the time integrator and the energy terms (the terms hold
    // references to the integrator)
    initialize_time_integrator(
        m_initial_poses, m_initial_velocities, m_time_integrator->dt());
    initialize_terms();

    m_kinematics->clear_candidates();
    m_normal_collisions->clear();
    m_tangential_collisions->clear();
    m_friction_vertex_history.resize(0, 0);
    m_last_momentum_balance = -1;

    // Restore the body state and kinematic drivers (kinematic bodies mutate
    // to STATIC on expiry).
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        rigid::RigidBody& body = (*m_bodies)[i];
        body.set_type(m_initial_body_types[i]);
        body.set_is_dof_fixed(m_initial_is_dof_fixed[i]);
    }
    m_kinematic_drivers = m_initial_kinematic_drivers;
    rebuild_dof_mask();
}

} // namespace ipc::demo
