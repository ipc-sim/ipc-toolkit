#include "simulator.hpp"

#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/dynamics/affine/joints.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/time_integration/bdf.hpp>
#include <ipc/geometry/normal.hpp> // cross_product_matrix
#include <ipc/potentials/barrier_potential.hpp>
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

    std::vector<rigid::Pose> zero_velocities(const size_t num_bodies)
    {
        return std::vector<rigid::Pose>(
            num_bodies,
            rigid::Pose(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()));
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
          zero_velocities(bodies->num_bodies()),
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
          zero_velocities(bodies->num_bodies()),
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
    assert(m_bodies->dim() == 3); // TODO: Support 2D bodies
    assert(initial_poses.size() == m_bodies->num_bodies());
    assert(initial_velocities.size() == m_bodies->num_bodies());

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

    initialize_time_integrator(initial_poses, initial_velocities, dt);
    initialize_terms();

    if (m_joints != nullptr && !m_joints->empty()) {
        m_joints->finalize();
    }

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
        "polysolve",
        std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
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
    // The integrator state is affine-shaped ([p (3); vec(Q) column-major (9)]
    // per body) for both dynamics models.
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        const Eigen::Matrix3d Q_t0 = poses[i].rotation_matrix();
        x0.segment<3>(12 * i) = poses[i].position;
        x0.segment<9>(12 * i + 3) = Q_t0.reshaped();

        // Q̇ = [ω]× Q with ω the world-frame angular velocity
        v0.segment<3>(12 * i) = velocities[i].position;
        v0.segment<9>(12 * i + 3) =
            (cross_product_matrix(velocities[i].rotation) * Q_t0).reshaped();
    }
    const Eigen::VectorXd a0 =
        Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());

    auto bdf = std::make_shared<BDF>(m_settings.bdf_order);
    bdf->set_dt(dt);
    bdf->init(x0, v0, a0, m_bodies->num_bodies());
    m_time_integrator = bdf;
}

// ----------------------------------------------------------------------------
// Energy terms
// ----------------------------------------------------------------------------

void Simulator::initialize_terms()
{
    m_rigid_inertial.reset();
    m_rigid_body_forces.reset();
    m_affine_inertial.reset();
    m_affine_body_forces.reset();
    m_orthogonality.reset();

    switch (m_settings.body_dynamics) {
    case BodyDynamics::RIGID:
        m_rigid_inertial.emplace(m_time_integrator);
        m_rigid_body_forces.emplace(m_time_integrator);
        m_rigid_body_forces->set_gravity(m_settings.gravity);
        break;
    case BodyDynamics::AFFINE:
        m_affine_inertial.emplace(*m_bodies, m_time_integrator);
        m_affine_body_forces.emplace(m_time_integrator);
        m_affine_body_forces->set_gravity(m_settings.gravity);
        m_orthogonality.emplace(m_settings.orthogonality_stiffness);
        break;
    }
}

void Simulator::update_terms()
{
    if (is_affine()) {
        m_affine_inertial->update(*m_bodies);
        m_affine_body_forces->update(*m_bodies);
    } else {
        m_rigid_inertial->update(*m_bodies);
        m_rigid_body_forces->update(*m_bodies);
    }
}

void Simulator::set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
{
    m_settings.gravity = gravity;
    if (m_rigid_body_forces) {
        m_rigid_body_forces->set_gravity(gravity);
    }
    if (m_affine_body_forces) {
        m_affine_body_forces->set_gravity(gravity);
    }
}

// ----------------------------------------------------------------------------
// Full-space (DOF) energy evaluation
// ----------------------------------------------------------------------------

double Simulator::energy(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    double energy;
    if (is_affine()) {
        // The orthogonality potential enters the incremental potential
        // scaled by (βΔt)² [Lan et al. 2022, Eq. 9]; β changes as the BDF
        // order ramps up, so query the integrator on every evaluation.
        energy = (*m_affine_inertial)(*m_bodies, x)
            + (*m_affine_body_forces)(*m_bodies, x)
            + m_time_integrator->acceleration_scaling()
                * (*m_orthogonality)(*m_bodies, x);
    } else {
        energy = (*m_rigid_inertial)(*m_bodies, x)
            + (*m_rigid_body_forces)(*m_bodies, x);
    }

    return energy
        + (*m_barrier_potential)(
               *m_normal_collisions, *m_bodies, m_kinematics->world_vertices(x));
}

Eigen::VectorXd
Simulator::non_barrier_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    if (is_affine()) {
        return m_affine_inertial->gradient(*m_bodies, x)
            + m_affine_body_forces->gradient(*m_bodies, x)
            + m_time_integrator->acceleration_scaling()
            * m_orthogonality->gradient(*m_bodies, x);
    } else {
        return m_rigid_inertial->gradient(*m_bodies, x)
            + m_rigid_body_forces->gradient(*m_bodies, x);
    }
}

Eigen::VectorXd
Simulator::barrier_gradient(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    return m_kinematics->map_vertex_gradient(
        x,
        m_barrier_potential->gradient(
            *m_normal_collisions, *m_bodies, m_kinematics->world_vertices(x)));
}

Eigen::SparseMatrix<double>
Simulator::full_hessian(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    // NOTE: The body-force energies are linear in the DOFs (zero Hessian).
    Eigen::SparseMatrix<double> hess(x.size(), x.size());
    if (is_affine()) {
        // The affine inertial Hessian is the constant mass matrix (PSD by
        // construction).
        hess = m_affine_inertial->mass_matrix()
            + Eigen::SparseMatrix<double>(
                   m_time_integrator->acceleration_scaling()
                   * m_orthogonality->hessian(
                       *m_bodies, x,
                       m_project_to_psd ? PSDProjectionMethod::CLAMP
                                        : PSDProjectionMethod::NONE));
    } else {
        hess = m_rigid_inertial->hessian(
            *m_bodies, x,
            m_project_to_psd ? PSDProjectionMethod::ABS
                             : PSDProjectionMethod::NONE);
    }

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

    return hess
        + m_kinematics->map_vertex_hessian(x, barrier_grad, barrier_hess);
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
    grad = non_barrier_gradient(x) + barrier_gradient(x);
    if (has_joints()) {
        // The pinned entries never move: zero their gradient so the search
        // direction leaves them fixed.
        grad = m_joints->U().transpose() * grad;
        grad.head(m_joints->num_constraints()).setZero();
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
        hess = H.sparseView();
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

    // Sum the gradients of every term except the barrier
    const Eigen::VectorXd grad_energy = non_barrier_gradient(x);

    // initial_barrier_stiffness works in the units of the *raw* barrier, but
    // the potential may scale the barrier to physical units by
    // dhat/units(dhat²). Remove the stiffness and physical scaling from the
    // gradient (both are linear), and convert the resulting stiffness back to
    // the potential's units.
    const double dhat = m_barrier_potential->dhat();
    const double physical_scale = m_barrier_potential->use_physical_barrier()
        ? (dhat / m_barrier_potential->barrier().units(dhat * dhat))
        : 1.0;

    const Eigen::VectorXd grad_barrier = barrier_gradient(x)
        / (m_barrier_potential->stiffness() * physical_scale);

    const double kappa = initial_barrier_stiffness(
        m_bbox_diagonal, m_barrier_potential->barrier(), dhat, average_mass,
        grad_energy, grad_barrier, m_max_barrier_stiffness,
        m_settings.min_barrier_stiffness_scale);
    m_barrier_potential->set_stiffness(kappa / physical_scale);
    m_max_barrier_stiffness /= physical_scale;

    m_prev_min_distance =
        m_normal_collisions->compute_minimum_distance(*m_bodies, V);

    logger().debug(
        "initial barrier stiffness: κ={:g} κ_max={:g} min_distance={:g}",
        m_barrier_potential->stiffness(), m_max_barrier_stiffness,
        m_prev_min_distance);
}

bool Simulator::step()
{
    initialize_step();

    // For the rigid model this maps the rotation matrices back to canonical
    // rotation vectors (‖θ‖ ≤ π), keeping the optimization variable bounded.
    Eigen::VectorXd x =
        m_kinematics->from_integrator_state(m_time_integrator->x_prev(0));

    update_collisions(x, /*update_candidates=*/true);
    initialize_barrier_stiffness(x);

    // Solve in the reduced coordinates z = Vx with the constrained entries
    // pinned when joints are present (this also projects x exactly onto the
    // constraints).
    Eigen::VectorXd v = x;
    if (has_joints()) {
        v = m_joints->to_reduced(x);
        v.head(m_joints->num_constraints()) = m_joints->rhs();
    }
    m_solver_v = v; // base point for step_norm() until solution_changed
    m_solver_stepped = false; // arm the velocity criterion after one step

    bool converged = true;
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
        logger().warn("step failed to converge: t={:g}", m_t);
        if (m_settings.abort_on_convergence_failure) {
            return false;
        }
    }

    m_pose_history.push_back(m_kinematics->poses(x));

    m_time_integrator->update(m_kinematics->to_integrator_state(x));

    m_t += m_time_integrator->dt();

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

// ----------------------------------------------------------------------------
// Pose history
// ----------------------------------------------------------------------------

namespace {
    std::vector<rigid::Pose>
    to_rigid_poses(const std::vector<rigid::AffinePose>& affine_poses)
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
}

} // namespace ipc::demo
