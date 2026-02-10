#include "simulator.hpp"

#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/ccd/additive_ccd.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/dynamics/rigid/body_forces.hpp>
#include <ipc/dynamics/rigid/ground_contact.hpp>
#include <ipc/dynamics/rigid/inertial_term.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/time_integrator.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <igl/PI.h>
#include <igl/writePLY.h>

namespace ipc::rigid {

Simulator::Simulator(
    const std::shared_ptr<RigidBodies>& _bodies,
    const std::vector<Pose>& initial_poses,
    const double dt)
    : m_bodies(_bodies)
{
    assert(initial_poses.size() == m_bodies->num_bodies());
    m_pose_history.push_back(initial_poses);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x0.segment<3>(12 * i) = initial_poses[i].position;
        x0.segment<9>(12 * i + 3) =
            initial_poses[i].rotation_matrix().reshaped();
    }
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    // for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
    //     // v0.segment<3>(12 * i).y() = 10;

    //     // ω = R₀ᵀω₀ (ω₀ expressed in body coordinates)
    //     Eigen::Vector3d omega(0, -100 * igl::PI / 180, 0);
    //     omega = (*m_bodies)[i].R0().transpose() * omega;
    //     Eigen::Matrix3d Q_t0 = initial_poses[i].rotation_matrix();
    //     v0.segment<9>(12 * i + 3) =
    //         (Q_t0 * cross_product_matrix(omega)).reshaped();
    // }
    Eigen::VectorXd a0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());

    // Initialize the time integrator
    m_time_integrator =
        std::make_shared<ImplicitEuler>(x0, v0, a0, dt, m_bodies->num_bodies());

    // Initialize the inertial term
    m_inertial_term = std::make_unique<InertialTerm>(m_time_integrator);

    // Iniaizlie the body force term
    m_body_forces = std::make_unique<BodyForces>(m_time_integrator);
    m_body_forces->set_gravity(Eigen::Vector3d(0, -9.81, 0));

    // Initialize the ground contact handler
    m_ground_contact = std::make_unique<GroundContact>(0.0);
    m_ground_contact->set_dhat(0.1);

    const double dhat = 0.1, kappa = 10.0;
    m_barrier_potential = std::make_unique<BarrierPotential>(dhat, kappa, true);
    m_normal_collisions = std::make_unique<NormalCollisions>();
    m_candidates = std::make_unique<Candidates>();
    m_broad_phase = std::make_unique<LBVH>();
    m_additive_ccd = std::make_unique<AdditiveCCD>();
}

Simulator::~Simulator() = default;

bool Simulator::run(
    // const double dt,
    const double t_end,
    const std::function<void(bool)>& callback)
{
    if (!has_time_remaining(m_time_integrator->dt, t_end)) {
        logger().warn(
            "simulation already complete: t={:g} t_end={:g}", m_t, t_end);
        return false; // Simulation already complete
    }

    // m_time_integrator->dt = dt;

    bool step_succeeded = true;
    while (step_succeeded && has_time_remaining(m_time_integrator->dt, t_end)) {
        step_succeeded = step();
        callback(step_succeeded);
    }
    return step_succeeded;
}

void Simulator::initialize_step()
{
    m_inertial_term->update(*m_bodies);
    m_body_forces->update(*m_bodies);
}

bool Simulator::step()
{
    std::vector<Pose> poses = m_pose_history.back();

    initialize_step();

    Eigen::VectorXd x = Pose::from_poses(m_pose_history.back());

    update_collisions(x, true);

    double dx, grad_norm;
    int iter = 0;
    do {
        Eigen::VectorXd grad = gradient(x);
        if ((grad_norm = grad.norm()) < 1e-6) {
            break;
        }
        Eigen::MatrixXd hess = hessian(x);
        Eigen::VectorXd step = -hess.llt().solve(grad);
        if (step.dot(grad) > 0) {
            logger().warn(
                "not a descent direction: step.dot(grad)={:g} grad.norm()={:g} hess.norm()={:g}",
                step.dot(grad), grad_norm, hess.norm());
            step = -grad;
        }
        dx = step.norm();
        double alpha = 1.0;
        double Ex = energy(x);
        logger().trace(
            "iter={} dx={:g} norm(grad)={:g} norm(hess)={:g} energy={:g}", iter,
            dx, grad_norm, hess.norm(), Ex);
        update_candidates(x, x + alpha * step);
        alpha = m_candidates->compute_collision_free_stepsize(
            *m_bodies, m_bodies->vertices(x),
            m_bodies->vertices(x + alpha * step), 0.0, *m_additive_ccd);
        logger().trace(
            "initial alpha={:g} energy={:g} num_candidates={}", alpha, Ex,
            m_candidates->size());
        update_collisions(x + alpha * step, false);
        while (energy(x + alpha * step) >= Ex) {
            alpha *= 0.5;
            update_collisions(x + alpha * step, false);
            logger().trace(
                "line search: alpha={:g} energy={:g} num_collisions={}", alpha,
                energy(x + alpha * step), m_normal_collisions->size());
            if (alpha < 1e-12) {
                logger().error(
                    "line search failed: alpha={:g} Ex={:g} E(x+alpha*step)={:g}",
                    alpha, Ex, energy(x + alpha * step));
                return false;
            }
        }
        x += alpha * step;
        double min_distance = m_normal_collisions->compute_minimum_distance(
            *m_bodies, m_bodies->vertices(x));
        logger().debug(
            "step: dx={:g} norm(grad)={:g} norm(hess)={:g} alpha={:g} min_distance={:g} num_collisions={}",
            dx, grad_norm, hess.norm(), alpha, min_distance,
            m_normal_collisions->size());
    } while (dx > m_time_integrator->dt * 1e-3 && ++iter < 100);

    if (iter == 100) {
        logger().warn(
            "failed to converge: iter={} dx={:g} norm(grad)={:g}", iter, dx,
            grad_norm);
        // TODO: May want to exit here without pushing the new poses and
        // updating the time integrator, since the new poses may be unstable and
        // cause the simulation to blow up. However, this may cause the
        // simulation to get stuck if it fails to converge at the first step.
        // May want to add a flag to allow exiting early on convergence failure.
    } else {
        logger().info(
            "converged: iter={} dx={:g} norm(grad)={:g}", iter, dx, grad_norm);
    }

    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        poses[i].position = x.segment<3>(6 * i);
        poses[i].rotation = x.segment<3>(6 * i + 3);
    }

    m_pose_history.push_back(poses);

    Eigen::VectorXd X = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        X.segment<3>(12 * i) = x.segment<3>(6 * i);
        X.segment<9>(12 * i + 3) =
            rotation_vector_to_matrix(x.segment<3>(6 * i + 3)).reshaped();
    }

    m_time_integrator->update(X);

    m_t += m_time_integrator->dt;

    return iter < 100; // Return false if failed to converge, true otherwise
}

double Simulator::energy(Eigen::ConstRef<Eigen::VectorXd> x)
{
    return m_inertial_term->operator()(*m_bodies, x)
        + m_body_forces->operator()(*m_bodies, x)
        + m_ground_contact->operator()(*m_bodies, x)
        + m_barrier_potential->operator()(
            *m_normal_collisions, *m_bodies, m_bodies->vertices(x));
}

Eigen::VectorXd Simulator::gradient(Eigen::ConstRef<Eigen::VectorXd> x)
{
    const std::vector<Pose> poses = Pose::to_poses(x, m_bodies->dim());
    return m_inertial_term->gradient(*m_bodies, x)
        + m_body_forces->gradient(*m_bodies, x)
        + m_ground_contact->gradient(*m_bodies, x)
        + m_bodies->to_rigid_dof(
            poses,
            m_barrier_potential->gradient(
                *m_normal_collisions, *m_bodies, m_bodies->vertices(poses)));
}

Eigen::MatrixXd
Simulator::hessian(Eigen::ConstRef<Eigen::VectorXd> x, bool project_to_psd)
{
    const std::vector<Pose> poses = Pose::to_poses(x, m_bodies->dim());

    Eigen::VectorXd tmp = m_barrier_potential->gradient(
        *m_normal_collisions, *m_bodies, m_bodies->vertices(x));
    Eigen::SparseMatrix<double> tmp_hess = m_barrier_potential->hessian(
        *m_normal_collisions, *m_bodies, m_bodies->vertices(x),
        project_to_psd ? PSDProjectionMethod::CLAMP
                       : PSDProjectionMethod::NONE);

    return m_inertial_term->hessian(
               *m_bodies, x,
               project_to_psd ? PSDProjectionMethod::ABS
                              : PSDProjectionMethod::NONE)
        + m_ground_contact->hessian(*m_bodies, x, PSDProjectionMethod::CLAMP)
        + m_body_forces->hessian(
            *m_bodies, x,
            project_to_psd ? PSDProjectionMethod::CLAMP
                           : PSDProjectionMethod::NONE)
        + m_bodies->to_rigid_dof(poses, tmp, tmp_hess);
}

namespace {
    // TODO: Why does this not work with 0.5?
    inline constexpr double INFLATION_RADIUS_MULTIPLIER = 1.0;
} // namespace

void Simulator::update_candidates(
    Eigen::ConstRef<Eigen::VectorXd> x_t0,
    Eigen::ConstRef<Eigen::VectorXd> x_t1)
{
    m_candidates->build(
        *m_bodies, m_bodies->vertices(x_t0), m_bodies->vertices(x_t1),
        INFLATION_RADIUS_MULTIPLIER * m_barrier_potential->dhat());
}

void Simulator::update_collisions(
    Eigen::ConstRef<Eigen::VectorXd> x, bool update_candidates)
{
    Eigen::MatrixXd vertices = m_bodies->vertices(x);

    if (update_candidates) {
        m_candidates->build(
            *m_bodies, vertices,
            INFLATION_RADIUS_MULTIPLIER * m_barrier_potential->dhat());
    }

    m_normal_collisions->build(
        *m_candidates, *m_bodies, vertices, m_barrier_potential->dhat());
}

void Simulator::reset()
{
    m_t = 0.0;
    // Reset pose history to only contain the initial poses
    m_pose_history.resize(1);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x0.segment<3>(12 * i) = m_pose_history.front()[i].position;
        x0.segment<9>(12 * i + 3) =
            m_pose_history.front()[i].rotation_matrix().reshaped();
    }
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    Eigen::VectorXd a0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());

    m_time_integrator = std::make_unique<ImplicitEuler>(
        x0, v0, a0, m_time_integrator->dt, m_bodies->num_bodies());

    m_inertial_term = std::make_unique<InertialTerm>(m_time_integrator);
    m_candidates->clear();
    m_normal_collisions->clear();
}

} // namespace ipc::rigid