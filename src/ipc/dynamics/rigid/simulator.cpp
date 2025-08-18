#include "simulator.hpp"

#include <ipc/dynamics/rigid/inertial_term.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/time_integrator.hpp>

#include <igl/PI.h>

namespace ipc::rigid {

Simulator::Simulator(
    const std::shared_ptr<RigidBodies> _bodies,
    const std::vector<Pose>& initial_poses,
    const double dt)
    : m_bodies(std::move(_bodies))
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
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        // v0.segment<3>(12 * i).y() = 10;

        // ω = R₀ᵀω₀ (ω₀ expressed in body coordinates)
        Eigen::Vector3d omega(-100 * igl::PI / 180, 0, 0);
        omega = (*m_bodies)[i].R0().transpose() * omega;
        Eigen::Matrix3d Q_t0 = initial_poses[i].rotation_matrix();
        v0.segment<9>(12 * i + 3) =
            (Q_t0 * cross_product_matrix(omega)).reshaped();
    }
    Eigen::VectorXd a0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());

    // Initialize the time integrator
    m_time_integrator =
        std::make_shared<ImplicitEuler>(x0, v0, a0, dt, m_bodies->num_bodies());

    // Initialize the inertial term
    m_inertial_term = std::make_shared<InertialTerm>(m_time_integrator);
}

void Simulator::run(
    // const double dt,
    const double t_end,
    const std::function<void(void)> callback)
{
    if (t_end <= m_t) {
        logger().warn(
            "simulation already complete: t={:g} t_end={:g}", m_t, t_end);
        return;
    }

    // m_time_integrator->dt = dt;

    while (m_t < t_end) {
        step(); // increments time by dt
        callback();
    }
}

void Simulator::step()
{
    std::vector<Pose> poses = m_pose_history.back();

    m_inertial_term->update(*m_bodies);

    Eigen::VectorXd x = Eigen::VectorXd(6 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x.segment<3>(6 * i) = m_pose_history.back()[i].position;
        x.segment<3>(6 * i + 3) = m_pose_history.back()[i].rotation;
    }

    double dx, grad_norm;
    int iter = 0;
    do {
        Eigen::VectorXd grad = m_inertial_term->gradient(*m_bodies, x);
        if ((grad_norm = grad.norm()) < 1e-6) {
            break;
        }
        Eigen::MatrixXd hess =
            m_inertial_term->hessian(*m_bodies, x, PSDProjectionMethod::ABS);
        Eigen::VectorXd step = -hess.lu().solve(grad);
        dx = step.norm();
        double alpha = 1.0;
        double Ex = (*m_inertial_term)(*m_bodies, x);
        while ((*m_inertial_term)(*m_bodies, x + alpha * step) >= Ex) {
            alpha *= 0.5;
        }
        x += alpha * step;
        logger().debug(
            "step: dx={:g} norm(grad)={:g} norm(hess)={:g} alpha={:g}", dx,
            grad_norm, hess.norm(), alpha);
    } while (dx > m_time_integrator->dt * 1e-3 && ++iter < 100);
    logger().info(
        "converged: iter={} dx={:g} norm(grad)={:g}", iter, dx, grad_norm);

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

    m_time_integrator = std::make_shared<ImplicitEuler>(
        x0, v0, a0, m_time_integrator->dt, m_bodies->num_bodies());

    m_inertial_term = std::make_shared<InertialTerm>(m_time_integrator);
}

} // namespace ipc::rigid