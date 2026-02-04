#include "simulator.hpp"

#include <ipc/dynamics/rigid/body_forces.hpp>
#include <ipc/dynamics/rigid/ground_contact.hpp>
#include <ipc/dynamics/rigid/inertial_term.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/time_integrator.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <igl/PI.h>

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
    m_inertial_term = std::make_shared<InertialTerm>(m_time_integrator);

    // Iniaizlie the body force term
    m_body_forces = std::make_shared<BodyForces>(m_time_integrator);
    m_body_forces->set_gravity(Eigen::Vector3d(0, -9.81, 0));

    // Initialize the ground contact handler
    m_ground_contact = std::make_shared<GroundContact>(0.0);
    m_ground_contact->set_dhat(0.1);
}

void Simulator::run(
    // const double dt,
    const double t_end,
    const std::function<void(void)>& callback)
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
    m_body_forces->update(*m_bodies);

    Eigen::VectorXd x = Eigen::VectorXd(6 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x.segment<3>(6 * i) = m_pose_history.back()[i].position;
        x.segment<3>(6 * i + 3) = m_pose_history.back()[i].rotation;
    }

    double dx, grad_norm;
    int iter = 0;
    do {
        Eigen::VectorXd grad = gradient(x);
        if ((grad_norm = grad.norm()) < 1e-6) {
            break;
        }
        Eigen::MatrixXd hess = hessian(x);
        Eigen::VectorXd step = -hess.llt().solve(grad);
        dx = step.norm();
        double alpha = 1.0;
        double Ex = energy(x);
        while (energy(x + alpha * step) >= Ex) {
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

double Simulator::energy(Eigen::ConstRef<Eigen::VectorXd> x)
{
    return m_inertial_term->operator()(*m_bodies, x)
        + m_body_forces->operator()(*m_bodies, x)
        + m_ground_contact->operator()(*m_bodies, x);
}

Eigen::VectorXd Simulator::gradient(Eigen::ConstRef<Eigen::VectorXd> x)
{
    return m_inertial_term->gradient(*m_bodies, x)
        + m_body_forces->gradient(*m_bodies, x)
        + m_ground_contact->gradient(*m_bodies, x);
}

Eigen::MatrixXd Simulator::hessian(Eigen::ConstRef<Eigen::VectorXd> x)
{
    return m_inertial_term->hessian(*m_bodies, x, PSDProjectionMethod::ABS)
        + m_body_forces->hessian(*m_bodies, x, PSDProjectionMethod::CLAMP)
        + m_ground_contact->hessian(*m_bodies, x, PSDProjectionMethod::CLAMP);
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