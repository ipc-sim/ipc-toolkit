#include "simulator.hpp"

#include <ipc/dynamics/rigid/inertial_term.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/time_integrator.hpp>

namespace ipc::rigid {

Simulator::Simulator(
    const std::shared_ptr<RigidBodies> _bodies,
    const std::vector<Pose>& initial_poses,
    const double dt)
    : m_bodies(std::move(_bodies))
{
    m_pose_history.push_back(initial_poses);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x0.segment<3>(12 * i) = initial_poses[i].position;
        x0.segment<9>(12 * i) = initial_poses[i].rotation_matrix().reshaped();
    }
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
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

    // TODO: Solve the implicit time integration problem here
    // poses = m_solver->solve(poses);
    Pose delta(Eigen::Vector3d::Zero(), Eigen::Vector3d(0.05, 0, 0));
    for (int i = 0; i < m_bodies->num_bodies(); ++i) {
        poses[i] = delta * poses[i];
    }

    m_pose_history.push_back(poses);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x.segment<3>(12 * i) = poses[i].position;
        x.segment<9>(12 * i) = poses[i].rotation_matrix().reshaped();
    }

    m_time_integrator->update(x);

    m_t += m_time_integrator->dt;
}

void Simulator::reset()
{
    m_t = 0.0;
    // Reset pose history to only contain the initial poses
    m_pose_history.resize(1);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    for (size_t i = 0; i < m_bodies->num_bodies(); ++i) {
        x0.segment<3>(12 * i) = m_pose_history[0][i].position;
        x0.segment<9>(12 * i) =
            m_pose_history[0][i].rotation_matrix().reshaped();
    }
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());
    Eigen::VectorXd a0 = Eigen::VectorXd::Zero(12 * m_bodies->num_bodies());

    m_time_integrator = std::make_shared<ImplicitEuler>(
        x0, v0, a0, m_time_integrator->dt, m_bodies->num_bodies());

    m_inertial_term = std::make_shared<InertialTerm>(m_time_integrator);
}

} // namespace ipc::rigid