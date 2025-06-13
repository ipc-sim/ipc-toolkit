#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc::rigid {

class RigidBody {
public:
    RigidBody(
        Eigen::Ref<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        Pose<>& initial_pose);

    double mass() const { return m_mass; }
    const VectorMax3d& moment_of_inertia() const { return m_moment_of_inertia; }
    const Eigen::Matrix3d& J() const { return m_J; }
    const Pose<>& external_force() const { return m_external_force; }

private:
    /// @brief Total mass of the rigid body
    double m_mass;

    /// @brief Moment of inertia measured with respect to the principal axes
    VectorMax3d m_moment_of_inertia;

    Eigen::Matrix3d m_J;

    /// @brief External force and torque applied to the rigid body
    Pose<> m_external_force;
};

} // namespace ipc::rigid