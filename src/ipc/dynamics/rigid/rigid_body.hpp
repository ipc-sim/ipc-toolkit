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
        const double density,
        Pose& initial_pose);

    double mass() const { return m_mass; }
    const VectorMax3d& moment_of_inertia() const { return m_moment_of_inertia; }
    const Eigen::DiagonalMatrix<double, 3>& J() const { return m_J; }
    const MatrixMax3d& R0() const { return m_R0; }
    const Pose& external_force() const { return m_external_force; }

private:
    /// @brief Total mass of the rigid body
    double m_mass;

    /// @brief Moment of inertia measured with respect to the principal axes
    VectorMax3d m_moment_of_inertia;

    Eigen::DiagonalMatrix<double, 3> m_J;

    /// @brief Rotation matrix from the principal axes to the world frame
    /// This also stored in initial_pose.rotation upon construction.
    /// This is useful for converting to and from input world coordinates.
    /// @note Maybe this should be a rotation vector instead?
    /// @note Maybe we should store position as well?
    MatrixMax3d m_R0;

    /// @brief External force and torque applied to the rigid body
    Pose m_external_force;
};

} // namespace ipc::rigid