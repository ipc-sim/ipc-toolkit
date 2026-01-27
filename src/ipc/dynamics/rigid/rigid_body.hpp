#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
class BVH;
}

namespace ipc::rigid {

class RigidBody {
public:
    enum class Type : uint8_t {
        /// @brief Static rigid body, does not move
        STATIC,
        /// @brief Kinematic rigid body, moves but does not respond to forces
        KINEMATIC,
        /// @brief Dynamic rigid body, moves and responds to forces
        DYNAMIC
    };

public:
    RigidBody(
        Eigen::Ref<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double density,
        Pose& initial_pose);

    double mass() const { return m_mass; }
    const VectorMax3d& moment_of_inertia() const { return m_moment_of_inertia; }
    // NOLINTNEXTLINE(readability-identifier-naming)
    const Eigen::DiagonalMatrix<double, 3>& J() const { return m_J; }
    // NOLINTNEXTLINE(readability-identifier-naming)
    const MatrixMax3d& R0() const { return m_R0; }
    const Pose& external_force() const { return m_external_force; }
    std::shared_ptr<const BVH> bvh() const { return m_bvh; }
    double bounding_radius() const { return m_bounding_radius; }

private:
    /// @brief Total mass of the rigid body
    double m_mass;

    /// @brief Moment of inertia measured with respect to the principal axes
    VectorMax3d m_moment_of_inertia;

    // NOLINTNEXTLINE(readability-identifier-naming)
    Eigen::DiagonalMatrix<double, 3> m_J;

    /// @brief Rotation matrix from the principal axes to the world frame
    /// This also stored in initial_pose.rotation upon construction.
    /// This is useful for converting to and from input world coordinates.
    /// @note Maybe this should be a rotation vector instead?
    /// @note Maybe we should store position as well?
    MatrixMax3d m_R0; // NOLINT(readability-identifier-naming)

    /// @brief External force and torque applied to the rigid body
    Pose m_external_force;

    /// @brief Statically constructed bounding volume hierarchy for collision detection
    /// @note This is defined in the inertial reference frame
    std::shared_ptr<BVH> m_bvh;

    /// @brief Bounding radius of the rigid body
    double m_bounding_radius;
};

} // namespace ipc::rigid