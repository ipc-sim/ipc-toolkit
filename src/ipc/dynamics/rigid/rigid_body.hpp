#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
class LBVH;
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
    /// @brief Construct a rigid body from a mesh and density.
    /// @param vertices Vertices of the mesh (will be modified to be centered at the center of mass)
    /// @param edges Edges of the mesh (2D)
    /// @param faces Faces of the mesh (3D) or empty (2D)
    /// @param density Density of the rigid body
    /// @param initial_pose Initial pose of the rigid body (will be modified to account for the center of mass)
    /// @param is_dof_fixed Per-DOF fixed flags ([position | rotation]; empty means all free). Positional flags refer to world axes; rotational flags to the components of the rotation vector.
    RigidBody(
        Eigen::Ref<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double density,
        Pose& initial_pose,
        const VectorMax6b& is_dof_fixed = VectorMax6b());

    // ---- Getters ------------------------------------------------------------
    double volume() const { return m_volume; }
    double mass() const { return m_mass; }
    double density() const { return m_mass / m_volume; }
    const VectorMax3d& moment_of_inertia() const { return m_moment_of_inertia; }
    /// @brief The diagonal second-moment matrix ∫ρ x̄x̄ᵀ in the principal frame (2×2 in 2D, 3×3 in 3D).
    const auto& J() const { return m_J; }          // NOLINT
    const MatrixMax3d& R0() const { return m_R0; } // NOLINT
    const Pose& external_force() const { return m_external_force; }
    std::shared_ptr<const LBVH> bvh() const { return m_bvh; }
    double bounding_radius() const { return m_bounding_radius; }

    Type type() const { return m_type; }
    bool is_dynamic() const { return m_type == Type::DYNAMIC; }

    /// @brief Per-DOF fixed flags ([position | rotation] of size ndof).
    const VectorMax6b& is_dof_fixed() const { return m_is_dof_fixed; }

    // ---- Setters ------------------------------------------------------------

    void set_external_force(const Pose& external_force)
    {
        assert(
            external_force.position.size() == m_external_force.position.size());
        assert(
            external_force.rotation.size() == m_external_force.rotation.size());
        m_external_force = external_force;
    }

    void set_type(const Type type) { m_type = type; }

    /// @brief Restore the fixed-DOF flags (used by Simulator::reset()).
    /// @warning Must match the construction-time flags: fixing rotational
    /// DOFs changes the inertia handling, which is baked at construction.
    void set_is_dof_fixed(const VectorMax6b& is_dof_fixed)
    {
        assert(is_dof_fixed.size() == m_is_dof_fixed.size());
        m_is_dof_fixed = is_dof_fixed;
    }

    /// @brief Convert this body to a STATIC body (all DOF fixed).
    /// @note Velocities need no explicit zeroing here: they live in the time
    /// integrator's history and decay to zero once the DOF are pinned.
    void convert_to_static()
    {
        m_type = Type::STATIC;
        m_is_dof_fixed.setOnes();
    }

private:
    /// @brief Volume of the rigid body
    double m_volume;

    /// @brief Total mass of the rigid body
    double m_mass;

    /// @brief Moment of inertia measured with respect to the principal axes
    VectorMax3d m_moment_of_inertia;

    // NOLINTNEXTLINE(readability-identifier-naming)
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, 3> m_J;

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
    std::shared_ptr<LBVH> m_bvh;

    /// @brief Bounding radius of the rigid body
    double m_bounding_radius;

    /// @brief How this body is simulated (STATIC/KINEMATIC/DYNAMIC).
    Type m_type = Type::DYNAMIC;

    /// @brief Per-DOF fixed flags ([position | rotation] of size ndof).
    VectorMax6b m_is_dof_fixed;
};

} // namespace ipc::rigid