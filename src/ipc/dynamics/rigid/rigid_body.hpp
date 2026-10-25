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

    /// @brief The measure of the body that its mass is computed from.
    ///
    /// This is the enclosed volume for a closed 3D triangle mesh. Otherwise
    /// the mass is lumped, and this is the total lumped measure: the edge
    /// length for 2D and edge meshes, or one unit per vertex for point clouds
    /// and open surfaces.
    double volume() const { return m_volume; }

    /// @brief The total mass (density × volume()).
    double mass() const { return m_mass; }

    /// @brief The uniform density (mass() / volume()).
    double density() const { return m_mass / m_volume; }

    /// @brief The principal moments of inertia about the center of mass.
    ///
    /// In 3D these are the inertia tensor's eigenvalues, ordered to match the
    /// axes of principal_rotation(). In 2D this is a single value, the moment
    /// about the out-of-plane axis. When exactly two rotational DOFs are
    /// fixed in 3D, these are instead the diagonal of the world-frame inertia
    /// tensor.
    const VectorMax3d& moment_of_inertia() const { return m_moment_of_inertia; }

    /// @brief The diagonal second-moment matrix ∫ρ x̄x̄ᵀ in the principal frame (2×2 in 2D, 3×3 in 3D).
    ///
    /// This is the J in the affine mass matrix blkdiag(m I, J ⊗ I). When
    /// exactly two rotational DOFs are fixed in 3D, it is the diagonal of the
    /// world-frame second moment instead.
    const auto& second_moment() const { return m_second_moment; }

    /// @brief The rotation from the principal axes to the input frame.
    ///
    /// At construction we express the rest vertices in the principal frame
    /// and fold this rotation into the initial pose, so world positions are
    /// unchanged. We keep it to convert back to the input frame (e.g., when
    /// writing glTF). It is the identity when exactly two rotational DOFs are
    /// fixed in 3D, so the free rotation component stays on its world axis.
    const MatrixMax3d& principal_rotation() const
    {
        return m_principal_rotation;
    }

    /// @brief The constant external load applied to the body.
    ///
    /// The position holds the world-space force and the rotation holds the
    /// world-space torque (a scalar in 2D). The simulator folds it into the
    /// incremental potential each step, for DYNAMIC bodies only.
    const Pose& external_force() const { return m_external_force; }

    /// @brief The static bounding volume hierarchy of the body's elements.
    ///
    /// We build it once at construction from the rest vertices, in the body
    /// frame, so the rigid broad phase can query a moving body without
    /// rebuilding it.
    std::shared_ptr<const LBVH> bvh() const { return m_bvh; }

    /// @brief The largest distance from the center of mass to a rest vertex.
    ///
    /// The body stays inside a sphere of this radius about its position in
    /// every orientation, which is what the body-level broad phase uses.
    double bounding_radius() const { return m_bounding_radius; }

    /// @brief How the body is simulated (STATIC, KINEMATIC, or DYNAMIC).
    Type type() const { return m_type; }

    /// @brief Whether the body is DYNAMIC, i.e., responds to forces and
    /// contact.
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

    /// @brief Change the body's (uniform) density, rescaling its mass and inertia.
    ///
    /// Mass and inertia are both linear in density, and the principal axes do
    /// not depend on it, so we rescale them instead of recomputing the mass
    /// properties from the mesh.
    ///
    /// @note A running simulator picks up the new mass matrix at its next time step.
    /// @param density New density (must be positive).
    void set_density(const double density)
    {
        assert(density > 0);
        const double scale = density / this->density();
        m_mass *= scale;
        m_moment_of_inertia *= scale;
        m_second_moment.diagonal() *= scale;
    }

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

    /// @brief Diagonal second-moment matrix ∫ρ x̄x̄ᵀ in the principal frame
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, 3> m_second_moment;

    /// @brief Rotation matrix from the principal axes to the world frame
    /// This also stored in initial_pose.rotation upon construction.
    /// This is useful for converting to and from input world coordinates.
    /// @note Maybe this should be a rotation vector instead?
    /// @note Maybe we should store position as well?
    MatrixMax3d m_principal_rotation;

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