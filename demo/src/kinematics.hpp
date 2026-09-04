#pragma once

#include <ipc/dynamics/affine/pose.hpp>
#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/dynamics/to_affine.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <memory>
#include <vector>

namespace ipc {
class BroadPhase;
class Candidates;
} // namespace ipc

namespace ipc::rigid {
class RigidBodies;
} // namespace ipc::rigid

namespace ipc::demo {

/// @brief The kinematics of a discretization: how the simulation state x
/// maps to the world.
///
/// Encapsulates the DOF layout, the world-vertex map and its chain rule, the
/// bridge to the time-integrator state, and the discretization-appropriate
/// collision candidates and CCD — while all forces/energies live elsewhere.
/// Today this hides the differences between the rigid (6-DOF [p; θ] per
/// body) and affine (12-DOF [p; vec(A) column-major] per body)
/// parameterizations; future discretizations (e.g., FEM nodal DOFs or SPH
/// particles, where x is the vertex/particle positions and the world map is
/// the identity) slot in as additional implementations.
class Kinematics {
public:
    virtual ~Kinematics();

    /// @brief Create the rigid (6-DOF [p; θ] per body) DOF model.
    static std::unique_ptr<Kinematics>
    create_rigid(const std::shared_ptr<const rigid::RigidBodies>& bodies);

    /// @brief Create the affine (12-DOF [p; vec(A)] per body) DOF model.
    static std::unique_ptr<Kinematics>
    create_affine(const std::shared_ptr<const rigid::RigidBodies>& bodies);

    /// @brief Total number of optimization DOFs.
    virtual int ndof() const = 0;

    // -----------------------------------------------------------------------
    // World map and chain rule
    // -----------------------------------------------------------------------

    /// @brief World-space collision-mesh vertices at x.
    virtual Eigen::MatrixXd
    world_vertices(Eigen::ConstRef<Eigen::VectorXd> x) const = 0;

    /// @brief Map a vertex-space gradient to the DOFs: Jᵀ ∇V.
    virtual Eigen::VectorXd map_vertex_gradient(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> grad) const = 0;

    /// @brief Map a vertex-space Hessian to the DOFs: Jᵀ H J (+ the
    /// second-order term Σ ∇V·d²V/dx² for the rigid parameterization; the
    /// affine map is linear so its second-order term is zero).
    virtual Eigen::SparseMatrix<double> map_vertex_hessian(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> grad,
        const Eigen::SparseMatrix<double>& hess) const = 0;

    // -----------------------------------------------------------------------
    // Pose conversions
    // -----------------------------------------------------------------------

    /// @brief The to-affine map (DOFs → affine pose coordinates) and its chain
    /// rule for this parameterization.
    std::shared_ptr<const dynamics::ToAffine> to_affine() const
    {
        return m_to_affine;
    }

    /// @brief Per-body affine poses at x (exact for rigid: A = R(θ)).
    virtual std::vector<affine::Pose>
    poses(Eigen::ConstRef<Eigen::VectorXd> x) const = 0;

    /// @brief Build the DOF vector from rigid poses.
    virtual Eigen::VectorXd
    dof(const std::vector<rigid::Pose>& poses) const = 0;

    // -----------------------------------------------------------------------
    // Time-integrator bridge
    // -----------------------------------------------------------------------
    // The integrator state is affine-shaped ([p; vec(Q) column-major] per
    // body) for both models; only the rigid model needs a conversion.

    /// @brief Convert the optimization DOFs to the integrator state.
    virtual Eigen::VectorXd
    to_integrator_state(Eigen::ConstRef<Eigen::VectorXd> x) const = 0;

    /// @brief Convert the integrator state to the optimization DOFs.
    /// @note For the rigid model the log map yields the canonical rotation
    /// vector (‖θ‖ ≤ π), which also keeps the optimization variable bounded.
    virtual Eigen::VectorXd
    from_integrator_state(Eigen::ConstRef<Eigen::VectorXd> X) const = 0;

    // -----------------------------------------------------------------------
    // Collision candidates and CCD
    // -----------------------------------------------------------------------

    /// @brief Build the discrete collision candidates at x.
    virtual void update_candidates(
        Eigen::ConstRef<Eigen::VectorXd> x,
        const double inflation_radius,
        BroadPhase* broad_phase) = 0;

    /// @brief Build the continuous collision candidates for the step x0 → x1.
    virtual void update_candidates(
        Eigen::ConstRef<Eigen::VectorXd> x0,
        Eigen::ConstRef<Eigen::VectorXd> x1,
        const double inflation_radius,
        BroadPhase* broad_phase) = 0;

    /// @brief Compute a collision-free maximum step size from x0 to x1 ∈ [0, 1].
    /// Also updates the collision candidates for the swept interval.
    virtual double max_step_size(
        Eigen::ConstRef<Eigen::VectorXd> x0,
        Eigen::ConstRef<Eigen::VectorXd> x1,
        const double inflation_radius,
        BroadPhase* broad_phase,
        const double min_distance = 0.0) = 0;

    /// @brief The current collision candidates.
    virtual const Candidates& candidates() const = 0;

    /// @brief Clear the collision candidates.
    virtual void clear_candidates() = 0;

protected:
    explicit Kinematics(
        const std::shared_ptr<const rigid::RigidBodies>& bodies);

    /// @brief The bodies in the simulation.
    std::shared_ptr<const rigid::RigidBodies> m_bodies;

    /// @brief The to-affine map for this parameterization (set by subclasses).
    std::shared_ptr<const dynamics::ToAffine> m_to_affine;
};

} // namespace ipc::demo
