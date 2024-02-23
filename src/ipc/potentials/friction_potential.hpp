#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/friction/friction_collisions.hpp>

namespace ipc {

/// @brief The friction dissipative potential.
class FrictionPotential : public Potential<FrictionCollisions> {
    using Super = Potential;

public:
    /// @brief Construct a friction potential.
    /// @param epsv The smooth friction mollifier parameter \f$\epsilon_v\f$.
    explicit FrictionPotential(const double epsv);

    /// @brief Get the smooth friction mollifier parameter \f$\epsilon_v\f$.
    double epsv() const { return m_epsv; }

    /// @brief Set the smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @param epsv The smooth friction mollifier parameter \f$\epsilon_v\f$.
    void set_epsv(const double epsv)
    {
        assert(epsv > 0);
        m_epsv = epsv;
    }

    // -- Cumulative methods ---------------------------------------------------

    // NOTE: X in this context are vertex velocities.

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

    /// @brief Variable to differentiate the friction force with respect to.
    enum class DiffWRT {
        REST_POSITIONS,       ///< Differentiate w.r.t. rest positions
        LAGGED_DISPLACEMENTS, ///< Differentiate w.r.t. lagged displacements
        VELOCITIES            ///< Differentiate w.r.t. current velocities
    };

    /// @brief Compute the friction force from the given velocities.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param barrier_potential Barrier potential (used for normal force magnitude).
    /// @param barrier_stiffness Barrier stiffness (used for normal force magnitude).
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @param no_mu whether to not multiply by mu
    /// @return The friction force.
    Eigen::VectorXd force(
        const FrictionCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& lagged_displacements,
        const Eigen::MatrixXd& velocities,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const double dmin = 0,
        const bool no_mu = false) const;

    /// @brief Compute the Jacobian of the friction force wrt the velocities.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param barrier_potential Barrier potential (used for normal force magnitude).
    /// @param barrier_stiffness Barrier stiffness (used for normal force magnitude).
    /// @param wrt The variable to take the derivative with respect to.
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @return The Jacobian of the friction force wrt the velocities.
    Eigen::SparseMatrix<double> force_jacobian(
        const FrictionCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& lagged_displacements,
        const Eigen::MatrixXd& velocities,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const DiffWRT wrt,
        const double dmin = 0) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision
    /// @param velocities The collision stencil's velocities.
    /// @return The potential.
    double operator()(
        const FrictionCollision& collision,
        const VectorMax12d& velocities) const override;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision
    /// @param velocities The collision stencil's velocities.
    /// @return The gradient of the potential.
    VectorMax12d gradient(
        const FrictionCollision& collision,
        const VectorMax12d& velocities) const override;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision
    /// @param velocities The collision stencil's velocities.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const FrictionCollision& collision,
        const VectorMax12d& velocities,
        const bool project_hessian_to_psd = false) const override;

    /// @brief Compute the friction force.
    /// @param collision The collision
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param barrier_potential Barrier potential (used for normal force magnitude).
    /// @param barrier_stiffness Barrier stiffness (used for normal force magnitude).
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @param no_mu Whether to not multiply by mu
    /// @return Friction force
    VectorMax12d force(
        const FrictionCollision& collision,
        const VectorMax12d& rest_positions,
        const VectorMax12d& lagged_displacements,
        const VectorMax12d& velocities,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const double dmin = 0,
        const bool no_mu = false) const; //< whether to not multiply by mu

    /// @brief Compute the friction force Jacobian.
    /// @param collision The collision
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param barrier_potential Barrier potential (used for normal force magnitude).
    /// @param barrier_stiffness Barrier stiffness (used for normal force magnitude).
    /// @param wrt Variable to differentiate the friction force with respect to.
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @return Friction force Jacobian
    MatrixMax12d force_jacobian(
        const FrictionCollision& collision,
        const VectorMax12d& rest_positions,
        const VectorMax12d& lagged_displacements,
        const VectorMax12d& velocities,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const DiffWRT wrt,
        const double dmin = 0) const;

protected:
    /// @brief The smooth friction mollifier parameter \f$\epsilon_v\f$.
    double m_epsv;
};

} // namespace ipc