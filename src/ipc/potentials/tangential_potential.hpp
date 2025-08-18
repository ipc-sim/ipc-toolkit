#pragma once

#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/potential.hpp>

namespace ipc {

/// @brief A tangential dissipative potential.
class TangentialPotential : public Potential<TangentialCollisions> {
    using Super = Potential;

public:
    virtual ~TangentialPotential() { }

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
    /// @param normal_potential Normal potential (used for normal force magnitude).
    /// @param normal_stiffness Normal stiffness (used for normal force magnitude).
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @param no_mu whether to not multiply by mu
    /// @return The friction force.
    Eigen::VectorXd force(
        const TangentialCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
        Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
        Eigen::ConstRef<Eigen::MatrixXd> velocities,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
        const double dmin = 0,
        const bool no_mu = false) const;

    /// @brief Compute the Jacobian of the friction force wrt the velocities.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param normal_potential Normal potential (used for normal force magnitude).
    /// @param normal_stiffness Normal stiffness (used for normal force magnitude).
    /// @param wrt The variable to take the derivative with respect to.
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @return The Jacobian of the friction force wrt the velocities.
    Eigen::SparseMatrix<double> force_jacobian(
        const TangentialCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
        Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
        Eigen::ConstRef<Eigen::MatrixXd> velocities,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
        const DiffWRT wrt,
        const double dmin = 0) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision
    /// @param velocities The collision stencil's velocities.
    /// @return The potential.
    double operator()(
        const TangentialCollision& collision,
        Eigen::ConstRef<VectorMax12d> velocities) const override;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision
    /// @param velocities The collision stencil's velocities.
    /// @return The gradient of the potential.
    VectorMax12d gradient(
        const TangentialCollision& collision,
        Eigen::ConstRef<VectorMax12d> velocities) const override;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision
    /// @param velocities The collision stencil's velocities.
    /// @param project_hessian_to_psd Whether to project the hessian to the positive semi-definite cone.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const TangentialCollision& collision,
        Eigen::ConstRef<VectorMax12d> velocities,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override;

    /// @brief Compute the friction force.
    /// @param collision The collision
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param normal_potential Normal potential (used for normal force magnitude).
    /// @param normal_stiffness Normal stiffness (used for normal force magnitude).
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @param no_mu Whether to not multiply by mu
    /// @return Friction force
    VectorMax12d force(
        const TangentialCollision& collision,
        Eigen::ConstRef<VectorMax12d> rest_positions,
        Eigen::ConstRef<VectorMax12d> lagged_displacements,
        Eigen::ConstRef<VectorMax12d> velocities,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
        const double dmin = 0,
        const bool no_mu = false) const; //< whether to not multiply by mu

    /// @brief Compute the friction force Jacobian.
    /// @param collision The collision
    /// @param rest_positions Rest positions of the vertices (rowwise).
    /// @param lagged_displacements Previous displacements of the vertices (rowwise).
    /// @param velocities Current displacements of the vertices (rowwise).
    /// @param normal_potential Normal potential (used for normal force magnitude).
    /// @param normal_stiffness Noraml stiffness (used for normal force magnitude).
    /// @param wrt Variable to differentiate the friction force with respect to.
    /// @param dmin Minimum distance (used for normal force magnitude).
    /// @return Friction force Jacobian
    MatrixMax12d force_jacobian(
        const TangentialCollision& collision,
        Eigen::ConstRef<VectorMax12d> rest_positions,
        Eigen::ConstRef<VectorMax12d> lagged_displacements,
        Eigen::ConstRef<VectorMax12d> velocities,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
        const DiffWRT wrt,
        const double dmin = 0) const;

protected:
    /// @brief Compute the value of the ∫ μ(y) f₁(y) dy, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static friction.
    /// @param mu_k Coefficient of kinetic friction.
    /// @return The value of the integral at x.
    virtual double
    mu_f0(const double x, const double mu_s, const double mu_k) const = 0;

    /// @brief Compute the value of the [μ(y) f₁(y)] / x, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static friction.
    /// @param mu_k Coefficient of kinetic friction.
    /// @return The value of the product at x.
    virtual double mu_f1_over_x(
        const double x, const double mu_s, const double mu_k) const = 0;

    /// @brief Compute the value of [(d/dx (μ(y) f₁(y))) x - μ(y) f₁(y)] / x³, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static friction.
    /// @param mu_k Coefficient of kinetic friction.
    /// @return The value of the derivative at x.
    virtual double mu_f2_x_minus_mu_f1_over_x3(
        const double x, const double mu_s, const double mu_k) const = 0;

    /// @brief Check if the speed is dynamic.
    /// @param speed The tangential relative speed.
    /// @return True if the speed is dynamic, false otherwise.
    virtual bool is_dynamic(const double speed) const = 0;
};

} // namespace ipc