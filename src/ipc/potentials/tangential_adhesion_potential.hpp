#pragma once

#include <ipc/potentials/tangential_potential.hpp>

namespace ipc {

/// @brief The tangential adhesion potential.
class TangentialAdhesionPotential : public TangentialPotential {
    using Super = TangentialPotential;

public:
    /// @brief Construct a tangential adhesion potential.
    /// @param eps_a The tangential adhesion mollifier parameter \f$\epsilon_a\f$.
    explicit TangentialAdhesionPotential(const double eps_a);

    /// @brief Get the tangential adhesion mollifier parameter \f$\epsilon_a\f$.
    double eps_a() const { return m_eps_a; }

    /// @brief Set the tangential adhesion mollifier parameter \f$\epsilon_v\f$.
    /// @param eps_a The tangential adhesion mollifier parameter \f$\epsilon_v\f$.
    void set_eps_a(const double eps_a)
    {
        assert(eps_a > 0);
        m_eps_a = eps_a;
    }

protected:
    /// @brief Compute the value of the ∫ μ(y) f₁(y) dy, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static adhesion.
    /// @param mu_k Coefficient of kinetic adhesion.
    /// @return The value of the integral at x.
    double
    mu_f0(const double x, const double mu_s, const double mu_k) const override;

    /// @brief Compute the value of the [μ(y) f₁(y)] / x, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static adhesion.
    /// @param mu_k Coefficient of kinetic adhesion.
    /// @return The value of the product at x.
    double mu_f1_over_x(
        const double x, const double mu_s, const double mu_k) const override;

    /// @brief Compute the value of [(d/dx (μ(y) f₁(y))) x - μ(y) f₁(y)] / x³, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static adhesion.
    /// @param mu_k Coefficient of kinetic adhesion.
    /// @return The value of the derivative at x.
    double mu_f2_x_minus_mu_f1_over_x3(
        const double x, const double mu_s, const double mu_k) const override;

    /// @brief Check if the potential is dynamic.
    /// @param speed The tangential relative speed.
    /// @return True if the potential is dynamic, false otherwise.
    bool is_dynamic(const double speed) const override
    {
        return speed > eps_a();
    }

    /// @brief The tangential adhesion mollifier parameter \f$\epsilon_a\f$.
    double m_eps_a;
};

} // namespace ipc