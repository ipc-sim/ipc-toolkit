#pragma once

#include <ipc/potentials/tangential_potential.hpp>

namespace ipc {

/// @brief The friction dissipative potential.
class FrictionPotential : public TangentialPotential {
    using Super = TangentialPotential;

public:
    /// @brief Construct a friction potential.
    /// @param eps_v The smooth friction mollifier parameter \f$\epsilon_v\f$.
    explicit FrictionPotential(const double eps_v);

    /// @brief Get the smooth friction mollifier parameter \f$\epsilon_v\f$.
    double eps_v() const { return m_eps_v; }

    /// @brief Set the smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @param eps_v The smooth friction mollifier parameter \f$\epsilon_v\f$.
    void set_eps_v(const double eps_v)
    {
        assert(eps_v > 0);
        m_eps_v = eps_v;
    }

protected:
    /// @brief Compute the value of the ∫ μ(y) f₁(y) dy, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static friction.
    /// @param mu_k Coefficient of kinetic friction.
    /// @return The value of the integral at x.
    double
    mu_f0(const double x, const double mu_s, const double mu_k) const override;

    /// @brief Compute the value of the [μ(y) f₁(y)] / x, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static friction.
    /// @param mu_k Coefficient of kinetic friction.
    /// @return The value of the product at x.
    double mu_f1_over_x(
        const double x, const double mu_s, const double mu_k) const override;

    /// @brief Compute the value of [(d/dx (μ(y) f₁(y))) x - μ(y) f₁(y)] / x³, where f₁ is the first derivative of the smooth mollifier.
    /// @param x The tangential relative speed.
    /// @param mu_s Coefficient of static friction.
    /// @param mu_k Coefficient of kinetic friction.
    /// @return The value of the derivative at x.
    double mu_f2_x_minus_mu_f1_over_x3(
        const double x, const double mu_s, const double mu_k) const override;

    /// @brief Check if the potential is dynamic.
    /// @param speed The tangential relative speed.
    /// @return True if the potential is dynamic, false otherwise.
    bool is_dynamic(const double speed) const override
    {
        return speed > eps_v();
    }

    /// @brief The smooth friction mollifier parameter \f$\epsilon_v\f$.
    double m_eps_v;
};

} // namespace ipc