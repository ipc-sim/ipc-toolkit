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
    double f0(const double x) const override;
    double f1_over_x(const double x) const override;
    double f2_x_minus_f1_over_x3(const double x) const override;

    bool is_dynamic(const double speed) const override
    {
        return speed > eps_a();
    }

    /// @brief The tangential adhesion mollifier parameter \f$\epsilon_a\f$.
    double m_eps_a;
};

} // namespace ipc