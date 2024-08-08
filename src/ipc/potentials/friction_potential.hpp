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
    double f0(const double x) const override;
    double f1_over_x(const double x) const override;
    double f2_x_minus_f1_over_x3(const double x) const override;

    bool is_dynamic(const double speed) const override
    {
        return speed > eps_v();
    }

    /// @brief The smooth friction mollifier parameter \f$\epsilon_v\f$.
    double m_eps_v;
};

} // namespace ipc