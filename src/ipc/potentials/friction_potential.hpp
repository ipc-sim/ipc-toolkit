#pragma once

#include <ipc/potentials/tangential_potential.hpp>

namespace ipc {

/// @brief The friction dissipative potential.
class FrictionPotential : public TangentialPotential {
    using Super = TangentialPotential;

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

protected:
    double f0(const double x) const override;
    double f1_over_x(const double x) const override;
    double df1_x_minus_f1_over_x3(const double x) const override;

    bool is_dynamic(const double speed) const override
    {
        return speed > epsv();
    }

    /// @brief The smooth friction mollifier parameter \f$\epsilon_v\f$.
    double m_epsv;
};

} // namespace ipc