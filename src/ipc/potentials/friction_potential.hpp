#pragma once

#include <ipc/potentials/tangential_potential.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>

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
    
    /// @brief Set the material friction coefficient lookup table
    /// @param table Table of friction coefficients indexed by material IDs
    void set_material_friction_table(const std::shared_ptr<Eigen::MatrixXd>& table)
    {
        material_friction_table = table;
    }
    
    /// @brief Get the material friction coefficient lookup table
    /// @return Shared pointer to the friction table or nullptr if not set
    std::shared_ptr<Eigen::MatrixXd> get_material_friction_table() const
    {
        return material_friction_table;
    }

protected:
    double f0(const double x) const override;
    double f1_over_x(const double x) const override;
    double f2_x_minus_f1_over_x3(const double x) const override;
    
    // Additional functions for handling static/kinetic friction
    double f0_mus(const double x, const double mu_s, const double mu_k) const;
    double f1_over_x_mus(const double x, const double mu_s, const double mu_k) const; 
    double f2_x_minus_f1_over_x3_mus(const double x, const double mu_s, const double mu_k) const;

    bool is_dynamic(const double speed) const override
    {
        return speed > eps_v();
    }

    /// @brief The smooth friction mollifier parameter \f$\epsilon_v\f$.
    double m_eps_v;
    
    /// @brief Table of friction coefficients indexed by material IDs
    std::shared_ptr<Eigen::MatrixXd> material_friction_table;
};

} // namespace ipc