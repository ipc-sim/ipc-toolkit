#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/friction/friction_constraints.hpp>

namespace ipc {

class FrictionPotential : public Potential<FrictionConstraints> {
    using Super = Potential;
    using Contacts = FrictionConstraints;
    using Contact = Super::Contact;

public:
    /// @brief Construct a friction potential.
    /// @param epsv The smooth friction mollifier parameter \f$\epsilon_v\f$.
    FrictionPotential(const double epsv);

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

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
    /// @brief Compute the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The potential.
    double
    operator()(const Contact& contact, const VectorMax12d& x) const override;

    /// @brief Compute the gradient of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The gradient of the potential.
    VectorMax12d
    gradient(const Contact& contact, const VectorMax12d& x) const override;

    /// @brief Compute the hessian of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const Contact& contact,
        const VectorMax12d& x,
        const bool project_hessian_to_psd = true) const override;

    /// @brief The smooth friction mollifier parameter \f$\epsilon_v\f$.
    double m_epsv;
};

} // namespace ipc