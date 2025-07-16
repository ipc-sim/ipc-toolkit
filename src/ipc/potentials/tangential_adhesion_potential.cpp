#include "tangential_adhesion_potential.hpp"

#include <ipc/adhesion/adhesion.hpp>

namespace ipc {

TangentialAdhesionPotential::TangentialAdhesionPotential(const double eps_a)
    : Super()
{
    set_eps_a(eps_a);
}

double TangentialAdhesionPotential::mu_f0(
    const double x, const double mu_s, const double mu_k) const
{
    return smooth_mu_a0(x, mu_s, mu_k, eps_a());
}

double TangentialAdhesionPotential::mu_f1_over_x(
    const double x, const double mu_s, const double mu_k) const
{
    return smooth_mu_a1_over_x(x, mu_s, mu_k, eps_a());
}

double TangentialAdhesionPotential::mu_f2_x_minus_mu_f1_over_x3(
    const double x, const double mu_s, const double mu_k) const
{
    return smooth_mu_a2_x_minus_mu_a1_over_x3(x, mu_s, mu_k, eps_a());
}

} // namespace ipc