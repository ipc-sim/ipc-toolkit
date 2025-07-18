#include "friction_potential.hpp"

#include <ipc/friction/smooth_mu.hpp>

namespace ipc {

FrictionPotential::FrictionPotential(const double eps_v) : Super()
{
    set_eps_v(eps_v);
}

double FrictionPotential::mu_f0(
    const double x, const double mu_s, const double mu_k) const
{
    return smooth_mu_f0(x, mu_s, mu_k, eps_v());
}

double FrictionPotential::mu_f1_over_x(
    const double x, const double mu_s, const double mu_k) const
{
    return smooth_mu_f1_over_x(x, mu_s, mu_k, eps_v());
}

double FrictionPotential::mu_f2_x_minus_mu_f1_over_x3(
    const double x, const double mu_s, const double mu_k) const
{
    return smooth_mu_f2_x_minus_mu_f1_over_x3(x, mu_s, mu_k, eps_v());
}

} // namespace ipc