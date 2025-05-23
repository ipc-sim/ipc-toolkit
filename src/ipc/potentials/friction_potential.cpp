#include "friction_potential.hpp"

namespace ipc {

FrictionPotential::FrictionPotential(const double eps_v) : Super()
{
    set_eps_v(eps_v);
}

double FrictionPotential::f0(const double x) const
{
    return smooth_friction_f0(x, eps_v());
}

double FrictionPotential::f1_over_x(const double x) const
{
    return smooth_friction_f1_over_x(x, eps_v());
}

double FrictionPotential::f2_x_minus_f1_over_x3(const double x) const
{
    return smooth_friction_f2_x_minus_f1_over_x3(x, eps_v());
}

// Functions for static/kinetic friction coefficients
double FrictionPotential::f0_mus(const double x, const double mu_s, const double mu_k) const
{
    return smooth_friction_f0_mus(x, eps_v(), mu_s, mu_k);
}

double FrictionPotential::f1_over_x_mus(const double x, const double mu_s, const double mu_k) const
{
    return smooth_friction_f1_over_x_mus(x, eps_v(), mu_s, mu_k);
}

double FrictionPotential::f2_x_minus_f1_over_x3_mus(const double x, const double mu_s, const double mu_k) const
{
    return smooth_friction_f2_x_minus_f1_over_x3_mus(x, eps_v(), mu_s, mu_k);
}

} // namespace ipc