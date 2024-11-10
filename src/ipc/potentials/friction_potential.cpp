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

} // namespace ipc