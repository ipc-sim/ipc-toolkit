#include "friction_potential.hpp"

namespace ipc {

FrictionPotential::FrictionPotential(const double epsv) : Super()
{
    set_epsv(epsv);
}

double FrictionPotential::f0(const double x) const
{
    return smooth_friction_f0(x, epsv());
}

double FrictionPotential::f1_over_x(const double x) const
{
    return smooth_friction_f1_over_x(x, epsv());
}

double FrictionPotential::f2_x_minus_f1_over_x3(const double x) const
{
    return smooth_friction_f2_x_minus_f1_over_x3(x, epsv());
}

} // namespace ipc