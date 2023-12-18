#include "friction_potential.hpp"

namespace ipc {

FrictionPotential::FrictionPotential(const double epsv) : Super()
{
    set_epsv(epsv);
}

double FrictionPotential::f0(const double x) const { return f0_SF(x, epsv()); }

double FrictionPotential::f1_over_x(const double x) const
{
    return f1_SF_over_x(x, epsv());
}

double FrictionPotential::df1_x_minus_f1_over_x3(const double x) const
{
    return df1_SF_x_minus_f1_SF_over_x3(x, epsv());
}

} // namespace ipc