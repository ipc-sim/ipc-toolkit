#include "tangential_adhesion_potential.hpp"

#include <ipc/adhesion/adhesion.hpp>

namespace ipc {

TangentialadhesionPotential::TangentialadhesionPotential(const double eps_a)
    : Super()
{
    set_eps_a(eps_a);
}

double TangentialadhesionPotential::f0(const double x) const
{
    return tangential_adhesion_f0(x, eps_a());
}

double TangentialadhesionPotential::f1_over_x(const double x) const
{
    return tangential_adhesion_f1_over_x(x, eps_a());
}

double TangentialadhesionPotential::f2_x_minus_f1_over_x3(const double x) const
{
    return tangential_adhesion_f2_x_minus_f1_over_x3(x, eps_a());
}

} // namespace ipc