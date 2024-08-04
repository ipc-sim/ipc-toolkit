#include "normal_adhesion_potential.hpp"

#include <ipc/adhesion/adhesion.hpp>

namespace ipc {

NormalAdhesionPotential::NormalAdhesionPotential(
    const double _dhat_p,
    const double _dhat_a,
    const double _Y,
    const double _eps_c)
    : dhat_p(_dhat_p)
    , dhat_a(_dhat_a)
    , Y(_Y)
    , eps_c(_eps_c)
{
}

double NormalAdhesionPotential::operator()(
    const double distance_sqr, const double dmin) const
{
    const double arg_dhat_p = 2 * dmin * dhat_p + dhat_p * dhat_p;
    const double arg_dhat_a = 2 * dmin * dhat_a + dhat_a * dhat_a;
    const double a2 =
        Y * eps_c / (4 * (dhat_p + dmin) * (arg_dhat_p - arg_dhat_a));
    return normal_adhesion_potential(
        distance_sqr - dmin * dmin, arg_dhat_p, arg_dhat_a, a2);
}

double NormalAdhesionPotential::gradient(
    const double distance_sqr, const double dmin) const
{
    const double arg_dhat_p = 2 * dmin * dhat_p + dhat_p * dhat_p;
    const double arg_dhat_a = 2 * dmin * dhat_a + dhat_a * dhat_a;
    const double a2 =
        Y * eps_c / (4 * (dhat_p + dmin) * (arg_dhat_p - arg_dhat_a));
    return normal_adhesion_potential_first_derivative(
        distance_sqr - dmin * dmin, arg_dhat_p, arg_dhat_a, a2);
}

double NormalAdhesionPotential::hessian(
    const double distance_sqr, const double dmin) const
{
    const double arg_dhat_p = 2 * dmin * dhat_p + dhat_p * dhat_p;
    const double arg_dhat_a = 2 * dmin * dhat_a + dhat_a * dhat_a;
    const double a2 =
        Y * eps_c / (4 * (dhat_p + dmin) * (arg_dhat_p - arg_dhat_a));
    return normal_adhesion_potential_second_derivative(
        distance_sqr - dmin * dmin, arg_dhat_p, arg_dhat_a, a2);
}

} // namespace ipc