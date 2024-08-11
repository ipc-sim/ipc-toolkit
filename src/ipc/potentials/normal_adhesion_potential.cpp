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

double NormalAdhesionPotential::force_magnitude(
    const double, const double dmin, const double) const
{
    const auto [arg_dhat_p, arg_dhat_a, a2] =
        normal_adhesion_potential_args(dmin);
    return max_normal_adhesion_force_magnitude(arg_dhat_p, arg_dhat_a, a2);
}

VectorMax12d NormalAdhesionPotential::force_magnitude_gradient(
    const double,
    const VectorMax12d& distance_squared_gradient,
    const double,
    const double) const
{
    // The force magnitude is constant wrt positions, so the gradient is zero.
    return VectorMax12d::Zero(distance_squared_gradient.size());
}

double NormalAdhesionPotential::operator()(
    const double distance_squared, const double dmin) const
{
    const auto [arg_dhat_p, arg_dhat_a, a2] =
        normal_adhesion_potential_args(dmin);
    return normal_adhesion_potential(
        distance_squared - dmin * dmin, arg_dhat_p, arg_dhat_a, a2);
}

double NormalAdhesionPotential::gradient(
    const double distance_squared, const double dmin) const
{
    const auto [arg_dhat_p, arg_dhat_a, a2] =
        normal_adhesion_potential_args(dmin);
    return normal_adhesion_potential_first_derivative(
        distance_squared - dmin * dmin, arg_dhat_p, arg_dhat_a, a2);
}

double NormalAdhesionPotential::hessian(
    const double distance_squared, const double dmin) const
{
    const auto [arg_dhat_p, arg_dhat_a, a2] =
        normal_adhesion_potential_args(dmin);
    return normal_adhesion_potential_second_derivative(
        distance_squared - dmin * dmin, arg_dhat_p, arg_dhat_a, a2);
}

std::array<double, 3>
NormalAdhesionPotential::normal_adhesion_potential_args(const double dmin) const
{
    const double arg_dhat_p = 2 * dmin * dhat_p + dhat_p * dhat_p;
    const double arg_dhat_a = 2 * dmin * dhat_a + dhat_a * dhat_a;
    const double a2 =
        Y * eps_c / (4 * (dhat_p + dmin) * (arg_dhat_p - arg_dhat_a));
    return { { arg_dhat_p, arg_dhat_a, a2 } };
}

} // namespace ipc