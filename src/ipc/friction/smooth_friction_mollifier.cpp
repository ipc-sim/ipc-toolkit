// smooth_friction_mollifier.cpp

#include "smooth_friction_mollifier.hpp"

#include <cassert>
#include <cmath>

namespace ipc {

double smooth_friction_f0(const double y, const double eps_v)
{
    assert(eps_v > 0);
    if (std::abs(y) >= eps_v) {
        return y;
    }
    return y * y * (1 - y / (3 * eps_v)) / eps_v + eps_v / 3;
}

double smooth_friction_f1(const double y, const double eps_v)
{
    assert(eps_v > 0);
    if (std::abs(y) >= eps_v) {
        return 1;
    }
    return y * (2 - y / eps_v) / eps_v;
}

double smooth_friction_f2(const double y, const double eps_v)
{
    assert(eps_v > 0);
    if (std::abs(y) >= eps_v) {
        return 0;
    }
    return (2 - 2 * y / eps_v) / eps_v;
}

double smooth_friction_f1_over_x(const double y, const double eps_v)
{
    assert(eps_v > 0);
    if (std::abs(y) >= eps_v) {
        return 1 / y;
    }
    return (2 - y / eps_v) / eps_v;
}

double smooth_friction_f2_x_minus_f1_over_x3(const double y, const double eps_v)
{
    assert(eps_v > 0);
    if (std::abs(y) >= eps_v) {
        return -1 / (y * y * y);
    }
    return -1 / (y * eps_v * eps_v);
}

double smooth_friction_mus(const double y, const double eps_v, const double mu_s, const double mu_k)
{
    assert(eps_v > 0);
    assert(mu_s > 0);
    assert(mu_k > 0);

    double abs_y = std::abs(y);

    if (abs_y >= eps_v) {
        return mu_k;
    }

    double normalized_y = abs_y / eps_v;
    // Cubic Hermite polynomial for smooth transition between mu_s and mu_k
    double transition = 3.0 * std::pow(normalized_y, 2) - 2.0 * std::pow(normalized_y, 3);
    return (mu_k - mu_s) * transition + mu_s;
}

double smooth_friction_f0_mus(const double y, const double eps_v, const double mu_s, const double mu_k)
{
    double mu = smooth_friction_mus(y, eps_v, mu_s, mu_k);
    double f0 = smooth_friction_f0(y, eps_v);
    return mu * f0;
}

double smooth_friction_f1_mus(const double y, const double eps_v, const double mu_s, const double mu_k)
{
    double mu = smooth_friction_mus(y, eps_v, mu_s, mu_k);
    double f1 = smooth_friction_f1(y, eps_v);
    return mu * f1;
}

double smooth_friction_f2_mus(const double y, const double eps_v, const double mu_s, const double mu_k)
{
    double mu = smooth_friction_mus(y, eps_v, mu_s, mu_k);
    double f2 = smooth_friction_f2(y, eps_v);
    return mu * f2;
}

double smooth_friction_f1_over_x_mus(const double y, const double eps_v, const double mu_s, const double mu_k)
{
    if (y == 0.0) {
        return smooth_friction_mus(y, eps_v, mu_s, mu_k) * (2.0 / eps_v);
    }

    double mu = smooth_friction_mus(y, eps_v, mu_s, mu_k);
    double f1_over_x = smooth_friction_f1_over_x(y, eps_v);
    return mu * f1_over_x;
}

double smooth_friction_f2_x_minus_f1_over_x3_mus(const double y, const double eps_v, const double mu_s, const double mu_k)
{
    if (y == 0.0) {
        return smooth_friction_mus(y, eps_v, mu_s, mu_k) * (-1.0 / (eps_v * eps_v));
    }

    double mu = smooth_friction_mus(y, eps_v, mu_s, mu_k);
    double f2_x_minus_f1_over_x3 = smooth_friction_f2_x_minus_f1_over_x3(y, eps_v);
    return mu * f2_x_minus_f1_over_x3;
}

} // namespace ipc
