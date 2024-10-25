#include "smooth_friction_mollifier.hpp"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace ipc {

// ----------------------------------------------------------------------------
// Original friction mollifier functions

double f0_SF(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return s;
    }
    return s * s * (-s / (3 * epsv) + 1) / epsv + epsv / 3;
}

double f1_SF_over_x(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return 1 / s;
    }
    return (-s / epsv + 2) / epsv;
}

double df1_x_minus_f1_over_x3(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return -1 / (s * s * s);
    }
    return -1 / (s * epsv * epsv);
}

// ----------------------------------------------------------------------------
// Pairwise friction mollifier functions

double blend_mu(const double mu1, const double mu2)
{
    return (mu1 + mu2) / 2;
}

double f0_SF_pairwise(
    const double s, const double epsv, const double mu1, const double mu2)
{
    assert(epsv > 0);
    double blended_mu = blend_mu(mu1, mu2);

    if (std::abs(s) >= epsv) {
        return s * blended_mu;
    }
    return (s * s * (-s / (3 * epsv) + 1) / epsv + epsv / 3) * blended_mu;
}

double f1_SF_over_x_pairwise(
    const double s, const double epsv, const double mu1, const double mu2)
{
    assert(epsv > 0);
    double blended_mu = blend_mu(mu1, mu2);

    if (std::abs(s) >= epsv) {
        return blended_mu / s;
    }
    return (-s / epsv + 2) / epsv * blended_mu;
}

double df1_x_minus_f1_over_x3_pairwise(
    const double s, const double epsv, const double mu1, const double mu2)
{
    assert(epsv > 0);
    double blended_mu = blend_mu(mu1, mu2);

    if (std::abs(s) >= epsv) {
        return -blended_mu / (s * s * s);
    }
    return -blended_mu / (s * epsv * epsv);
}

// ----------------------------------------------------------------------------
// Pairwise friction mollifier functions with transition (no blending)

inline double select_mu(const double s, const double epsv, const double static_mu, const double kinetic_mu)
{
    // Use kinetic friction before transition (i.e., s < epsv), static friction after
    return (std::abs(s) < epsv) ? kinetic_mu : static_mu;
}

double f0_SF_pairwise_transition(
    const double s, const double epsv, const double static_mu, const double kinetic_mu)
{
    assert(epsv > 0);
    double mu = select_mu(s, epsv, static_mu, kinetic_mu);

    if (std::abs(s) >= epsv) {
        return s * mu;
    }
    return (s * s * (-s / (3 * epsv) + 1) / epsv + epsv / 3) * mu;
}

double f1_SF_over_x_pairwise_transition(
    const double s, const double epsv, const double static_mu, const double kinetic_mu)
{
    assert(epsv > 0);
    double mu = select_mu(s, epsv, static_mu, kinetic_mu);

    if (std::abs(s) >= epsv) {
        return mu / s;
    }
    return (-s / epsv + 2) / epsv * mu;
}

double df1_x_minus_f1_over_x3_pairwise_transition(
    const double s, const double epsv, const double static_mu, const double kinetic_mu)
{
    assert(epsv > 0);
    double mu = select_mu(s, epsv, static_mu, kinetic_mu);

    if (std::abs(s) >= epsv) {
        return -mu / (s * s * s);
    }
    return -mu / (s * epsv * epsv);
}

} // namespace ipc
