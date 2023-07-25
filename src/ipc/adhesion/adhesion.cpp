#include "adhesion.hpp"

#include <cassert>

namespace ipc {

// RCC: Pₙₖ(d) = ½Cₙ βₖ² dₖ²
// return 0.5 * adhession_stiffness * adhession_intensity *
// adhession_intensity
//     * distance_sqr;

// Fang and Li et al. [2023]:

double normal_adhesion_potential(
    const double d, const double dhat_p, const double dhat_a, const double a2)
{
    assert(d >= 0);
    assert(dhat_p < dhat_a);
    assert(a2 < 0);
    if (d < dhat_p) {
        const double a1 = a2 * (1 - dhat_a / dhat_p);
        const double c1 =
            a2 * (dhat_a - dhat_p) * (dhat_a - dhat_p) - dhat_p * dhat_p * a1;
        return a1 * d * d + c1;
    } else if (d < dhat_a) {
        const double b2 = -2 * a2 * dhat_a;
        const double c2 = a2 * dhat_a * dhat_a;
        return (a2 * d + b2) * d + c2;
    } else {
        return 0;
    }
}

double normal_adhesion_potential_gradient(
    const double d, const double dhat_p, const double dhat_a, const double a2)
{
    assert(d >= 0);
    assert(dhat_p < dhat_a);
    assert(a2 < 0);
    if (d < dhat_p) {
        const double a1 = a2 * (1 - dhat_a / dhat_p);
        return 2 * a1 * d;
    } else if (d < dhat_a) {
        const double b2 = -2 * a2 * dhat_a;
        return 2 * a2 * d + b2;
    } else {
        return 0;
    }
}

double normal_adhesion_potential_hessian(
    const double d, const double dhat_p, const double dhat_a, const double a2)
{
    assert(d >= 0);
    assert(dhat_p < dhat_a);
    assert(a2 < 0);
    if (d < dhat_p) {
        const double a1 = a2 * (1 - dhat_a / dhat_p);
        return 2 * a1;
    } else if (d < dhat_a) {
        return 2 * a2;
    } else {
        return 0;
    }
}

// RCC: Pₜₖ(‖uₖ‖) = ½Cₜ βₖ² ‖uₖ‖²
// return 0.5 * adhession_stiffness * adhession_intensity * adhession_intensity
//     * relative_sliding_displacement * relative_sliding_displacement;

// Fang and Li et al. [2023]:
// double tangential_adhesion_potential(
//     const double relative_sliding_displacement,
//     const double adhession_stiffness,
//     const double adhession_intensity)
// {
// }
} // namespace ipc
