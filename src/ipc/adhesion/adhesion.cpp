// Adhesion model of Fang and Li et al. [2023].

#include "adhesion.hpp"

#include <cassert>

namespace ipc {

// -- Normal Adhesion ----------------------------------------------------------

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

double normal_adhesion_potential_first_derivative(
    const double d, const double dhat_p, const double dhat_a, const double a2)
{
    assert(d >= 0);
    assert(dhat_p < dhat_a);
    assert(a2 < 0);
    if (d < dhat_p) {
        const double a1 = a2 * (1 - dhat_a / dhat_p);
        return 2 * a1 * d;
    } else if (d < dhat_a) {
        // const double b2 = -2 * a2 * dhat_a;
        return 2 * a2 * (d - dhat_a);
    } else {
        return 0;
    }
}

double normal_adhesion_potential_second_derivative(
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

double max_normal_adhesion_force_magnitude(
    const double dhat_p, const double dhat_a, const double a2)
{
    assert(dhat_p < dhat_a);
    assert(a2 < 0);
    // max_d a' = a'(d̂ₚ) = 2a₂ (d̂ₚ - d̂ₐ)
    return 2 * a2 * (dhat_p - dhat_a);
}

// -- Tangential Adhesion ------------------------------------------------------

double tangential_adhesion_f0(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y <= 0) {
        return 0;
    } else if (y >= 2 * eps_a) {
        return 4 * eps_a / 3;
    }
    return y * y / eps_a * (1 - y / (3 * eps_a)); // -y³/(3ϵ²) + y²/ϵ
}

double tangential_adhesion_f1(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y >= 2 * eps_a || y <= 0) {
        return 0;
    }

    const double y_over_eps_a = y / eps_a;
    return y_over_eps_a * (2 - y_over_eps_a); // -y²/ϵ² + 2y/ϵ
}

double tangential_adhesion_f2(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y >= 2 * eps_a || y <= 0) {
        return 0;
    }

    return (-y / eps_a + 1) * 2 / eps_a; // -2y/ϵ² + 2/ϵ
}

double tangential_adhesion_f1_over_x(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y >= 2 * eps_a || y <= 0) {
        return 0;
    }

    return (2 - y / eps_a) / eps_a; // -y/ϵ² + 2/ϵ
}

double
tangential_adhesion_f2_x_minus_f1_over_x3(const double y, const double eps_a)
{
    assert(eps_a > 0);
    assert(y >= 0);
    if (y >= 2 * eps_a) {
        return 0;
    }
    return -1 / (y * eps_a * eps_a);
}

} // namespace ipc
