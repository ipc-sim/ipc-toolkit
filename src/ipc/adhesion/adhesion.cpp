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
        const double b2 = -2 * a2 * dhat_a;
        return 2 * a2 * d + b2;
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

// -- Tangential Adhesion ------------------------------------------------------

double f0_t(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y <= 0) {
        return 0;
    } else if (y >= 2 * eps_a) {
        return 4 * eps_a / 3;
    }
    return y * y / eps_a * (1 - y / (3 * eps_a)); // -y³/(3ϵ²) + y²/ϵ
}

double f1_t(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y >= 2 * eps_a || y <= 0) {
        return 0;
    }

    return y / eps_a * (2 - y / eps_a); // -y²/ϵ² + 2y/ϵ
}

double df1_t(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y >= 2 * eps_a || y <= 0) {
        return 0;
    }

    return (-y / eps_a + 1) * 2 / eps_a; // -2y/ϵ² + 2/ϵ
}

double f1_t_over_x(const double y, const double eps_a)
{
    assert(eps_a > 0);
    if (y >= 2 * eps_a || y <= 0) {
        return 0;
    }

    return (2 - y / eps_a) / eps_a; // -y/ϵ² + 2/ϵ
}

double df1_t_x_minus_f1_t_over_x3(const double y, const double eps_a)
{
    assert(eps_a > 0);
    assert(y >= 0);
    if (y >= 2 * eps_a) {
        return 0;
    }
    return -1 / (y * eps_a * eps_a);
}

} // namespace ipc
