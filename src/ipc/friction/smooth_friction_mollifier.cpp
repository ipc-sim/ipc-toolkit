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

} // namespace ipc
