#include "smooth_friction_mollifier.hpp"

#include <cassert>
#include <cmath>

namespace ipc {

double f0_SF(const double x, const double epsv_times_h)
{
    assert(epsv_times_h > 0);
    if (std::abs(x) >= epsv_times_h) {
        return x;
    }
    return x * x * (-x / (3 * epsv_times_h) + 1) / epsv_times_h
        + epsv_times_h / 3;
}

double f1_SF_over_x(const double x, const double epsv_times_h)
{
    assert(epsv_times_h > 0);
    if (std::abs(x) >= epsv_times_h) {
        return 1 / x;
    }
    return (-x / epsv_times_h + 2) / epsv_times_h;
}

double df1_x_minus_f1_over_x3(const double x, const double epsv_times_h)
{
    assert(epsv_times_h > 0);
    if (std::abs(x) >= epsv_times_h) {
        return -1 / (x * x * x);
    }
    return -1 / (x * epsv_times_h * epsv_times_h);
}

} // namespace ipc
