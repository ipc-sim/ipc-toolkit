#include "smooth_friction_mollifier.hpp"

#include <cassert>
#include <cmath>

namespace ipc {

double smooth_friction_f0(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return s;
    }
    return s * s * (1 - s / (3 * epsv)) / epsv + epsv / 3;
}

double smooth_friction_f1_over_x(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return 1 / s;
    }
    return (2 - s / epsv) / epsv;
}

double
smooth_friction_f2_x_minus_f1_over_x3(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return -1 / (s * s * s);
    }
    return -1 / (s * epsv * epsv);
}

} // namespace ipc
