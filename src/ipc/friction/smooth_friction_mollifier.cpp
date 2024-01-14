#include "smooth_friction_mollifier.hpp"

#include <cassert>
#include <cmath>

namespace ipc {

double f0_SF(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return s;
    }
    return s * s * (1 - s / (3 * epsv)) / epsv + epsv / 3;
}

double f1_SF_over_x(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return 1 / s;
    }
    return (2 - s / epsv) / epsv;
}

double df1_SF_x_minus_f1_SF_over_x3(const double s, const double epsv)
{
    assert(epsv > 0);
    if (std::abs(s) >= epsv) {
        return -1 / (s * s * s);
    }
    return -1 / (s * epsv * epsv);
}

} // namespace ipc
