// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <cmath>
#include <limits>

namespace ipc {

double barrier(const double d, const double dhat)
{
    if (d <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (d >= dhat) {
        return 0;
    }
    // b(d) = -(d-d̂)²ln(d / d̂)
    const double d_minus_dhat = (d - dhat);
    return -d_minus_dhat * d_minus_dhat * log(d / dhat);
}

double barrier_first_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    return (dhat - d) * (2 * log(d / dhat) - dhat / d + 1);
}

double barrier_second_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double dhat_d = dhat / d;
    return (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;
}

} // namespace ipc
