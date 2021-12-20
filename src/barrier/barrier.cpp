// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequlity constraints on a function.
#include <ipc/barrier/barrier.hpp>

namespace ipc {

double barrier_gradient(double d, double dhat)
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

double barrier_hessian(double d, double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }

    double dhat_d = dhat / d;

    return (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;
}

double physical_barrier_gradient(double d, double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    // b(d) = -d̂(d/d̂ - 1)²ln(d / d̂)
    // b'(d) = -2d̂(d/d̂ - 1)/d̂ ln(d / d̂) + -d̂(d/d̂ - 1)² (1/d)
    //       = -2(d/d̂ - 1)ln(d/d̂) - d̂(d/d̂ - 1)² (1/d)
    //       = -(d/d̂ - 1)(2ln(d/d̂) + d̂/d(d/d̂ - 1))
    const double tmp = (d / dhat - 1);
    return -tmp * (2 * log(d / dhat) + dhat / d * tmp);
}

} // namespace ipc
