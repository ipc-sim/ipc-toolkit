// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequlity constraints on a function.
#include <ipc/barrier/barrier.hpp>

namespace ipc {

double barrier_gradient(double d, double dhat)
{
#ifdef IPC_TOOLKIT_CONVERGENT
    return physical_barrier_gradient(d, dhat);
#else
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    return (dhat - d) * (2 * log(d / dhat) - dhat / d + 1);
#endif
}

double barrier_hessian(double d, double dhat)
{
#ifdef IPC_TOOLKIT_CONVERGENT
    return physical_barrier_hessian(d, dhat);
#else
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }

    double dhat_d = dhat / d;

    return (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;
#endif
}

double physical_barrier_gradient(double d, double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    return (d - dhat) * (-2 * d * log(d / dhat) - d + dhat) / (d * dhat);
}

double physical_barrier_hessian(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double d_minus_dhat = d - dhat;
    return (-2 * log(d / dhat) + (d_minus_dhat * d_minus_dhat) / (d * d) - 4)
        / dhat
        + 4 / d;
}

} // namespace ipc
