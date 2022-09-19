// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequlity constraints on a function.
#include <ipc/barrier/barrier.hpp>

#include <ipc/config.hpp>

namespace ipc {

double barrier_gradient(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    double b_grad = (dhat - d) * (2 * log(d / dhat) - dhat / d + 1);

#ifdef IPC_TOOLKIT_CONVERGENT
    b_grad /= dhat;
#endif

    return b_grad;
}

double barrier_hessian(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }

    const double dhat_d = dhat / d;

    double b_hess = (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;

#ifdef IPC_TOOLKIT_CONVERGENT
    b_hess /= dhat;
#endif

    return b_hess;
}

} // namespace ipc
