// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <cmath>
#include <limits>

namespace ipc {

std::shared_ptr<Barrier> Barrier::get(Type type, double dhat)
{
    switch (type) {
    case IPC:
        return std::make_shared<IPCBarrier>(dhat);
    case NORMALIZED:
        return std::make_shared<NormalizedBarrier>(dhat);
    case PHYISCAL:
        return std::make_shared<PhysicalBarrier>(dhat);
    default:
        throw std::runtime_error("Unknown barrier type");
    }
}

// =============================================================================

double NormalizedBarrier::value(const double d, const double dhat)
{
    if (d <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (d >= dhat) {
        return 0;
    }

    // b(d) = -(d/d̂ - 1)²ln(d / d̂)
    const auto t0 = d / dhat;
    return -std::pow(1 - t0, 2) * std::log(t0);
}

double NormalizedBarrier::first_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double t0 = 1.0 / dhat;
    const double t1 = d * t0;
    const double t2 = 1 - t1;
    return t2 * (2 * t0 * std::log(t1) - t2 / d);
}

double NormalizedBarrier::second_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }

    const double t0 = 1.0 / dhat;
    const double t1 = d * t0;
    const double t2 = 1 - t1;
    return 4 * t0 * t2 / d + std::pow(t2, 2) / std::pow(d, 2)
        - 2 * std::log(t1) / std::pow(dhat, 2);
}

// =============================================================================

double IPCBarrier::barrier(const double d, const double dhat)
{
    return dhat * dhat * NormalizedBarrier::value(d, dhat);
}

double IPCBarrier::barrier_gradient(const double d, const double dhat)
{
    return dhat * dhat * NormalizedBarrier::first_derivative(d, dhat);
}

double IPCBarrier::barrier_hessian(const double d, const double dhat)
{
    return dhat * dhat * NormalizedBarrier::second_derivative(d, dhat);
}

// =============================================================================

double PhysicalBarrier::barrier(const double d, const double dhat)
{
    return dhat * NormalizedBarrier::value(d, dhat);
}

double PhysicalBarrier::barrier_gradient(const double d, const double dhat)
{
    return dhat * NormalizedBarrier::first_derivative(d, dhat);
}

double PhysicalBarrier::barrier_hessian(const double d, const double dhat)
{
    return dhat * NormalizedBarrier::second_derivative(d, dhat);
}

} // namespace ipc
