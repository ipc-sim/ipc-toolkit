// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace ipc {

const Barrier& Barrier::get(Type type)
{
    static IPCBarrier ipc_barrier;
    static NormalizedBarrier normalized_barrier;
    static PhysicalBarrier physical_barrier;
    switch (type) {
    case IPC:
        return ipc_barrier;
    case NORMALIZED:
        return normalized_barrier;
    case PHYSICAL:
        return physical_barrier;
    default:
        throw std::runtime_error("Unknown barrier type");
    }
}

// =============================================================================

double normalized_barrier(const double d, const double dhat)
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

double normalized_barrier_first_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double t0 = 1.0 / dhat;
    const double t1 = d * t0;
    const double t2 = 1 - t1;
    return t2 * (2 * t0 * std::log(t1) - t2 / d);
}

double normalized_barrier_second_derivative(const double d, const double dhat)
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

double ipc_barrier(const double d, const double dhat)
{
    return dhat * dhat * normalized_barrier(d, dhat);
}

double ipc_barrier_first_derivative(const double d, const double dhat)
{
    return dhat * dhat * normalized_barrier_first_derivative(d, dhat);
}

double ipc_barrier_second_derivative(const double d, const double dhat)
{
    return dhat * dhat * normalized_barrier_second_derivative(d, dhat);
}

// =============================================================================

double physical_barrier(const double d, const double dhat)
{
    return dhat * normalized_barrier(d, dhat);
}

double physical_barrier_first_derivative(const double d, const double dhat)
{
    return dhat * normalized_barrier_first_derivative(d, dhat);
}

double physical_barrier_second_derivative(const double d, const double dhat)
{
    return dhat * normalized_barrier_second_derivative(d, dhat);
}

} // namespace ipc
