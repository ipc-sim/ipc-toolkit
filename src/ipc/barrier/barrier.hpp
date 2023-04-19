// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.

#pragma once

namespace ipc {

/// @brief Function that grows to infinity as x approaches 0 from the right.
///
/// \f\[
///     b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)
/// \f\]
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The value of the barrier function at d.
double ipc_barrier(const double d, const double dhat);

/// @brief Derivative of the barrier function.
///
/// \f\[
///     b'(d) = (\hat{d}-d) \left( 2\ln\left( \frac{d}{\hat{d}} \right) -
///     \frac{\hat{d}}{d} + 1\right)
/// \f\]
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The derivative of the barrier wrt d.
double ipc_barrier_first_derivative(const double d, const double dhat);

/// @brief Second derivative of the barrier function.
///
/// \f\[
///     b''(d) = \left( \frac{\hat{d}}{d} + 2 \right) \frac{\hat{d}}{d} -
///     2\ln\left( \frac{d}{\hat{d}} \right) - 3
/// \f\]
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The second derivative of the barrier wrt d.
double ipc_barrier_second_derivative(const double d, const double dhat);

double normalized_barrier(const double d, const double dhat);
double normalized_barrier_first_derivative(const double d, const double dhat);
double normalized_barrier_second_derivative(const double d, const double dhat);

double physical_barrier(const double d, const double dhat);
double physical_barrier_first_derivative(const double d, const double dhat);
double physical_barrier_second_derivative(const double d, const double dhat);

/// Base class for barrier functions.
class Barrier {
public:
    Barrier() = default;

    virtual ~Barrier() = default;

    virtual double operator()(const double d, const double dhat) const = 0;
    virtual double
    first_derivative(const double d, const double dhat) const = 0;
    virtual double
    second_derivative(const double d, const double dhat) const = 0;

    enum Type { IPC, NORMALIZED, PHYSICAL };

    static const Barrier& get(Type type);
};

/// Barrier functions from [Li et al. 2020].
class IPCBarrier : public Barrier {
public:
    IPCBarrier() = default;

    double operator()(const double d, const double dhat) const override
    {
        return ipc_barrier(d, dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        return ipc_barrier_first_derivative(d, dhat);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        return ipc_barrier_second_derivative(d, dhat);
    }
};

class NormalizedBarrier : public Barrier {
public:
    NormalizedBarrier() = default;

    double operator()(const double d, const double dhat) const override
    {
        return normalized_barrier(d, dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        return normalized_barrier_first_derivative(d, dhat);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        return normalized_barrier_second_derivative(d, dhat);
    }
};

class PhysicalBarrier : public Barrier {
public:
    PhysicalBarrier() = default;

    double operator()(const double d, const double dhat) const override
    {
        return physical_barrier(d, dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        return physical_barrier_first_derivative(d, dhat);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        return physical_barrier_second_derivative(d, dhat);
    }
};

} // namespace ipc
