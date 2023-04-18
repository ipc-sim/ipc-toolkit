// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.

#pragma once

namespace ipc {

/// Base class for barrier functions.
class Barrier {
public:
    Barrier(double dhat) : dhat(dhat) { }
    virtual ~Barrier() = default;

    virtual double operator()(const double d) const = 0;
    virtual double first_derivative(const double d) const = 0;
    virtual double second_derivative(const double d) const = 0;

    enum Type { IPC, NORMALIZED, PHYISCAL };

    static std::shared_ptr<Barrier> get(Type type, double dhat);

    double dhat;
};

/// Barrier functions from [Li et al. 2020].
class IPCBarrier : public Barrier {
public:
    OriginalBarrier(double dhat) : Barrier(dhat) { }

    double operator()(const double d) const override { return value(d, dhat); }

    double first_derivative(const double d) const override
    {
        return first_derivative(d, dhat);
    }

    double second_derivative(const double d) const override
    {
        return second_derivative(d, dhat);
    }

    /// @brief Function that grows to infinity as x approaches 0 from the right.
    ///
    /// \f\[
    ///     b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    static double value(const double d, const double dhat);

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
    static double first_derivative(const double d, const double dhat);

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
    static double second_derivative(const double d, const double dhat);
};

class NormalizedBarrier : public Barrier {
public:
    NormalizedBarrier(double dhat) : Barrier(dhat) { }

    double operator()(const double d) const override { return value(d, dhat); }

    double first_derivative(const double d) const override
    {
        return first_derivative(d, dhat);
    }

    double second_derivative(const double d) const override
    {
        return second_derivative(d, dhat);
    }

    static double value(const double d, const double dhat);
    static double first_derivative(const double d, const double dhat);
    static double second_derivative(const double d, const double dhat);
};

class PhysicalBarrier : public Barrier {
public:
    PhysicalBarrier(double dhat) : Barrier(dhat) { }

    double operator()(const double d) const override { return value(d, dhat); }

    double first_derivative(const double d) const override
    {
        return first_derivative(d, dhat);
    }

    double second_derivative(const double d) const override
    {
        return second_derivative(d, dhat);
    }

    static double value(const double d, const double dhat);
    static double first_derivative(const double d, const double dhat);
    static double second_derivative(const double d, const double dhat);
};

} // namespace ipc
