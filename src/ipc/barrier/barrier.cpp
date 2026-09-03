#include <ipc/config.hpp>

// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <ipc/math/math.hpp>

#include <cmath>
#include <limits>

namespace ipc {

template <typename T> T barrier(const T d, const T dhat)
{
    if (d <= T(0)) {
        return std::numeric_limits<T>::infinity();
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = -(d-d̂)²ln(d / d̂)
    const T d_minus_dhat = (d - dhat);
    return -d_minus_dhat * d_minus_dhat * log(d / dhat);
}

template <typename T> T barrier_first_derivative(const T d, const T dhat)
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    return (dhat - d) * (2 * log(d / dhat) - dhat / d + 1);
}

template <typename T> T barrier_second_derivative(const T d, const T dhat)
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    const T dhat_d = dhat / d;
    return (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;
}

// ============================================================================

template <typename T>
T ClampedLogSqBarrier<T>::operator()(const T d, const T dhat) const
{
    if (d <= T(0)) {
        return std::numeric_limits<T>::infinity();
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = (d-d̂)²ln²(d / d̂)
    const T d_minus_dhat = (d - dhat);
    const T log_d_dhat = log(d / dhat);
    return d_minus_dhat * d_minus_dhat * log_d_dhat * log_d_dhat;
}

template <typename T>
T ClampedLogSqBarrier<T>::first_derivative(const T d, const T dhat) const
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    // b(d) = (d - d̂)²ln²(d / d̂)
    // b'(d) = 2 (d - d̂) ln²(d / d̂) + 2 (d - d̂)² ln(d / d̂) / d
    //       = 2 (d - d̂) ln(d / d̂) [ln(d / d̂) + (d - d̂) / d]
    const T d_minus_dhat = (d - dhat);
    const T log_d_dhat = log(d / dhat);
    return T(2) * d_minus_dhat * log_d_dhat * (log_d_dhat + d_minus_dhat / d);
}

template <typename T>
T ClampedLogSqBarrier<T>::second_derivative(const T d, const T dhat) const
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    const T t0 = dhat - d;
    const T t1 = log(d / dhat);
    const T t2 = (t0 * t0) / (d * d);
    return T(2) * ((t1 * t1) - (t1 - T(1)) * t2 - T(4) * t1 * t0 / d);
}

// ============================================================================

template <typename T>
T CubicBarrier<T>::operator()(const T d, const T dhat) const
{
    if (d < dhat) {
        // b(d) = (d - d̂)³
        const T d_minus_dhat = (d - dhat);
        return -T(2) / T(3) / dhat * d_minus_dhat * d_minus_dhat * d_minus_dhat;
    } else {
        return T(0);
    }
}

template <typename T>
T CubicBarrier<T>::first_derivative(const T d, const T dhat) const
{
    if (d < dhat) {
        const T d_minus_dhat = (d - dhat);
        return T(-2) / dhat * d_minus_dhat * d_minus_dhat;
    } else {
        return T(0);
    }
}

template <typename T>
T CubicBarrier<T>::second_derivative(const T d, const T dhat) const
{
    if (d < dhat) {
        return T(4) * (T(1) - d / dhat);
    } else {
        return T(0);
    }
}

// ============================================================================

template <typename T>
T TwoStageBarrier<T>::operator()(const T d, const T dhat) const
{
    if (d >= dhat) {
        return T(0);
    } else if (d >= T(0.5) * dhat) {
        return T(0.5) * (dhat - d) * (dhat - d);
    } else {
        return T(-0.25) * dhat * dhat * (std::log(T(2) * d / dhat) - T(0.5));
    }
}

template <typename T>
T TwoStageBarrier<T>::first_derivative(const T d, const T dhat) const
{
    if (d >= dhat) {
        return T(0);
    } else if (d >= T(0.5) * dhat) {
        return d - dhat;
    } else {
        return T(-0.25) * dhat * dhat / d;
    }
}

template <typename T>
T TwoStageBarrier<T>::second_derivative(const T d, const T dhat) const
{
    if (d >= dhat) {
        return T(0);
    } else if (d >= T(0.5) * dhat) {
        return T(1);
    } else {
        return (T(0.25) * dhat * dhat) / (d * d);
    }
}

// ============================================================================

void InversePowerBarrier::h_and_derivs(
    const double d, const double dhat, double& h, double& dh, double& ddh)
{
    const double t = 2.0 * d / dhat;
    double B, dB, ddB;
    if (t < 1.0) {
        B = 2.0 / 3.0 - t * t + 0.5 * t * t * t;
        dB = -2.0 * t + 1.5 * t * t;
        ddB = -2.0 + 3.0 * t;
    } else if (t < 2.0) {
        const double s = 2.0 - t;
        B = s * s * s / 6.0;
        dB = -s * s / 2.0;
        ddB = s;
    } else {
        h = dh = ddh = 0.0;
        return;
    }
    // h = 2*B(t), dh/dd = 2*B'(t)*(dt/dd) = 2*B'*(2/dhat) = 4/dhat * B'
    h = 2.0 * B;
    dh = 4.0 / dhat * dB;
    ddh = 8.0 / (dhat * dhat) * ddB;
}

double InversePowerBarrier::operator()(const double d, const double dhat) const
{
    if (d <= 0.0)
        return std::numeric_limits<double>::infinity();
    double h, dh, ddh;
    h_and_derivs(d, dhat, h, dh, ddh);
    if (h == 0.0)
        return 0.0;
    return h / std::pow(d, m_power);
}

double
InversePowerBarrier::first_derivative(const double d, const double dhat) const
{
    if (d <= 0.0 || d >= dhat)
        return 0.0;
    double h, dh, ddh;
    h_and_derivs(d, dhat, h, dh, ddh);
    // b'(d) = (dh·d − p·h) / d^(p+1)
    return (dh * d - m_power * h) / std::pow(d, m_power + 1.0);
}

double
InversePowerBarrier::second_derivative(const double d, const double dhat) const
{
    if (d <= 0.0 || d >= dhat)
        return 0.0;
    double h, dh, ddh;
    h_and_derivs(d, dhat, h, dh, ddh);
    // b''(d) = (ddh·d² − 2p·dh·d + p(p+1)·h) / d^(p+2)
    const double d2 = d * d;
    return (ddh * d2 - 2.0 * m_power * dh * d + m_power * (m_power + 1.0) * h)
        / std::pow(d, m_power + 2.0);
}

// ============================================================================

double NearFarBarrier::near(const double d, const double dhat) const
{
    const double dhat_end = m_alpha * dhat;
    const double dhat_start = -dhat_end / 2.0;
    return (*m_base_barrier)(d, dhat)
        * (1.0 - Math<double>::smooth_heaviside(d, dhat_start, dhat_end));
}

double NearFarBarrier::far(const double d, const double dhat) const
{
    const double dhat_end = m_alpha * dhat;
    const double dhat_start = -dhat_end / 2.0;
    return (*m_base_barrier)(d, dhat)
        * Math<double>::smooth_heaviside(d, dhat_start, dhat_end);
}

double
NearFarBarrier::first_derivative_near(const double d, const double dhat) const
{
    const double dhat_end = m_alpha * dhat;
    const double dhat_start = -dhat_end / 2.0;
    const double b = (*m_base_barrier)(d, dhat);
    const double bp = m_base_barrier->first_derivative(d, dhat);
    const double w = Math<double>::smooth_heaviside(d, dhat_start, dhat_end);
    const double wp =
        Math<double>::smooth_heaviside_grad(d, dhat_start, dhat_end);
    return bp * (1.0 - w) - b * wp;
}

double
NearFarBarrier::first_derivative_far(const double d, const double dhat) const
{
    const double dhat_end = m_alpha * dhat;
    const double dhat_start = -dhat_end / 2.0;
    const double b = (*m_base_barrier)(d, dhat);
    const double bp = m_base_barrier->first_derivative(d, dhat);
    const double w = Math<double>::smooth_heaviside(d, dhat_start, dhat_end);
    const double wp =
        Math<double>::smooth_heaviside_grad(d, dhat_start, dhat_end);
    return bp * w + b * wp;
}

double
NearFarBarrier::second_derivative_near(const double d, const double dhat) const
{
    const double dhat_end = m_alpha * dhat;
    const double dhat_start = -dhat_end / 2.0;
    const double b = (*m_base_barrier)(d, dhat);
    const double bp = m_base_barrier->first_derivative(d, dhat);
    const double bpp = m_base_barrier->second_derivative(d, dhat);
    const double w = Math<double>::smooth_heaviside(d, dhat_start, dhat_end);
    const double wp =
        Math<double>::smooth_heaviside_grad(d, dhat_start, dhat_end);
    const double wpp =
        Math<double>::smooth_heaviside_hess(d, dhat_start, dhat_end);
    return bpp * (1.0 - w) - 2.0 * bp * wp - b * wpp;
}

double
NearFarBarrier::second_derivative_far(const double d, const double dhat) const
{
    const double dhat_end = m_alpha * dhat;
    const double dhat_start = -dhat_end / 2.0;
    const double b = (*m_base_barrier)(d, dhat);
    const double bp = m_base_barrier->first_derivative(d, dhat);
    const double bpp = m_base_barrier->second_derivative(d, dhat);
    const double w = Math<double>::smooth_heaviside(d, dhat_start, dhat_end);
    const double wp =
        Math<double>::smooth_heaviside_grad(d, dhat_start, dhat_end);
    const double wpp =
        Math<double>::smooth_heaviside_hess(d, dhat_start, dhat_end);
    return bpp * w + 2.0 * bp * wp + b * wpp;
}

// ============================================================================
// Explicit template instantiations
/// @cond DOXYGEN_SKIP
template class BarrierBase<float>;
template class BarrierBase<double>;
template class ClampedLogBarrier<float>;
template class ClampedLogBarrier<double>;
template class ClampedLogSqBarrier<float>;
template class ClampedLogSqBarrier<double>;
template class CubicBarrier<float>;
template class CubicBarrier<double>;
template class TwoStageBarrier<float>;
template class TwoStageBarrier<double>;
template float barrier(const float d, const float dhat);
template double barrier(const double d, const double dhat);
template float barrier_first_derivative(const float d, const float dhat);
template double barrier_first_derivative(const double d, const double dhat);
template float barrier_second_derivative(const float d, const float dhat);
template double barrier_second_derivative(const double d, const double dhat);
/// @endcond
// ============================================================================

} // namespace ipc
