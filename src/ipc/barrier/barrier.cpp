#include <ipc/config.hpp>

// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

namespace ipc {

// ============================================================================
// Free barrier functions -- shared between host C++ and CUDA device code.
// ============================================================================
//
// Each barrier is one select_lazy cascade, ordered by increasing d so it
// reads like the piecewise definition in the header. A scalar evaluates only
// the case it lands in -- so the log below is never reached for d <= 0 -- while
// a batch evaluates every case and blends per-lane, earlier cases winning.

template <typename T> IPC_TOOLKIT_HOST_DEVICE T barrier(const T d, const T dhat)
{
    using namespace ipc::numext; // log
    // b(d) = -(d-d̂)²ln(d / d̂)
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < dhat, [&] { return -sqr(d - dhat) * log(d / dhat); },
        [&] { return T(0); });
}

template <typename T>
IPC_TOOLKIT_HOST_DEVICE T barrier_first_derivative(const T d, const T dhat)
{
    using namespace ipc::numext; // log
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] { return (dhat - d) * (2 * log(d / dhat) - dhat / d + 1); },
        [&] { return T(0); });
}

template <typename T>
IPC_TOOLKIT_HOST_DEVICE T barrier_second_derivative(const T d, const T dhat)
{
    using namespace ipc::numext; // log
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T dhat_d = dhat / d;
            return (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;
        },
        [&] { return T(0); });
}
// ============================================================================
// Explicit template instantiations
/// @cond DOXYGEN_SKIP
#if IPC_TOOLKIT_INSTANTIATE_DEVICE_SCALARS
template float barrier(const float d, const float dhat);
template double barrier(const double d, const double dhat);
template float barrier_first_derivative(const float d, const float dhat);
template double barrier_first_derivative(const double d, const double dhat);
template float barrier_second_derivative(const float d, const float dhat);
template double barrier_second_derivative(const double d, const double dhat);
#endif
#ifdef IPC_TOOLKIT_WITH_SIMD
template SimdBatch<float>
barrier(const SimdBatch<float> d, const SimdBatch<float> dhat);
template SimdBatch<double>
barrier(const SimdBatch<double> d, const SimdBatch<double> dhat);
template SimdBatch<float>
barrier_first_derivative(const SimdBatch<float> d, const SimdBatch<float> dhat);
template SimdBatch<double> barrier_first_derivative(
    const SimdBatch<double> d, const SimdBatch<double> dhat);
template SimdBatch<float> barrier_second_derivative(
    const SimdBatch<float> d, const SimdBatch<float> dhat);
template SimdBatch<double> barrier_second_derivative(
    const SimdBatch<double> d, const SimdBatch<double> dhat);
#endif
/// @endcond

// ============================================================================
// Barrier class hierarchy -- host only.
// ============================================================================
//
// The classes are a virtual dispatch layer over the free functions above, and
// virtual dispatch cannot cross the host/device boundary: a vtable built on the
// host holds host code addresses, CUDA forbids passing an object of a class
// with virtual functions to a __global__ function, and BarrierPotential owns
// its barrier through a host-only std::shared_ptr. Skipping the hierarchy in
// the device pass also keeps every class symbol -- including the float and
// double ones -- in the host object, so each is emitted exactly once.
#ifndef __CUDACC__

template <typename T>
T ClampedLogSqBarrier<T>::operator()(const T d, const T dhat) const
{
    using namespace ipc::numext; // log
    // b(d) = (d-d̂)²ln²(d / d̂)
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < dhat,
        [&] {
            const T log_d_dhat = log(d / dhat);
            return sqr(d - dhat) * sqr(log_d_dhat);
        },
        [&] { return T(0); });
}

template <typename T>
T ClampedLogSqBarrier<T>::first_derivative(const T d, const T dhat) const
{
    using namespace ipc::numext; // log
    // b(d) = (d - d̂)²ln²(d / d̂)
    // b'(d) = 2 (d - d̂) ln²(d / d̂) + 2 (d - d̂)² ln(d / d̂) / d
    //       = 2 (d - d̂) ln(d / d̂) [ln(d / d̂) + (d - d̂) / d]
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T d_minus_dhat = (d - dhat);
            const T log_d_dhat = log(d / dhat);
            return T(2) * d_minus_dhat * log_d_dhat
                * (log_d_dhat + d_minus_dhat / d);
        },
        [&] { return T(0); });
}

template <typename T>
T ClampedLogSqBarrier<T>::second_derivative(const T d, const T dhat) const
{
    using namespace ipc::numext; // log
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T t0 = dhat - d;
            const T t1 = log(d / dhat);
            const T t2 = sqr(t0) / sqr(d);
            return T(2) * (sqr(t1) - (t1 - T(1)) * t2 - T(4) * t1 * t0 / d);
        },
        [&] { return T(0); });
}

// ============================================================================

template <typename T>
T CubicBarrier<T>::operator()(const T d, const T dhat) const
{
    // b(d) = (d - d̂)³
    //
    // A polynomial: finite at d <= 0, so unlike the log barriers it needs no
    // penetration guard.
    return select_lazy(
        d < dhat, [&] { return -T(2) / T(3) / dhat * cubic(d - dhat); },
        [&] { return T(0); });
}

template <typename T>
T CubicBarrier<T>::first_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d < dhat, [&] { return T(-2) / dhat * sqr(d - dhat); },
        [&] { return T(0); });
}

template <typename T>
T CubicBarrier<T>::second_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d < dhat, [&] { return T(4) * (T(1) - d / dhat); },
        [&] { return T(0); });
}

// ============================================================================

template <typename T>
T TwoStageBarrier<T>::operator()(const T d, const T dhat) const
{
    using namespace ipc::numext; // log
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < T(0.5) * dhat,
        [&] { return T(-0.25) * sqr(dhat) * (log(T(2) * d / dhat) - T(0.5)); },
        d < dhat, [&] { return T(0.5) * sqr(dhat - d); }, //
        [&] { return T(0); });
}

template <typename T>
T TwoStageBarrier<T>::first_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d <= T(0), [&] { return T(0); },                             //
        d < T(0.5) * dhat, [&] { return T(-0.25) * sqr(dhat) / d; }, //
        d < dhat, [&] { return d - dhat; },                          //
        [&] { return T(0); });
}

template <typename T>
T TwoStageBarrier<T>::second_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d <= T(0), [&] { return T(0); },                                   //
        d < T(0.5) * dhat, [&] { return (T(0.25) * sqr(dhat)) / sqr(d); }, //
        d < dhat, [&] { return T(1); },                                    //
        [&] { return T(0); });
}

// ============================================================================

template <typename T>
void InversePowerBarrier<T>::h_and_derivs(
    const T d, const T dhat, T& h, T& dh, T& ddh)
{
    const T t = T(2) * d / dhat;
    const T s = T(2) - t; // only meaningful on the second piece, 1 <= t < 2
    const T B = select_lazy(
        t < T(1), [&] { return T(2) / T(3) - sqr(t) + T(0.5) * cubic(t); },
        t < T(2), [&] { return cubic(s) / T(6); }, //
        [&] { return T(0); });
    const T dB = select_lazy(
        t < T(1), [&] { return T(-2) * t + T(1.5) * sqr(t); }, t < T(2),
        [&] { return -sqr(s) / T(2); }, //
        [&] { return T(0); });
    const T ddB = select_lazy(
        t < T(1), [&] { return T(-2) + T(3) * t; }, //
        t < T(2), [&] { return s; },                //
        [&] { return T(0); });
    // h = 2*B(t), dh/dd = 2*B'(t)*(dt/dd) = 2*B'*(2/dhat) = 4/dhat * B'
    h = T(2) * B;
    dh = T(4) / dhat * dB;
    ddh = T(8) / sqr(dhat) * ddB;
}

template <typename T>
T InversePowerBarrier<T>::operator()(const T d, const T dhat) const
{
    using namespace ipc::numext; // pow
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < dhat,
        [&] {
            T h, dh, ddh;
            h_and_derivs(d, dhat, h, dh, ddh);
            return h / pow(d, literal<T>(m_power));
        },
        [&] { return T(0); });
}

template <typename T>
T InversePowerBarrier<T>::first_derivative(const T d, const T dhat) const
{
    using namespace ipc::numext; // pow
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T p = literal<T>(m_power);
            T h, dh, ddh;
            h_and_derivs(d, dhat, h, dh, ddh);
            // b'(d) = (dh·d − p·h) / d^(p+1)
            return (dh * d - p * h) / pow(d, p + T(1));
        },
        [&] { return T(0); });
}

template <typename T>
T InversePowerBarrier<T>::second_derivative(const T d, const T dhat) const
{
    using namespace ipc::numext; // pow
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T p = literal<T>(m_power);
            T h, dh, ddh;
            h_and_derivs(d, dhat, h, dh, ddh);
            // b''(d) = (ddh·d² − 2p·dh·d + p(p+1)·h) / d^(p+2)
            return (ddh * sqr(d) - T(2) * p * dh * d + p * (p + T(1)) * h)
                / pow(d, p + T(2));
        },
        [&] { return T(0); });
}

template <typename T> T InversePowerBarrier<T>::units(const T dhat) const
{
    using namespace ipc::numext; // pow
    return T(1) / pow(dhat, literal<T>(m_power));
}

// ============================================================================

namespace {
    // The C² piecewise-cubic step behind NearFarBarrier's far weight, rising
    // from 0 at x = -3 to 1 at x = 0, and its first two derivatives. These are
    // the same polynomials as Math::smooth_heaviside, but we write them out as
    // select_lazy cascades because Math branches with `if`, which a batch
    // cannot do since each lane may land in a different piece.

    template <typename T> T smooth_step(const T x)
    {
        return select_lazy(
            x <= T(-3), [&] { return T(0); },                   //
            x <= T(-2), [&] { return cubic(T(3) + x) / T(6); }, //
            x <= T(-1),
            [&] { return (((T(-2) * x - T(9)) * x - T(9)) * x + T(3)) / T(6); },
            x < T(0), [&] { return cubic(x) / T(6) + T(1); }, //
            [&] { return T(1); });
    }

    template <typename T> T smooth_step_first_derivative(const T x)
    {
        return select_lazy(
            x <= T(-3), [&] { return T(0); },                 //
            x <= T(-2), [&] { return sqr(T(3) + x) / T(2); }, //
            x <= T(-1), [&] { return -(sqr(x) + T(3) * x + T(1.5)); }, x < T(0),
            [&] { return sqr(x) / T(2); }, //
            [&] { return T(0); });
    }

    template <typename T> T smooth_step_second_derivative(const T x)
    {
        return select_lazy(
            x <= T(-3), [&] { return T(0); },             //
            x <= T(-2), [&] { return T(3) + x; },         //
            x <= T(-1), [&] { return T(-3) - T(2) * x; }, //
            x < T(0), [&] { return x; },                  //
            [&] { return T(0); });
    }
} // namespace

// We map the far weight's ramp [alpha * dhat / 2, alpha * dhat] onto the
// standard step's [-3, 0] with x = (d - alpha * dhat) * scale, where
// scale = 6 / (alpha * dhat). Each derivative of H then picks up one factor of
// that scale from the chain rule.

template <typename T>
auto NearFarBarrier<T>::is_purely_near(const T d, const T dhat) const -> Mask
{
    // The far weight H and all its derivatives are exactly zero at or below
    // alpha * dhat / 2, so near = b and far = 0 there. We return those values
    // directly instead of multiplying through, because the base barrier is +inf
    // at d <= 0 (and InversePowerBarrier overflows to +inf for tiny d > 0), and
    // inf * 0 would turn an exact contact into NaN.
    return d <= literal<T>(0.5 * m_alpha) * dhat;
}

template <typename T>
T NearFarBarrier<T>::far_weight(const T d, const T dhat) const
{
    const T end = literal<T>(m_alpha) * dhat;
    return smooth_step((d - end) * (T(6) / end));
}

template <typename T>
T NearFarBarrier<T>::far_weight_first_derivative(const T d, const T dhat) const
{
    const T end = literal<T>(m_alpha) * dhat;
    const T scale = T(6) / end;
    return smooth_step_first_derivative((d - end) * scale) * scale;
}

template <typename T>
T NearFarBarrier<T>::far_weight_second_derivative(const T d, const T dhat) const
{
    const T end = literal<T>(m_alpha) * dhat;
    const T scale = T(6) / end;
    return smooth_step_second_derivative((d - end) * scale) * sqr(scale);
}

template <typename T>
T NearFarBarrier<T>::near_value(const T d, const T dhat) const
{
    const T b = (*m_base_barrier)(d, dhat);
    return select_lazy(
        is_purely_near(d, dhat), [&] { return b; },
        [&] { return b * (T(1) - far_weight(d, dhat)); });
}

template <typename T>
T NearFarBarrier<T>::far_value(const T d, const T dhat) const
{
    return select_lazy(
        is_purely_near(d, dhat), [&] { return T(0); },
        [&] { return (*m_base_barrier)(d, dhat) * far_weight(d, dhat); });
}

template <typename T>
T NearFarBarrier<T>::first_derivative_near(const T d, const T dhat) const
{
    const T bp = m_base_barrier->first_derivative(d, dhat);
    return select_lazy(
        is_purely_near(d, dhat), [&] { return bp; },
        [&] {
            const T b = (*m_base_barrier)(d, dhat);
            const T w = far_weight(d, dhat);
            const T wp = far_weight_first_derivative(d, dhat);
            return bp * (T(1) - w) - b * wp;
        });
}

template <typename T>
T NearFarBarrier<T>::first_derivative_far(const T d, const T dhat) const
{
    return select_lazy(
        is_purely_near(d, dhat), [&] { return T(0); },
        [&] {
            const T b = (*m_base_barrier)(d, dhat);
            const T bp = m_base_barrier->first_derivative(d, dhat);
            const T w = far_weight(d, dhat);
            const T wp = far_weight_first_derivative(d, dhat);
            return bp * w + b * wp;
        });
}

template <typename T>
T NearFarBarrier<T>::second_derivative_near(const T d, const T dhat) const
{
    const T bpp = m_base_barrier->second_derivative(d, dhat);
    return select_lazy(
        is_purely_near(d, dhat), [&] { return bpp; },
        [&] {
            const T b = (*m_base_barrier)(d, dhat);
            const T bp = m_base_barrier->first_derivative(d, dhat);
            const T w = far_weight(d, dhat);
            const T wp = far_weight_first_derivative(d, dhat);
            const T wpp = far_weight_second_derivative(d, dhat);
            return bpp * (T(1) - w) - T(2) * bp * wp - b * wpp;
        });
}

template <typename T>
T NearFarBarrier<T>::second_derivative_far(const T d, const T dhat) const
{
    return select_lazy(
        is_purely_near(d, dhat), [&] { return T(0); },
        [&] {
            const T b = (*m_base_barrier)(d, dhat);
            const T bp = m_base_barrier->first_derivative(d, dhat);
            const T bpp = m_base_barrier->second_derivative(d, dhat);
            const T w = far_weight(d, dhat);
            const T wp = far_weight_first_derivative(d, dhat);
            const T wpp = far_weight_second_derivative(d, dhat);
            return bpp * w + T(2) * bp * wp + b * wpp;
        });
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
template class InversePowerBarrier<float>;
template class InversePowerBarrier<double>;
template class NearFarBarrier<float>;
template class NearFarBarrier<double>;
#ifdef IPC_TOOLKIT_WITH_SIMD
template class BarrierBase<SimdBatch<float>>;
template class BarrierBase<SimdBatch<double>>;
template class ClampedLogBarrier<SimdBatch<float>>;
template class ClampedLogBarrier<SimdBatch<double>>;
template class ClampedLogSqBarrier<SimdBatch<float>>;
template class ClampedLogSqBarrier<SimdBatch<double>>;
template class CubicBarrier<SimdBatch<float>>;
template class CubicBarrier<SimdBatch<double>>;
template class TwoStageBarrier<SimdBatch<float>>;
template class TwoStageBarrier<SimdBatch<double>>;
template class InversePowerBarrier<SimdBatch<float>>;
template class InversePowerBarrier<SimdBatch<double>>;
template class NearFarBarrier<SimdBatch<float>>;
template class NearFarBarrier<SimdBatch<double>>;
#endif
/// @endcond
// ============================================================================

#endif // !__CUDACC__

} // namespace ipc
