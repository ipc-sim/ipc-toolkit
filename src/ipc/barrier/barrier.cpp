#include <ipc/config.hpp>

// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

namespace ipc {

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
        d < dhat, [&] { return -sqr(d - dhat) * ipc::numext::log(d / dhat); },
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
        [&] {
            return (dhat - d) * (2 * ipc::numext::log(d / dhat) - dhat / d + 1);
        },
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
            return (dhat_d + 2) * dhat_d - 2 * ipc::numext::log(d / dhat) - 3;
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

} // namespace ipc
