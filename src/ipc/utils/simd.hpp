#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <Eigen/Core>
#include <xsimd/xsimd.hpp>

namespace Eigen {

/// @brief Let Eigen treat an xsimd batch as a scalar type.
///
/// This makes ``Eigen::Vector3<ipc::SimdBatch<double>>`` a valid argument to
/// the distance functions, evaluating one independent problem per SIMD lane.
template <typename T, typename A>
struct NumTraits<xsimd::batch<T, A>> : GenericNumTraits<xsimd::batch<T, A>> {
    using Real = xsimd::batch<T, A>;
    using NonInteger = xsimd::batch<T, A>;
    using Nested = xsimd::batch<T, A>;
    using Literal = xsimd::batch<T, A>;

    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        // A batch is not trivially default-constructible in the sense Eigen
        // means; make it initialize its coefficients.
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 1,
        MulCost = 1
    };

    static inline Real epsilon() { return Real(NumTraits<T>::epsilon()); }
    static inline Real dummy_precision()
    {
        return Real(NumTraits<T>::dummy_precision());
    }
    static inline int digits10() { return NumTraits<T>::digits10(); }
    static inline Real highest() { return Real(NumTraits<T>::highest()); }
    static inline Real lowest() { return Real(NumTraits<T>::lowest()); }
};

} // namespace Eigen

namespace ipc {

/// @brief A SIMD batch of ``T``, usable as the scalar type of a distance query.
///
/// Packing one independent problem per lane (a structure-of-arrays layout)
/// lets a single call evaluate ``SimdBatch<T>::size`` problems at once:
///
/// @code
/// Eigen::Vector3<ipc::SimdBatch<double>> p, e0, e1;
/// for (int k = 0; k < 3; ++k) {
///     p[k] = ipc::SimdBatch<double>::load_unaligned(&p_soa[k * n + i]);
///     // ... likewise for e0, e1
/// }
/// const auto d = ipc::point_edge_distance(p, e0, e1,
/// PointEdgeDistanceType::P_E);
/// @endcode
///
/// @warning **An explicit distance type is required.** ``AUTO`` and the
/// ``*_distance_type`` predicates are not available for batch scalars, and
/// throw ``std::invalid_argument`` if reached: the distance type is a per-lane
/// property, but the predicates return a single enum, so two lanes of the same
/// batch cannot report different closest features. Resolve the distance types
/// scalar-side and group the problems by type before batching, which is the
/// natural structure-of-arrays layout anyway.
///
/// @warning Only the distance *values* are instantiated. The gradients and
/// Hessians additionally require the ``autogen`` kernels to be instantiated for
/// the batch type.
///
/// @warning ``xsimd::default_arch`` is chosen from the compiler flags of each
/// translation unit, so a caller compiled without the SIMD flags the library
/// was built with names a *different* type than the one instantiated here, and
/// will fail to link. Build such callers with the same flags — CMake reports
/// them as ``SIMD_CXX_FLAGS``.
template <typename T> using SimdBatch = xsimd::batch<T, xsimd::default_arch>;

} // namespace ipc

#endif
