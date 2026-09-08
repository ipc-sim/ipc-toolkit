// Pieces of the per-scalar-type barrier potential benchmark that are shared
// between its host translation unit (benchmark_simd_barrier_potential.cpp) and
// its CUDA one (benchmark_cuda_barrier_potential.cu).
//
// The CUDA interface below takes raw pointers rather than the packed buffers'
// own container type, which keeps that type -- and the xsimd allocator behind
// it -- out of the CUDA translation unit. nvcc compiles xsimd, so this is no
// longer forced, but a batch scalar still has no device-callable operations:
// the .cu instantiates float and double only, and the narrower interface is
// what keeps it that way.

#pragma once

#include <ipc/config.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace ipc::tests::bench {

// -- Collision kinds ----------------------------------------------------------

enum class Kind : uint8_t { VV, EV, EE, FV };

constexpr int num_vertices(const Kind k)
{
    switch (k) {
    case Kind::VV:
        return 2;
    case Kind::EV:
        return 3;
    default:
        return 4;
    }
}

inline const char* kind_name(const Kind k)
{
    switch (k) {
    case Kind::VV:
        return "vertex-vertex";
    case Kind::EV:
        return "edge-vertex";
    case Kind::EE:
        return "edge-edge";
    default:
        return "face-vertex";
    }
}

/// @brief The distance query for a kind, for any scalar `T`.
///
/// Annotated for both execution spaces: the CPU variants call these from a
/// host loop and the CUDA variants from a kernel, on the same instantiations.
template <typename T, Kind K> struct Stencil;

template <typename T> struct Stencil<T, Kind::VV> {
    static constexpr int NV = 2;
    using V = Eigen::Vector3<T>;
    IPC_TOOLKIT_HOST_DEVICE static T distance(const V* x, int /*dtype*/)
    {
        return point_point_distance(x[0], x[1]);
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Vector<T, 6>
    gradient(const V* x, int /*dtype*/)
    {
        return point_point_distance_gradient(x[0], x[1]);
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Matrix<T, 6, 6>
    hessian(const V* x, int /*dtype*/)
    {
        return point_point_distance_hessian(x[0], x[1]);
    }
};

template <typename T> struct Stencil<T, Kind::EV> {
    static constexpr int NV = 3;
    using V = Eigen::Vector3<T>;
    IPC_TOOLKIT_HOST_DEVICE static T distance(const V* x, int dtype)
    {
        return point_edge_distance(
            x[0], x[1], x[2], PointEdgeDistanceType(dtype));
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Vector<T, 9>
    gradient(const V* x, int dtype)
    {
        return point_edge_distance_gradient(
            x[0], x[1], x[2], PointEdgeDistanceType(dtype));
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Matrix<T, 9, 9>
    hessian(const V* x, int dtype)
    {
        return point_edge_distance_hessian(
            x[0], x[1], x[2], PointEdgeDistanceType(dtype));
    }
};

template <typename T> struct Stencil<T, Kind::EE> {
    static constexpr int NV = 4;
    using V = Eigen::Vector3<T>;
    IPC_TOOLKIT_HOST_DEVICE static T distance(const V* x, int dtype)
    {
        return edge_edge_distance(
            x[0], x[1], x[2], x[3], EdgeEdgeDistanceType(dtype));
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Vector<T, 12>
    gradient(const V* x, int dtype)
    {
        return edge_edge_distance_gradient(
            x[0], x[1], x[2], x[3], EdgeEdgeDistanceType(dtype));
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Matrix<T, 12, 12>
    hessian(const V* x, int dtype)
    {
        return edge_edge_distance_hessian(
            x[0], x[1], x[2], x[3], EdgeEdgeDistanceType(dtype));
    }
};

template <typename T> struct Stencil<T, Kind::FV> {
    static constexpr int NV = 4;
    using V = Eigen::Vector3<T>;
    IPC_TOOLKIT_HOST_DEVICE static T distance(const V* x, int dtype)
    {
        return point_triangle_distance(
            x[0], x[1], x[2], x[3], PointTriangleDistanceType(dtype));
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Vector<T, 12>
    gradient(const V* x, int dtype)
    {
        return point_triangle_distance_gradient(
            x[0], x[1], x[2], x[3], PointTriangleDistanceType(dtype));
    }
    IPC_TOOLKIT_HOST_DEVICE static Eigen::Matrix<T, 12, 12>
    hessian(const V* x, int dtype)
    {
        return point_triangle_distance_hessian(
            x[0], x[1], x[2], x[3], PointTriangleDistanceType(dtype));
    }
};

// -- What is evaluated --------------------------------------------------------

struct Params {
    double dhat_sqr; ///< The potential is a function of squared distance.
    double kappa;
};

enum class Quantity : uint8_t { VALUE, GRADIENT, HESSIAN, HESSIAN_SUM };

inline const char* quantity_name(const Quantity q)
{
    switch (q) {
    case Quantity::VALUE:
        return "value";
    case Quantity::GRADIENT:
        return "gradient";
    case Quantity::HESSIAN:
        return "hessian";
    default:
        return "hessian_sum";
    }
}

constexpr Quantity QUANTITIES[] = { Quantity::VALUE, Quantity::GRADIENT,
                                    Quantity::HESSIAN, Quantity::HESSIAN_SUM };

// -- The CUDA variants --------------------------------------------------------

#ifdef IPC_TOOLKIT_WITH_CUDA

/// @brief Collisions per packed block for the CUDA variants: a warp.
///
/// The packed layout stores each of the 3·NV coordinates contiguously across
/// the block's lanes, so with 32 lanes the 32 threads of a warp read 32
/// consecutive `R`s per coordinate -- one coalesced transaction, the GPU's
/// analogue of the aligned batch load an `xsimd` lane count buys on the CPU.
constexpr int CUDA_LANES = 32;

/// @brief Whether a CUDA device is present (kernels compile either way).
bool cuda_device_available();

/// @brief The name of CUDA device 0, or "" if there is none.
std::string cuda_device_name();

/// @brief Create the CUDA context up front.
///
/// Otherwise the first thing to touch the device pays ~70 ms for it, which
/// would land on whichever variant happens to be packed first.
void cuda_initialize();

/// @brief One packed group resident on the device, and the kernels over it.
///
/// The device is handed the very same array-of-structures-of-arrays buffers
/// the CPU variants evaluate, only packed at `CUDA_LANES` instead of the
/// architecture's SIMD width, so the two run the same arithmetic on the same
/// numbers in the same order. One thread evaluates one collision -- padding
/// lanes included, as on the CPU, so the work is identical.
template <typename R> class CudaGroup {
public:
    /// @param kind Which stencil the group's collisions are.
    /// @param dtype The group's (single) distance type.
    /// @param n Real collisions; the buffers hold `ceil(n / CUDA_LANES)`
    ///          blocks, whose trailing lanes carry zero weight.
    /// @param x Packed coordinates, `nblocks * 3*NV * CUDA_LANES` of them.
    /// @param w Packed weights, `nblocks * CUDA_LANES` of them.
    CudaGroup(Kind kind, int dtype, size_t n, const R* x, const R* w);

    ~CudaGroup();
    CudaGroup(CudaGroup&&) noexcept;
    CudaGroup& operator=(CudaGroup&&) noexcept;
    CudaGroup(const CudaGroup&) = delete;
    CudaGroup& operator=(const CudaGroup&) = delete;

    /// @brief Seconds the constructor spent uploading `x` and `w`.
    double upload_seconds() const;

    /// @brief Allocate the device buffer `q` writes to (a no-op for the
    /// reducing quantities, and for a buffer that is already the right size).
    void allocate_outputs(Quantity q);

    /// @brief Evaluate `q` on the device, blocking until it is idle.
    /// @return The reduced total for VALUE and HESSIAN_SUM, else 0.
    double run(const Params& p, Quantity q);

    /// @brief Copy the last GRADIENT or HESSIAN output back to the host.
    /// @param out A host buffer in the same packed layout as the input.
    /// @return Seconds spent in the copy.
    double download(Quantity q, R* out) const;

private:
    // Pimpl, so this header stays free of CUDA and Thrust types and can be
    // included by the host translation unit.
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

extern template class CudaGroup<float>;
extern template class CudaGroup<double>;

#endif

} // namespace ipc::tests::bench
