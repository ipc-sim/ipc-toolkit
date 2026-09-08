// The CUDA half of the per-scalar-type barrier potential benchmark: the same
// barrier chain the host variants evaluate, one thread per collision, reading
// the very same packed buffers.
//
//   d = dist²(x),  b = κ·w·f(d),  ∇b = κ·w·f'(d)·∇d,
//   ∇²b = κ·w·(f"(d)·∇d∇dᵀ + f'(d)·∇²d)
//
// The distance type is resolved host-side and passed as a launch argument, so
// a whole grid shares one instantiation and no thread branches on it. The
// edge-edge mollifier is left out here as it is on the host, keeping the
// variants on identical arithmetic.
//
// Only `float` and `double` appear; `ipc::detail`'s distance templates and the
// barrier trio are explicitly instantiated as device code for exactly those
// two (see ipc_toolkit_shared_device_sources.cmake).

#include "barrier_potential_bench.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/cuda/device_utils.cuh>

#include <thrust/device_vector.h>
#include <thrust/fill.h>

#include <chrono>
#include <cstdlib>
#include <type_traits>

namespace ipc::tests::bench {

namespace {

    /// @brief Threads per block.
    ///
    /// Deliberately not `ipc::cuda::KERNEL_BLOCK_SIZE` (256). Each group is
    /// its own launch and the groups are small next to the machine, so what
    /// binds is how many SMs a launch reaches, not steady-state occupancy:
    /// cloth-ball's largest group is 2592 padded threads, which is 11 blocks
    /// at 256 -- around a tenth of the SMs -- against 41 at 64. Measured on an
    /// RTX 5080, 64 is 1.1x to 2.0x faster than 256 across the scenes, the
    /// margin shrinking as a scene grows big enough to fill the GPU anyway.
    ///
    /// 512 and beyond cannot run at all: the edge-vertex Hessian kernel needs
    /// 250 registers a thread, and 250 x 512 exceeds an SM's register file.
    ///
    /// Must be a multiple of `CUDA_LANES` -- the reductions assume a block
    /// holds whole warps. `IPC_TOOLKIT_BENCH_CUDA_BLOCK` overrides it, so the
    /// launch configuration can be re-swept without a rebuild.
    int block_size()
    {
        constexpr int DEFAULT_BLOCK_SIZE = 64;
        static const int size = [] {
            const char* s = std::getenv("IPC_TOOLKIT_BENCH_CUDA_BLOCK");
            const int v = s != nullptr ? std::atoi(s) : 0;
            return (v >= CUDA_LANES && v <= 1024 && v % CUDA_LANES == 0)
                ? v
                : DEFAULT_BLOCK_SIZE;
        }();
        return size;
    }

    /// @brief Seconds since an arbitrary epoch, for the host-side timers.
    double now_seconds()
    {
        return std::chrono::duration<double>(
                   std::chrono::steady_clock::now().time_since_epoch())
            .count();
    }

    /// @brief Gather one thread's collision out of the packed layout.
    ///
    /// Thread `i` owns lane `i % L` of block `i / L`, whose coordinate `c`
    /// lives at `x[(i / L) * N * L + c * L + (i % L)]`. Consecutive threads
    /// differ only in the lane, so each of the `N` reads below is contiguous
    /// across the warp.
    template <typename R, int N, int NV>
    __device__ void load(
        const R* x,
        const R* w,
        const size_t i,
        Eigen::Vector3<R>* out,
        R& weight)
    {
        constexpr int L = CUDA_LANES;
        const R* block = x + (i / L) * (N * L) + (i % L);
        for (int v = 0; v < NV; v++) {
            for (int k = 0; k < 3; k++) {
                out[v][k] = block[(3 * v + k) * L];
            }
        }
        weight = w[(i / L) * L + (i % L)];
    }

    /// @brief Sum `v` across the warp and add one result per warp to `*total`.
    ///
    /// Summing in `R` first bounds a float's accumulation error to 32 terms --
    /// tighter than the host loop's 256-collision chunk -- and leaves one
    /// atomic per 32 collisions rather than one per collision. The order the
    /// warps land in differs from the host's, so the two totals agree only to
    /// rounding, which is what the accuracy column reports.
    template <typename R> __device__ void warp_reduce_add(R v, double* total)
    {
        static_assert(CUDA_LANES == 32, "A block's lanes are one warp.");
#pragma unroll
        for (int offset = CUDA_LANES / 2; offset > 0; offset >>= 1) {
            v += __shfl_down_sync(0xffffffff, v, offset);
        }
        if ((threadIdx.x % CUDA_LANES) == 0) {
            atomicAdd(total, double(v));
        }
    }

    template <typename R, Kind K> struct Kernels {
        using S = Stencil<R, K>;
        static constexpr int NV = S::NV;
        static constexpr int N = 3 * NV;

        /// @brief κ·w·(f"(d)·∇d∇dᵀ + f'(d)·∇²d), the host's `local_hessian`.
        __device__ static Eigen::Matrix<R, N, N> local_hessian(
            const Eigen::Vector3<R>* x,
            const R w,
            const int dtype,
            const R dhat,
            const R kappa)
        {
            const R d = S::distance(x, dtype);
            const Eigen::Vector<R, N> grad_d = S::gradient(x, dtype);
            const Eigen::Matrix<R, N, N> hess_d = S::hessian(x, dtype);
            const R grad_f = barrier_first_derivative(d, dhat);
            const R hess_f = barrier_second_derivative(d, dhat);
            const R kw = kappa * w;
            return (kw * hess_f) * (grad_d * grad_d.transpose())
                + (kw * grad_f) * hess_d;
        }
    };

    // -- The kernels ----------------------------------------------------------
    //
    // Every launch covers `nthreads = nblocks * CUDA_LANES` collisions, the
    // padded count the host variants also evaluate. A thread past the end
    // contributes zero rather than returning early, so the warp-wide shuffles
    // in warp_reduce_add() always see a converged warp.

    template <typename R, Kind K>
    __global__ void value_kernel(
        const R* x,
        const R* w,
        const int dtype,
        const R dhat,
        const R kappa,
        const size_t nthreads,
        double* total)
    {
        using Kn = Kernels<R, K>;
        const size_t i =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

        R v = R(0);
        if (i < nthreads) {
            Eigen::Vector3<R> xi[Kn::NV];
            R wi;
            load<R, Kn::N, Kn::NV>(x, w, i, xi, wi);
            const R d = Stencil<R, K>::distance(xi, dtype);
            v = (kappa * wi) * barrier(d, dhat);
        }
        warp_reduce_add(v, total);
    }

    template <typename R, Kind K>
    __global__ void gradient_kernel(
        const R* x,
        const R* w,
        const int dtype,
        const R dhat,
        const R kappa,
        const size_t nthreads,
        R* out)
    {
        using Kn = Kernels<R, K>;
        constexpr int L = CUDA_LANES;
        const size_t i =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (i >= nthreads) {
            return;
        }

        Eigen::Vector3<R> xi[Kn::NV];
        R wi;
        load<R, Kn::N, Kn::NV>(x, w, i, xi, wi);

        const R d = Stencil<R, K>::distance(xi, dtype);
        const Eigen::Vector<R, Kn::N> grad_d =
            Stencil<R, K>::gradient(xi, dtype);
        const R grad_f = barrier_first_derivative(d, dhat);
        const Eigen::Vector<R, Kn::N> grad = (kappa * wi * grad_f) * grad_d;

        R* dst = out + (i / L) * (Kn::N * L) + (i % L);
        for (int k = 0; k < Kn::N; k++) {
            dst[k * L] = grad[k];
        }
    }

    template <typename R, Kind K>
    __global__ void hessian_kernel(
        const R* x,
        const R* w,
        const int dtype,
        const R dhat,
        const R kappa,
        const size_t nthreads,
        R* out)
    {
        using Kn = Kernels<R, K>;
        constexpr int L = CUDA_LANES;
        const size_t i =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (i >= nthreads) {
            return;
        }

        Eigen::Vector3<R> xi[Kn::NV];
        R wi;
        load<R, Kn::N, Kn::NV>(x, w, i, xi, wi);
        const Eigen::Matrix<R, Kn::N, Kn::N> hess =
            Kn::local_hessian(xi, wi, dtype, dhat, kappa);

        R* dst = out + (i / L) * (Kn::N * Kn::N * L) + (i % L);
        for (int k = 0; k < Kn::N * Kn::N; k++) {
            dst[k * L] = hess.data()[k];
        }
    }

    /// @brief The Hessian's compute cost alone: the N² entries are summed into
    /// one accumulator instead of being written out, so the stores -- which
    /// dominate once the bus is saturated -- are not part of the measurement.
    template <typename R, Kind K>
    __global__ void hessian_sum_kernel(
        const R* x,
        const R* w,
        const int dtype,
        const R dhat,
        const R kappa,
        const size_t nthreads,
        double* total)
    {
        using Kn = Kernels<R, K>;
        const size_t i =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

        R v = R(0);
        if (i < nthreads) {
            Eigen::Vector3<R> xi[Kn::NV];
            R wi;
            load<R, Kn::N, Kn::NV>(x, w, i, xi, wi);
            v = Kn::local_hessian(xi, wi, dtype, dhat, kappa).sum();
        }
        warp_reduce_add(v, total);
    }

} // namespace

// ============================================================================

bool cuda_device_available()
{
    int device_count = 0;
    const cudaError_t err = cudaGetDeviceCount(&device_count);
    return err == cudaSuccess && device_count > 0;
}

std::string cuda_device_name()
{
    if (!cuda_device_available()) {
        return "";
    }
    cudaDeviceProp props {};
    if (cudaGetDeviceProperties(&props, 0) != cudaSuccess) {
        return "";
    }
    return props.name;
}

void cuda_initialize()
{
    // Any allocation forces the lazy context creation.
    IPC_TOOLKIT_CUDA_CHECK(cudaFree(nullptr));
    IPC_TOOLKIT_CUDA_CHECK(cudaDeviceSynchronize());
}

// ----------------------------------------------------------------------------
// CudaGroup

template <typename R> struct CudaGroup<R>::Impl {
    Kind kind = Kind::VV;
    int dtype = 0;
    size_t nblocks = 0;
    double upload_s = 0;

    thrust::device_vector<R> x;
    thrust::device_vector<R> w;
    thrust::device_vector<R> out;        ///< Gradient or Hessian entries.
    thrust::device_vector<double> total; ///< One scalar, for the reductions.

    int ndof() const { return 3 * num_vertices(kind); }
    /// @brief Collisions launched, i.e. the padded count.
    size_t nthreads() const { return nblocks * size_t(CUDA_LANES); }
};

template <typename R>
CudaGroup<R>::CudaGroup(
    const Kind kind, const int dtype, const size_t n, const R* x, const R* w)
    : m_impl(std::make_unique<Impl>())
{
    Impl& impl = *m_impl;
    impl.kind = kind;
    impl.dtype = dtype;
    impl.nblocks = (n + CUDA_LANES - 1) / CUDA_LANES;
    impl.total.resize(1);

    const size_t nx = impl.nthreads() * size_t(impl.ndof());
    const size_t nw = impl.nthreads();

    // Allocation is not part of the upload time: a simulator would allocate
    // once and re-upload per iteration.
    impl.x.resize(nx);
    impl.w.resize(nw);
    IPC_TOOLKIT_CUDA_CHECK(cudaDeviceSynchronize());

    const double start = now_seconds();
    IPC_TOOLKIT_CUDA_CHECK(cudaMemcpy(
        thrust::raw_pointer_cast(impl.x.data()), x, nx * sizeof(R),
        cudaMemcpyHostToDevice));
    IPC_TOOLKIT_CUDA_CHECK(cudaMemcpy(
        thrust::raw_pointer_cast(impl.w.data()), w, nw * sizeof(R),
        cudaMemcpyHostToDevice));
    impl.upload_s = now_seconds() - start;
}

template <typename R> CudaGroup<R>::~CudaGroup() = default;
template <typename R> CudaGroup<R>::CudaGroup(CudaGroup&&) noexcept = default;
template <typename R>
CudaGroup<R>& CudaGroup<R>::operator=(CudaGroup&&) noexcept = default;

template <typename R> double CudaGroup<R>::upload_seconds() const
{
    return m_impl->upload_s;
}

template <typename R> void CudaGroup<R>::allocate_outputs(const Quantity q)
{
    Impl& impl = *m_impl;
    size_t needed = 0;
    if (q == Quantity::GRADIENT) {
        needed = impl.nthreads() * size_t(impl.ndof());
    } else if (q == Quantity::HESSIAN) {
        needed = impl.nthreads() * size_t(impl.ndof()) * size_t(impl.ndof());
    }
    if (impl.out.size() == needed) {
        return;
    }
    impl.out.resize(needed);
    // Touch the buffer so a first write is not what pages it in, matching the
    // host variants' warm-up.
    thrust::fill(impl.out.begin(), impl.out.end(), R(0));
    IPC_TOOLKIT_CUDA_CHECK(cudaDeviceSynchronize());
}

template <typename R>
double CudaGroup<R>::run(const Params& p, const Quantity q)
{
    Impl& impl = *m_impl;
    const R dhat = R(p.dhat_sqr);
    const R kappa = R(p.kappa);
    const size_t nthreads = impl.nthreads();
    const int BLOCK = block_size();
    const int grid = int((nthreads + BLOCK - 1) / BLOCK);

    const R* x = thrust::raw_pointer_cast(impl.x.data());
    const R* w = thrust::raw_pointer_cast(impl.w.data());
    R* out = thrust::raw_pointer_cast(impl.out.data());
    double* total = thrust::raw_pointer_cast(impl.total.data());

    const bool reduces = q == Quantity::VALUE || q == Quantity::HESSIAN_SUM;
    if (reduces) {
        IPC_TOOLKIT_CUDA_CHECK(cudaMemset(total, 0, sizeof(double)));
    }

    // The kind is a runtime property of the group but a compile-time one of
    // the kernel, so it is dispatched here, once per launch, not per thread.
    const auto launch = [&](auto kind_tag) {
        constexpr Kind K = decltype(kind_tag)::value;
        switch (q) {
        case Quantity::VALUE:
            value_kernel<R, K><<<grid, BLOCK>>>(
                x, w, impl.dtype, dhat, kappa, nthreads, total);
            break;
        case Quantity::GRADIENT:
            gradient_kernel<R, K>
                <<<grid, BLOCK>>>(x, w, impl.dtype, dhat, kappa, nthreads, out);
            break;
        case Quantity::HESSIAN:
            hessian_kernel<R, K>
                <<<grid, BLOCK>>>(x, w, impl.dtype, dhat, kappa, nthreads, out);
            break;
        default:
            hessian_sum_kernel<R, K><<<grid, BLOCK>>>(
                x, w, impl.dtype, dhat, kappa, nthreads, total);
            break;
        }
    };
    switch (impl.kind) {
    case Kind::VV:
        launch(std::integral_constant<Kind, Kind::VV>());
        break;
    case Kind::EV:
        launch(std::integral_constant<Kind, Kind::EV>());
        break;
    case Kind::EE:
        launch(std::integral_constant<Kind, Kind::EE>());
        break;
    default:
        launch(std::integral_constant<Kind, Kind::FV>());
        break;
    }
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    if (!reduces) {
        IPC_TOOLKIT_CUDA_CHECK(cudaDeviceSynchronize());
        return 0.0;
    }
    // Reading the accumulator back is a synchronizing 8-byte copy, so this
    // also serves as the barrier the other quantities take explicitly.
    double host_total = 0;
    IPC_TOOLKIT_CUDA_CHECK(
        cudaMemcpy(&host_total, total, sizeof(double), cudaMemcpyDeviceToHost));
    return host_total;
}

template <typename R>
double CudaGroup<R>::download(const Quantity q, R* out) const
{
    if (q != Quantity::GRADIENT && q != Quantity::HESSIAN) {
        return 0.0;
    }
    const Impl& impl = *m_impl;
    const double start = now_seconds();
    IPC_TOOLKIT_CUDA_CHECK(cudaMemcpy(
        out, thrust::raw_pointer_cast(impl.out.data()),
        impl.out.size() * sizeof(R), cudaMemcpyDeviceToHost));
    return now_seconds() - start;
}

template class CudaGroup<float>;
template class CudaGroup<double>;

} // namespace ipc::tests::bench

#endif // IPC_TOOLKIT_WITH_CUDA
