#pragma once

#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/profiler.hpp>

#include <Eigen/Core>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <array>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ipc {

namespace detail {

    template <typename T> struct IsStdArray : std::false_type { };

    template <typename T, size_t N>
    struct IsStdArray<std::array<T, N>> : std::true_type { };

} // namespace detail

/// @brief Whether an assembly of this shape uses the two-pass strategy.
///
/// Exposed so that callers' tests can assert which strategy a case exercises
/// without restating the condition.
///
/// @param num_slots Total stencil-vertex contributions (stencil size × count).
/// @param out_ndof Size of the assembled gradient.
inline bool
gradient_assembly_is_two_pass(const size_t num_slots, const size_t out_ndof)
{
    return out_ndof > num_slots;
}

/// @brief Assemble a global gradient from per-stencil local gradients.
///
/// Evaluates the local gradients in parallel and sums them into a dense global
/// vector. Two strategies with complementary cost scaling are available, and
/// the crossover between them sits inside the operating range of real contact
/// scenes, so the strategy is chosen per call from the problem shape (see
/// gradient_assembly_is_two_pass). Measured on the sweep in
/// tests/src/tests/potential/benchmark_gradient_assembly.cpp:
///
///  - Two-pass costs num_slots adds in the serial scatter, plus zeroing the
///    output once, so it scales with the contact set. It wins when contact is
///    sparse relative to the DOF count (the typical case), by up to ~75x, and
///    is bitwise deterministic because the scatter order is fixed. Requires
///    bounded stencils (see @note).
///  - Thread-local accumulators cost out_ndof·nthreads to zero and reduce, so
///    they scale with the DOF count. They win in the dense self-contact
///    regime, by up to ~2x.
///
/// @note Two-pass applies only to bounded stencils: it buffers one local
/// gradient and one id set per stencil, and a fixed-max Eigen vector stores
/// them inline, so the buffer costs a single allocation. Callers whose local
/// gradients are dynamically sized (e.g. smooth contact, whose stencils are not
/// bounded at four vertices) always take the thread-local path.
///
/// @param out_ndof Size of the assembled gradient.
/// @param dim Spatial dimension (entries per vertex).
/// @param num_stencils Number of local gradients to assemble.
/// @param local_gradient Callable: stencil index → its local gradient. Must be
///     safe to call concurrently.
/// @param stencil_ids Callable: stencil index → the global vertex ids of that
///     stencil, already mapped to the output's indexing. Must be safe to call
///     concurrently.
/// @returns The assembled gradient, of size out_ndof.
template <typename LocalGradientFn, typename StencilIDsFn>
Eigen::VectorXd assemble_gradient(
    const int out_ndof,
    const int dim,
    const size_t num_stencils,
    LocalGradientFn&& local_gradient,
    StencilIDsFn&& stencil_ids)
{
    using LocalGrad =
        std::decay_t<std::invoke_result_t<LocalGradientFn&, size_t>>;
    using IDs = std::decay_t<std::invoke_result_t<StencilIDsFn&, size_t>>;

    constexpr bool CAN_TWO_PASS =
        LocalGrad::MaxSizeAtCompileTime != Eigen::Dynamic
        && detail::IsStdArray<IDs>::value;

    if (num_stencils == 0) {
        return Eigen::VectorXd::Zero(out_ndof);
    }

    if constexpr (CAN_TWO_PASS) {
        constexpr size_t STENCIL_SIZE = std::tuple_size_v<IDs>;

        if (gradient_assembly_is_two_pass(
                STENCIL_SIZE * num_stencils, size_t(out_ndof))) {
            std::vector<LocalGrad> local_grads(num_stencils);
            std::vector<IDs> local_ids(num_stencils);

            {
                IPC_TOOLKIT_PROFILE_BLOCK("Compute Local Gradients");
                tbb::parallel_for(size_t(0), num_stencils, [&](const size_t i) {
                    local_grads[i] = local_gradient(i);
                    local_ids[i] = stencil_ids(i);
                });
            }

            IPC_TOOLKIT_PROFILE_BLOCK("Assemble Gradient");
            Eigen::VectorXd grad = Eigen::VectorXd::Zero(out_ndof);
            for (size_t i = 0; i < num_stencils; i++) {
                local_gradient_to_global_gradient(
                    local_grads[i], local_ids[i], dim, grad);
            }
            return grad;
        }
    }

    tbb::enumerable_thread_specific<Eigen::VectorXd> partials([out_ndof] {
        return Eigen::VectorXd(Eigen::VectorXd::Zero(out_ndof));
    });

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Compute Local Gradients");
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), num_stencils),
            [&](const tbb::blocked_range<size_t>& r) {
                Eigen::VectorXd& local_sum = partials.local();
                for (size_t i = r.begin(); i < r.end(); i++) {
                    local_gradient_to_global_gradient(
                        local_gradient(i), stencil_ids(i), dim, local_sum);
                }
            });
    }

    IPC_TOOLKIT_PROFILE_BLOCK("Assemble Gradient");
    std::vector<const Eigen::VectorXd*> ptrs;
    for (const Eigen::VectorXd& p : partials) {
        ptrs.push_back(&p);
    }
    Eigen::VectorXd combined(out_ndof);
    tbb::parallel_for(
        tbb::blocked_range<int>(0, out_ndof),
        [&](const tbb::blocked_range<int>& r) {
            auto seg = combined.segment(r.begin(), r.size());
            seg.setZero();
            for (const Eigen::VectorXd* p : ptrs) {
                seg += p->segment(r.begin(), r.size());
            }
        });
    return combined;
}

} // namespace ipc
