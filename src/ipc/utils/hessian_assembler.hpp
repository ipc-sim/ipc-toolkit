#pragma once

#include <ipc/config.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <tbb/enumerable_thread_specific.h>

#include <array>
#include <functional>
#include <memory>

namespace ipc {

/// @brief Abstract sink for assembling local (per-collision) Hessians into a
/// global matrix.
///
/// The driver (e.g., Potential::hessian) evaluates one local Hessian per
/// collision stencil in parallel and hands each to add_local_hessian(). A
/// concrete assembler decides how the global matrix is stored and built (e.g.,
/// triplets + setFromTriplets, or a persistent block-sparse pattern), and
/// exposes its own accessors for the result.
///
/// Call sequence: begin(), then any number of add_local_hessian() calls
/// (possibly concurrent), then end().
class HessianAssembler {
public:
    /// @brief Callable returning the global vertex IDs of stencil i.
    /// The IDs match those later passed to add_local_hessian for the same
    /// stencil index (invalid trailing entries are negative). Must be safe to
    /// call concurrently.
    using StencilGetter = std::function<std::array<index_t, 4>(size_t)>;

    virtual ~HessianAssembler() = default;

    /// @brief Prepare for assembly. Called once before any add_local_hessian.
    /// @param ndof Number of global scalar DOF (rows == cols of the result).
    /// @param dim Spatial dimension (rows/cols per vertex block).
    /// @param num_stencils Number of local Hessians that will be added.
    /// @param stencil Enumerates the stencils' vertex IDs; pattern-based
    ///     assemblers use this to build their sparsity pattern up front. Only
    ///     valid for the duration of the begin() call.
    virtual void begin(
        int ndof,
        int dim,
        size_t num_stencils,
        const StencilGetter& stencil) = 0;

    /// @brief Add one local (stencil) Hessian to the global matrix.
    ///
    /// Must be safe to call concurrently from multiple threads between
    /// begin() and end().
    ///
    /// @param local_hess Local Hessian of size (n·dim)×(n·dim), where
    ///     n = local_hess.rows() / dim is the number of stencil vertices.
    /// @param vertex_ids Global vertex IDs of the stencil; the first n entries
    ///     are valid (remaining entries may be negative placeholders).
    virtual void add_local_hessian(
        const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) = 0;

    /// @brief Finish assembly. Called once after all add_local_hessian calls.
    virtual void end() = 0;
};

/// @brief The default HessianAssembler: thread-local triplet caches merged
/// into an Eigen::SparseMatrix via setFromTriplets.
///
/// This reproduces the historical behavior of Potential::hessian exactly.
class TripletHessianAssembler final : public HessianAssembler {
public:
    TripletHessianAssembler() = default;

    /// @note `stencil` is unused: the triplet path needs no sparsity pattern.
    void
    begin(int ndof, int dim, size_t num_stencils, const StencilGetter& stencil)
        override;

    void add_local_hessian(
        const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) override;

    void end() override { } // All work happens in get_matrix().

    /// @brief Merge the thread-local caches and build the global matrix.
    /// Call once, after end(); the internal caches are consumed.
    ///
    /// @throws std::runtime_error if called before the first assembly.
    Eigen::SparseMatrix<double> get_matrix();

private:
    int m_ndof = 0;
    int m_dim = 0;
    std::unique_ptr<tbb::enumerable_thread_specific<LocalThreadMatStorage>>
        m_storage;
};

} // namespace ipc
