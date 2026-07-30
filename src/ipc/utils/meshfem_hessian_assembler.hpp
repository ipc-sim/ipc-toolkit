#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE

#include <ipc/utils/hessian_assembler.hpp>

#include <memory>

namespace ipc {

// The MeshFEM backend scatters per-vertex d×d blocks, which requires
// derivatives to be laid out [x0, y0, z0, x1, ...].
static_assert(
    VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor,
    "The MeshFEMSparse assembly backend requires "
    "IPC_TOOLKIT_VERTEX_DERIVATIVE_LAYOUT=RowMajor.");

/// @brief HessianAssembler backed by MeshFEMSparse's block-CSC data
/// structures (Mohammadian et al. 2026).
///
/// begin() builds a block sparsity pattern (one dim×dim block per interacting
/// vertex pair) from the stencils; add_local_hessian() scatters each local
/// Hessian directly into the pattern's value array using a sorted
/// column-merge with per-column locks — no triplets, no setFromTriplets.
///
/// The assembler is designed to be **reused across assemblies** (e.g., one
/// instance per Newton solve): begin() compares the stencils against the
/// cached sparsity pattern and rebuilds it only if the contact set gained new
/// entries or lost more than stale_block_tolerance() blocks; otherwise the
/// pattern (and the cached Eigen structure) are reused and only the values
/// are recomputed. Stale blocks left in a reused pattern assemble to explicit
/// zeros, which do not affect the matrix's value.
///
/// The assembled matrix is symmetric and stored upper-triangle-only in block
/// CSC format; get_matrix() converts to a full symmetric Eigen matrix,
/// reusing the cached structure when the pattern is unchanged.
class MeshFEMHessianAssembler final : public HessianAssembler {
public:
    MeshFEMHessianAssembler();
    ~MeshFEMHessianAssembler() override;

    void
    begin(int ndof, int dim, size_t num_stencils, const StencilGetter& stencil)
        override;

    void add_local_hessian(
        const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) override;

    void end() override { } // Values are scattered in place; nothing to merge.

    /// @brief Convert the assembled matrix to a full symmetric Eigen matrix.
    ///
    /// The reference stays valid (and its values current) until the next
    /// begin() call; copy it to keep a snapshot.
    const Eigen::SparseMatrix<double>& get_matrix() const;

    /// @brief Number of vanished blocks tolerated before a pattern rebuild.
    ///
    /// When a stencil introduces a vertex pair absent from the cached pattern,
    /// the pattern is always rebuilt. When blocks merely disappear (e.g., a
    /// contact separates), the pattern is reused as long as at most this many
    /// blocks vanished. Mirrors MeshFEM's sparsityPatternUpdateThreshold.
    size_t stale_block_tolerance() const { return m_stale_block_tolerance; }

    /// @brief Set the number of vanished blocks tolerated before a rebuild.
    void set_stale_block_tolerance(const size_t tolerance)
    {
        m_stale_block_tolerance = tolerance;
    }

    /// @brief Whether the last begin() call reused the cached pattern.
    bool reused_pattern() const { return m_reused_pattern; }

    /// @brief Whether begin() may skip change detection entirely.
    bool assume_unchanged_stencils() const
    {
        return m_assume_unchanged_stencils;
    }

    /// @brief Allow begin() to skip change detection entirely.
    ///
    /// When enabled, begin() reuses the cached pattern without comparing the
    /// stencils against it, as long as a pattern exists and the stencil count
    /// is unchanged (a differing count falls back to normal detection). Use
    /// this when the caller knows the collision set is identical to the
    /// previous assembly (e.g., reassembling with a different PSD projection
    /// or stiffness): on large scenes, change detection costs as much as a
    /// pattern rebuild.
    ///
    /// @warning If the stencils did change (with an equal count), assembly
    /// reads out of bounds. Debug builds verify the assumption.
    void set_assume_unchanged_stencils(const bool assume)
    {
        m_assume_unchanged_stencils = assume;
    }

private:
    struct ImplBase;
    template <int dim> struct Impl;
    std::unique_ptr<ImplBase> m_impl;

    size_t m_stale_block_tolerance = 0;
    bool m_assume_unchanged_stencils = false;
    bool m_reused_pattern = false;
};

} // namespace ipc

#endif // IPC_TOOLKIT_WITH_MESHFEM_SPARSE
