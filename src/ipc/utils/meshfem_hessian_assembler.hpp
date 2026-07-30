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
/// The assembled matrix is symmetric and stored upper-triangle-only in block
/// CSC format; get_matrix() converts to a full symmetric Eigen matrix. This
/// conversion costs a full copy — direct (zero-copy) access to the block-CSC
/// format is planned as part of the persistent-pattern API (Phase 4).
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
    Eigen::SparseMatrix<double> get_matrix() const;

private:
    struct ImplBase;
    template <int dim> struct Impl;
    std::unique_ptr<ImplBase> m_impl;
};

} // namespace ipc

#endif // IPC_TOOLKIT_WITH_MESHFEM_SPARSE
