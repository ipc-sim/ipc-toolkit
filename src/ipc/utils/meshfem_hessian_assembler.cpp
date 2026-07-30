#include "meshfem_hessian_assembler.hpp"

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE

#include <ipc/utils/logger.hpp>
#include <ipc/utils/meshfem_eigen_compat.hpp> // must precede MeshFEM headers
#include <ipc/utils/profiler.hpp>

#include <MeshFEMSparse/SystemAssembler.hh>
#include <tbb/parallel_for.h>

#include <cassert>
#include <cstddef>
#include <numeric>
#include <vector>

namespace ipc {

struct MeshFEMHessianAssembler::ImplBase {
    virtual ~ImplBase() = default;
    virtual void
    add(const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) = 0;
    virtual Eigen::SparseMatrix<double> to_eigen() const = 0;
};

/// Dimension-specific implementation (dim ∈ {2, 3}), mirroring MeshFEM's own
/// IPC integration (IPCWrapper.cc): one dim-sized block variable per vertex.
template <int dim>
struct MeshFEMHessianAssembler::Impl final
    : public MeshFEMHessianAssembler::ImplBase {
    using VarStructure = MeshFEM::OptimizationVarStructure<dim>;
    /// Local (block) variables of a collision stencil: 1–4 vertices.
    using Stencil = MeshFEM::ElementBlockVarsWithSizeRange<1, 4>;

    Impl(
        const size_t num_block_vars,
        const size_t num_stencils,
        const StencilGetter& stencil)
        : m_assembler(num_block_vars)
        , m_vars(num_block_vars)
    {
        {
            IPC_TOOLKIT_PROFILE_BLOCK("MeshFEM block sparsity pattern");
            m_H = m_assembler.blockSparsityPattern(
                num_stencils,
                [&stencil](const size_t i) { return to_stencil(stencil(i)); });
        }
        m_H->setZero(); // Allocate (and zero) the value array.
        m_locks.init(num_block_vars);
    }

    void
    add(const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) override
    {
        assert(local_hess.rows() == local_hess.cols());
        assert(local_hess.rows() % dim == 0);

        const Stencil evars = to_stencil(vertex_ids);
        assert(evars.size() == size_t(local_hess.rows()) / dim);

        // Per-vertex-pair dim×dim block of the local Hessian; (a, b) are
        // scalar offsets computed by the assembler (multiples of dim).
        const auto local_block = [&local_hess](
                                     const size_t a, const size_t b,
                                     const size_t /*block_rows*/,
                                     const size_t /*block_cols*/) {
            return local_hess.template block<dim, dim>(a, b);
        };

        // Sorted column-merge scatter into the block-CSC value array, with
        // per-block-column spin locks for thread safety.
        MeshFEM::ElementHessianContribAssembler<
            /*UseBlockMergeAlgorithm=*/true>::
            template run</*InParallel=*/true>(
                m_H->Ax.data(), *m_H, local_block, evars, m_vars, m_locks);
    }

    // Convert the block-CSC upper-triangle matrix to a full symmetric Eigen
    // matrix.
    //
    // NOTE: We deliberately do not use BlockCSCHessian::toEigen here: its
    // toScalar step assumes every block column is non-empty (true for FE
    // Hessians, where every node belongs to an element, but false for contact
    // Hessians, where most vertices are collision-free) and reads out of
    // bounds otherwise. This implementation also skips toEigen's intermediate
    // upper-triangle scalar matrix, symmetrizing directly instead.
    //
    // Value layout (ContiguousBlocks=true, StoreFullDiagonalBlocks=true, the
    // library default): block entry ii occupies Ax[N²·ii, N²·(ii+1)),
    // column-major within the block; diagonal blocks are full N×N with a
    // zeroed strict lower triangle.
    Eigen::SparseMatrix<double> to_eigen() const override
    {
        IPC_TOOLKIT_PROFILE_BLOCK("MeshFEM block CSC to Eigen");
        static constexpr int N = dim;

        const auto& Ap = m_H->Ap; // block column pointers (size n_blocks + 1)
        const auto& Ai = m_H->Ai; // block row indices
        const auto& Ax = m_H->Ax; // scalar values (N² per block)

        const std::ptrdiff_t n_blocks = m_H->n;
        const std::ptrdiff_t n = N * n_blocks; // scalar dimension

        // --- Symmetrized *block* pattern -------------------------------
        // Every stored block (bi, bj) (bi ≤ bj) appears at (bi, bj) and, if
        // off-diagonal, mirrored at (bj, bi) in the full matrix.
        struct BlockEntry {
            std::ptrdiff_t row; ///< Block row in the full (symmetrized) matrix
            std::ptrdiff_t src; ///< Index of the stored block (into Ai/Ax)
            bool transposed;    ///< Whether the stored block is mirrored
        };

        std::vector<std::ptrdiff_t> sym_col_start(n_blocks + 1, 0);
        for (std::ptrdiff_t bj = 0; bj < n_blocks; bj++) {
            for (std::ptrdiff_t ii = Ap[bj]; ii < Ap[bj + 1]; ii++) {
                sym_col_start[bj + 1]++;
                if (Ai[ii] != bj) {
                    sym_col_start[Ai[ii] + 1]++;
                }
            }
        }
        std::partial_sum(
            sym_col_start.begin(), sym_col_start.end(), sym_col_start.begin());
        const std::ptrdiff_t num_sym_blocks = sym_col_start[n_blocks];

        // Sweeping bj in ascending order appends rows to every column in
        // ascending order: column c receives its upper entries (rows ≤ c) at
        // step bj = c and its mirrored entries (rows bj > c) at later steps.
        std::vector<BlockEntry> sym_entries(num_sym_blocks);
        {
            std::vector<std::ptrdiff_t> cursor(
                sym_col_start.begin(), sym_col_start.end() - 1);
            for (std::ptrdiff_t bj = 0; bj < n_blocks; bj++) {
                for (std::ptrdiff_t ii = Ap[bj]; ii < Ap[bj + 1]; ii++) {
                    const std::ptrdiff_t bi = Ai[ii];
                    sym_entries[cursor[bj]++] = { bi, ii, false };
                    if (bi != bj) {
                        sym_entries[cursor[bi]++] = { bj, ii, true };
                    }
                }
            }
        }

        // --- Expand to a scalar CSC matrix ------------------------------
        Eigen::SparseMatrix<double> M(n, n);
        M.makeCompressed();
        M.resizeNonZeros(N * N * num_sym_blocks);

        int* outer = M.outerIndexPtr();
        int* inner = M.innerIndexPtr();
        double* values = M.valuePtr();

        for (std::ptrdiff_t bj = 0; bj <= n_blocks; bj++) {
            const std::ptrdiff_t col_nnz = bj < n_blocks
                ? N * (sym_col_start[bj + 1] - sym_col_start[bj])
                : 0;
            for (int c = 0; c < ((bj < n_blocks) ? N : 1); c++) {
                outer[N * bj + c] =
                    int(N * N * sym_col_start[bj] + c * col_nnz);
            }
        }

        tbb::parallel_for(
            std::ptrdiff_t(0), n_blocks, [&](const std::ptrdiff_t bj) {
                for (int c = 0; c < N; c++) {
                    std::ptrdiff_t out = outer[N * bj + c];
                    for (std::ptrdiff_t ei = sym_col_start[bj];
                         ei < sym_col_start[bj + 1]; ei++) {
                        const BlockEntry& e = sym_entries[ei];
                        const double* block = Ax.data() + N * N * e.src;
                        const bool diagonal = e.row == bj;
                        for (int r = 0; r < N; r++) {
                            inner[out] = int(N * e.row + r);
                            if (e.transposed || (diagonal && r > c)) {
                                values[out] = block[r * N + c]; // transposed
                            } else {
                                values[out] = block[c * N + r];
                            }
                            out++;
                        }
                    }
                }
            });

        return M;
    }

private:
    /// Compact the (possibly -1-padded) vertex ID array into a MeshFEM
    /// stencil of block variables.
    static Stencil to_stencil(const std::array<index_t, 4>& vertex_ids)
    {
        Stencil bvars(0);
        size_t back = 0;
        for (const index_t id : vertex_ids) {
            if (id >= 0) {
                bvars[back++] = id;
            }
        }
        bvars.resize(back);
        return bvars;
    }

    MeshFEM::SystemAssembler<dim> m_assembler;
    VarStructure m_vars;
    std::unique_ptr<MeshFEM::BlockCSCHessian<VarStructure>> m_H;
    MeshFEM::VarLocks m_locks;
};

MeshFEMHessianAssembler::MeshFEMHessianAssembler() = default;
MeshFEMHessianAssembler::~MeshFEMHessianAssembler() = default;

void MeshFEMHessianAssembler::begin(
    const int ndof,
    const int dim,
    const size_t num_stencils,
    const StencilGetter& stencil)
{
    assert(ndof % dim == 0);
    const size_t num_block_vars = size_t(ndof) / dim;
    if (dim == 2) {
        m_impl =
            std::make_unique<Impl<2>>(num_block_vars, num_stencils, stencil);
    } else if (dim == 3) {
        m_impl =
            std::make_unique<Impl<3>>(num_block_vars, num_stencils, stencil);
    } else {
        log_and_throw_error(
            "MeshFEMHessianAssembler: unsupported dimension {}!", dim);
    }
}

void MeshFEMHessianAssembler::add_local_hessian(
    const MatrixMax12d& local_hess, const std::array<index_t, 4>& vertex_ids)
{
    assert(m_impl != nullptr);
    m_impl->add(local_hess, vertex_ids);
}

Eigen::SparseMatrix<double> MeshFEMHessianAssembler::get_matrix() const
{
    assert(m_impl != nullptr);
    return m_impl->to_eigen();
}

} // namespace ipc

#endif // IPC_TOOLKIT_WITH_MESHFEM_SPARSE
