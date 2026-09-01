#include "meshfem_hessian_assembler.hpp"

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE

#include <ipc/utils/logger.hpp>
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

    virtual int dimension() const = 0;
    virtual size_t num_block_vars() const = 0;

    /// Reuse the cached pattern if the stencils still fit (returns true) or
    /// rebuild it (returns false). Zeroes the values either way.
    virtual bool update_pattern(
        size_t num_stencils,
        const StencilGetter& stencil,
        size_t stale_block_tolerance,
        bool assume_unchanged) = 0;

    virtual void
    add(const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) = 0;

    virtual const Eigen::SparseMatrix<double>& to_eigen() const = 0;
    virtual Eigen::SparseMatrix<double> take_eigen() = 0;
    virtual const MeshFEM::BlockCSCHessianBase& block_matrix() const = 0;
};

/// Dimension-specific implementation (dim ∈ {2, 3}), mirroring MeshFEM's own
/// IPC integration (IPCWrapper.cc): one dim-sized block variable per vertex.
template <int dim>
struct MeshFEMHessianAssembler::Impl final
    : public MeshFEMHessianAssembler::ImplBase {
    using VarStructure = MeshFEM::OptimizationVarStructure<dim>;
    /// Local (block) variables of a collision stencil: 1–4 vertices.
    using Stencil = MeshFEM::ElementBlockVarsWithSizeRange<1, 4>;

    explicit Impl(const size_t num_block_vars)
        : m_assembler(num_block_vars)
        , m_vars(num_block_vars)
    {
        m_locks.init(num_block_vars);
    }

    int dimension() const override { return dim; }
    size_t num_block_vars() const override { return m_vars.numBlocks(); }

    bool update_pattern(
        const size_t num_stencils,
        const StencilGetter& stencil,
        const size_t stale_block_tolerance,
        const bool assume_unchanged) override
    {
        const auto block_vars_of_stencil = [&stencil](const size_t i) {
            return to_stencil(stencil(i));
        };

        if (m_H != nullptr) {
            // Caller-asserted fast path: skip change detection, which costs
            // as much as a pattern rebuild on large scenes. A differing
            // stencil count disproves the assumption, so fall through to
            // detection in that case.
            if (assume_unchanged && num_stencils == m_num_stencils) {
                // Verify the caller's claim where it is cheap to do so.
                assert(
                    m_assembler.detectChangedEntries(
                        *m_H, num_stencils, block_vars_of_stencil)
                    != MeshFEM::SystemAssembler<dim>::NEW_ENTRIES);
                m_H->setZero();
                return true;
            }

            IPC_TOOLKIT_PROFILE_BLOCK("MeshFEM detect pattern change");
            const size_t changed = m_assembler.detectChangedEntries(
                *m_H, num_stencils, block_vars_of_stencil);
            // NEW_ENTRIES (size_t max) must always trigger a rebuild:
            // scattering into a pattern that is missing an entry is undefined.
            if (changed != MeshFEM::SystemAssembler<dim>::NEW_ENTRIES
                && changed <= stale_block_tolerance) {
                m_num_stencils = num_stencils;
                m_H->setZero();
                return true;
            }
        }

        {
            IPC_TOOLKIT_PROFILE_BLOCK("MeshFEM block sparsity pattern");
            m_H = m_assembler.blockSparsityPattern(
                num_stencils, block_vars_of_stencil);
        }
        m_num_stencils = num_stencils;
        m_H->setZero(); // Allocate (and zero) the value array.
        m_eigen_structure_valid = false;
        return false;
    }

    void
    add(const MatrixMax12d& local_hess,
        const std::array<index_t, 4>& vertex_ids) override
    {
        assert(m_H != nullptr);
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
    // matrix. The structure (symmetrized pattern + index arrays) is cached
    // and reused while the block pattern is unchanged; only the values are
    // recomputed per call.
    //
    // NOTE: We deliberately do not use BlockCSCHessian::toEigen here: it
    // cannot reuse a cached structure across assemblies, and it goes through
    // an intermediate upper-triangle scalar matrix (serially) where this
    // implementation symmetrizes directly (in parallel) — roughly 2× faster
    // even cold.
    //
    // Value layout (ContiguousBlocks=true, StoreFullDiagonalBlocks=true, the
    // library default): block entry ii occupies Ax[N²·ii, N²·(ii+1)),
    // column-major within the block; diagonal blocks are full N×N with a
    // zeroed strict lower triangle.
    const Eigen::SparseMatrix<double>& to_eigen() const override
    {
        IPC_TOOLKIT_PROFILE_BLOCK("MeshFEM block CSC to Eigen");
        assert(m_H != nullptr);

        if (!m_eigen_structure_valid) {
            build_eigen_structure();
            m_eigen_structure_valid = true;
        }
        fill_eigen_values();
        return m_M;
    }

    Eigen::SparseMatrix<double> take_eigen() override
    {
        to_eigen();
        m_eigen_structure_valid = false; // m_M is about to be gutted
        return std::move(m_M);
    }

    const MeshFEM::BlockCSCHessianBase& block_matrix() const override
    {
        assert(m_H != nullptr);
        return *m_H;
    }

private:
    static constexpr int N = dim;

    /// A block of the full (symmetrized) matrix and the stored block backing
    /// it.
    struct BlockEntry {
        std::ptrdiff_t row; ///< Block row in the full (symmetrized) matrix
        std::ptrdiff_t src; ///< Index of the stored block (into Ai/Ax)
        bool transposed;    ///< Whether the stored block is mirrored
    };

    /// Whether the valid (non-negative) entries of `vertex_ids` are pairwise
    /// distinct. Only used to check to_stencil()'s precondition.
    static bool has_distinct_ids(const std::array<index_t, 4>& vertex_ids)
    {
        for (size_t i = 0; i < vertex_ids.size(); i++) {
            for (size_t j = i + 1; j < vertex_ids.size(); j++) {
                if (vertex_ids[i] >= 0 && vertex_ids[i] == vertex_ids[j]) {
                    return false;
                }
            }
        }
        return true;
    }

    /// Compact the (possibly -1-padded) vertex ID array into a MeshFEM
    /// stencil of block variables.
    ///
    /// @note The IDs must be pairwise distinct: we scatter with
    /// ElementHessianContribAssembler<true>, whose sorted column-merge assumes
    /// each block variable occurs once (MeshFEMSparse ships a separate
    /// ...SupportingStencilDuplicates for the duplicate case). A repeated ID
    /// would silently double-count the diagonal block and pollute its strict
    /// lower triangle. The broad phase rejects candidates that share a vertex,
    /// so this only guards hand-built collision sets.
    static Stencil to_stencil(const std::array<index_t, 4>& vertex_ids)
    {
        assert(has_distinct_ids(vertex_ids));
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

    void build_eigen_structure() const
    {
        const auto& Ap = m_H->Ap; // block column pointers (size n_blocks + 1)
        const auto& Ai = m_H->Ai; // block row indices

        const std::ptrdiff_t n_blocks = m_H->n;
        const std::ptrdiff_t n = N * n_blocks; // scalar dimension

        // --- Symmetrized *block* pattern --------------------------------
        // Every stored block (bi, bj) (bi ≤ bj) appears at (bi, bj) and, if
        // off-diagonal, mirrored at (bj, bi) in the full matrix.
        m_sym_col_start.assign(n_blocks + 1, 0);
        for (std::ptrdiff_t bj = 0; bj < n_blocks; bj++) {
            for (std::ptrdiff_t ii = Ap[bj]; ii < Ap[bj + 1]; ii++) {
                m_sym_col_start[bj + 1]++;
                if (Ai[ii] != bj) {
                    m_sym_col_start[Ai[ii] + 1]++;
                }
            }
        }
        std::partial_sum(
            m_sym_col_start.begin(), m_sym_col_start.end(),
            m_sym_col_start.begin());
        const std::ptrdiff_t num_sym_blocks = m_sym_col_start[n_blocks];

        // Sweeping bj in ascending order appends rows to every column in
        // ascending order: column c receives its upper entries (rows ≤ c) at
        // step bj = c and its mirrored entries (rows bj > c) at later steps.
        m_sym_entries.resize(num_sym_blocks);
        {
            std::vector<std::ptrdiff_t> cursor(
                m_sym_col_start.begin(), m_sym_col_start.end() - 1);
            for (std::ptrdiff_t bj = 0; bj < n_blocks; bj++) {
                for (std::ptrdiff_t ii = Ap[bj]; ii < Ap[bj + 1]; ii++) {
                    const std::ptrdiff_t bi = Ai[ii];
                    m_sym_entries[cursor[bj]++] = { bi, ii, false };
                    if (bi != bj) {
                        m_sym_entries[cursor[bi]++] = { bj, ii, true };
                    }
                }
            }
        }

        // --- Scalar CSC index arrays ------------------------------------
        m_M.resize(n, n);
        m_M.makeCompressed();
        m_M.resizeNonZeros(N * N * num_sym_blocks);

        int* outer = m_M.outerIndexPtr();
        int* inner = m_M.innerIndexPtr();

        for (std::ptrdiff_t bj = 0; bj <= n_blocks; bj++) {
            const std::ptrdiff_t col_nnz = bj < n_blocks
                ? N * (m_sym_col_start[bj + 1] - m_sym_col_start[bj])
                : 0;
            for (int c = 0; c < ((bj < n_blocks) ? N : 1); c++) {
                outer[N * bj + c] =
                    int(N * N * m_sym_col_start[bj] + c * col_nnz);
            }
        }

        tbb::parallel_for(
            std::ptrdiff_t(0), n_blocks, [&](const std::ptrdiff_t bj) {
                for (int c = 0; c < N; c++) {
                    std::ptrdiff_t out = outer[N * bj + c];
                    for (std::ptrdiff_t ei = m_sym_col_start[bj];
                         ei < m_sym_col_start[bj + 1]; ei++) {
                        for (int r = 0; r < N; r++) {
                            inner[out++] = int(N * m_sym_entries[ei].row + r);
                        }
                    }
                }
            });
    }

    void fill_eigen_values() const
    {
        const auto& Ax = m_H->Ax; // scalar values (N² per block)
        const std::ptrdiff_t n_blocks = m_H->n;

        const int* outer = m_M.outerIndexPtr();
        double* values = m_M.valuePtr();

        tbb::parallel_for(
            std::ptrdiff_t(0), n_blocks, [&](const std::ptrdiff_t bj) {
                for (int c = 0; c < N; c++) {
                    std::ptrdiff_t out = outer[N * bj + c];
                    for (std::ptrdiff_t ei = m_sym_col_start[bj];
                         ei < m_sym_col_start[bj + 1]; ei++) {
                        const BlockEntry& e = m_sym_entries[ei];
                        const double* block = Ax.data() + N * N * e.src;
                        const bool diagonal = e.row == bj;
                        for (int r = 0; r < N; r++) {
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
    }

    MeshFEM::SystemAssembler<dim> m_assembler;
    VarStructure m_vars;
    std::unique_ptr<MeshFEM::BlockCSCHessian<VarStructure>>
        m_H; // NOLINT(readability-identifier-naming): H = Hessian matrix
    MeshFEM::VarLocks m_locks;
    /// Stencil count of the cached pattern (assume-unchanged sanity check).
    size_t m_num_stencils = 0;

    // Cached Eigen conversion structure; valid while the pattern is reused.
    mutable bool m_eigen_structure_valid = false;
    mutable std::vector<std::ptrdiff_t> m_sym_col_start;
    mutable std::vector<BlockEntry> m_sym_entries;
    mutable Eigen::SparseMatrix<double>
        m_M; // NOLINT(readability-identifier-naming): M = matrix
};

MeshFEMHessianAssembler::MeshFEMHessianAssembler() = default;
MeshFEMHessianAssembler::~MeshFEMHessianAssembler() = default;

void MeshFEMHessianAssembler::begin(
    const int ndof,
    const int dim,
    const size_t num_stencils,
    const StencilGetter& stencil)
{
    // Validate before dividing by dim: dim == 0 would trap.
    if (dim != 2 && dim != 3) {
        log_and_throw_error(
            "MeshFEMHessianAssembler: unsupported dimension {}!", dim);
    }
    assert(ndof % dim == 0);
    const size_t num_block_vars = size_t(ndof) / dim;

    if (m_impl == nullptr || m_impl->dimension() != dim
        || m_impl->num_block_vars() != num_block_vars) {
        if (dim == 2) {
            m_impl = std::make_unique<Impl<2>>(num_block_vars);
        } else {
            m_impl = std::make_unique<Impl<3>>(num_block_vars);
        }
    }

    m_reused_pattern = m_impl->update_pattern(
        num_stencils, stencil, m_stale_block_tolerance,
        m_assume_unchanged_stencils);
}

void MeshFEMHessianAssembler::add_local_hessian(
    const MatrixMax12d& local_hess, const std::array<index_t, 4>& vertex_ids)
{
    assert(m_impl != nullptr);
    m_impl->add(local_hess, vertex_ids);
}

const Eigen::SparseMatrix<double>& MeshFEMHessianAssembler::get_matrix() const
{
    if (m_impl == nullptr) {
        log_and_throw_error(
            "MeshFEMHessianAssembler::get_matrix() called before the first "
            "assembly!");
    }
    return m_impl->to_eigen();
}

Eigen::SparseMatrix<double> MeshFEMHessianAssembler::take_matrix()
{
    if (m_impl == nullptr) {
        log_and_throw_error(
            "MeshFEMHessianAssembler::take_matrix() called before the first "
            "assembly!");
    }
    return m_impl->take_eigen();
}

const MeshFEM::BlockCSCHessianBase&
MeshFEMHessianAssembler::block_matrix() const
{
    if (m_impl == nullptr) {
        log_and_throw_error(
            "MeshFEMHessianAssembler::block_matrix() called before the first "
            "assembly!");
    }
    return m_impl->block_matrix();
}

} // namespace ipc

#endif // IPC_TOOLKIT_WITH_MESHFEM_SPARSE
