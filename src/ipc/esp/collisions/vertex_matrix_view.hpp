#pragma once

#include <ipc/math/math.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <cassert>

namespace ipc {

/// @brief A non-owning view that presents one or two column-major matrices
/// (with the same number of columns) as a single vertically concatenated
/// matrix.
template <int ncols = 3> class VertexMatrixView {
public:
    /// @brief Construct a view concatenating two matrices vertically.
    /// @param A The top matrix.
    /// @param B The bottom matrix.
    VertexMatrixView(
        Eigen::ConstRef<Eigen::MatrixXd> A, Eigen::ConstRef<Eigen::MatrixXd> B)
        : m_n_a_rows(A.rows())
        , m_n_b_rows(B.rows())
        , m_a(A.data())
        , m_b(B.data())
    {
        if (A.cols() != ncols || B.cols() != ncols) {
            log_and_throw_error("Incompatible matrix columns!");
        }
    }

    /// @brief Construct a view wrapping a single matrix (no concatenation).
    explicit VertexMatrixView(Eigen::ConstRef<Eigen::MatrixXd> A)
        : m_n_a_rows(A.rows())
        , m_n_b_rows(0)
        , m_a(A.data())
        , m_b(nullptr)
    {
        if (A.cols() != ncols) {
            log_and_throw_error("Incompatible matrix columns!");
        }
    }

    /// @brief Access row i of the concatenated matrix.
    Eigen::RowVector<double, ncols> operator()(index_t i) const
    {
        assert(i < rows());
        Eigen::RowVector<double, ncols> row;
        const double* src = (i < m_n_a_rows) ? m_a : m_b;
        const index_t nrows = (i < m_n_a_rows) ? m_n_a_rows : m_n_b_rows;
        const index_t li = (i < m_n_a_rows) ? i : (i - m_n_a_rows);
        for (int d = 0; d < ncols; ++d) {
            row[d] = src[li + d * nrows];
        }
        return row;
    }

    /// @brief Total number of rows (A rows + B rows).
    index_t rows() const { return m_n_a_rows + m_n_b_rows; }

    /// @brief Number of columns (compile-time constant).
    index_t cols() const { return ncols; }

    const index_t m_n_a_rows;
    const index_t m_n_b_rows;
    const double* const m_a;
    const double* const m_b;
};

} // namespace ipc
