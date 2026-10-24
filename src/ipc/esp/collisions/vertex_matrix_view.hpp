#pragma once

#include <ipc/math/math.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <cassert>
#include <type_traits>

namespace ipc {

/// @brief A non-owning view that presents one or two matrices (with the same
/// number of columns) as a single vertically concatenated matrix.
///
/// The view reads each matrix through its own row and column strides, so any
/// directly accessible double expression works without a copy: plain
/// matrices of either storage order, blocks such as `V.topRows(n)`, strided
/// `Eigen::Map`s, and `Eigen::Ref`s bound to any of these.
///
/// We take each argument as a template instead of an `Eigen::Ref<const
/// MatrixXd>` because a const Ref silently copies arguments it cannot map
/// (row-major matrices, arithmetic expressions) into a temporary that dies at
/// the end of the constructor call, which would leave the view dangling.
/// Instead, expressions without direct access fail to compile, and so do
/// temporary plain matrices, since a view of them would dangle the same way.
template <int ncols = 3> class VertexMatrixView {
public:
    /// @brief Construct a view concatenating two matrices vertically.
    /// @param A The top matrix.
    /// @param B The bottom matrix.
    template <typename DerivedA, typename DerivedB>
    VertexMatrixView(
        const Eigen::DenseBase<DerivedA>& A,
        const Eigen::DenseBase<DerivedB>& B)
        : m_n_a_rows(A.rows())
        , m_n_b_rows(B.rows())
        , m_a(A.derived().data())
        , m_b(B.derived().data())
        , m_a_row_stride(A.derived().rowStride())
        , m_a_col_stride(A.derived().colStride())
        , m_b_row_stride(B.derived().rowStride())
        , m_b_col_stride(B.derived().colStride())
    {
        check_viewable<DerivedA>();
        check_viewable<DerivedB>();
        if (A.cols() != ncols || B.cols() != ncols) {
            log_and_throw_error("Incompatible matrix columns!");
        }
    }

    /// @brief Construct a view wrapping a single matrix (no concatenation).
    template <typename DerivedA>
    explicit VertexMatrixView(const Eigen::DenseBase<DerivedA>& A)
        : m_n_a_rows(A.rows())
        , m_n_b_rows(0)
        , m_a(A.derived().data())
        , m_b(nullptr)
        , m_a_row_stride(A.derived().rowStride())
        , m_a_col_stride(A.derived().colStride())
        , m_b_row_stride(0)
        , m_b_col_stride(0)
    {
        check_viewable<DerivedA>();
        if (A.cols() != ncols) {
            log_and_throw_error("Incompatible matrix columns!");
        }
    }

    // A temporary matrix is destroyed at the end of the full-expression that
    // creates it, so a view of one would dangle as soon as it is constructed.
    template <typename DerivedA>
    explicit VertexMatrixView(const Eigen::PlainObjectBase<DerivedA>&&) =
        delete;
    template <typename DerivedA, typename DerivedB>
    VertexMatrixView(
        const Eigen::PlainObjectBase<DerivedA>&&,
        const Eigen::DenseBase<DerivedB>&) = delete;
    template <typename DerivedA, typename DerivedB>
    VertexMatrixView(
        const Eigen::DenseBase<DerivedA>&,
        const Eigen::PlainObjectBase<DerivedB>&&) = delete;

    /// @brief Access row i of the concatenated matrix.
    Eigen::RowVector<double, ncols> operator()(index_t i) const
    {
        assert(i < rows());
        const bool in_a = i < m_n_a_rows;
        const double* src = in_a ? m_a : m_b;
        const index_t li = in_a ? i : (i - m_n_a_rows);
        const Eigen::Index row_stride = in_a ? m_a_row_stride : m_b_row_stride;
        const Eigen::Index col_stride = in_a ? m_a_col_stride : m_b_col_stride;
        Eigen::RowVector<double, ncols> row;
        for (int d = 0; d < ncols; ++d) {
            row[d] = src[li * row_stride + d * col_stride];
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

private:
    template <typename Derived> static constexpr void check_viewable()
    {
        static_assert(
            std::is_same_v<typename Derived::Scalar, double>,
            "VertexMatrixView requires a matrix of doubles");
        static_assert(
            (int(Derived::Flags) & Eigen::DirectAccessBit) != 0,
            "VertexMatrixView requires direct access to the matrix "
            "coefficients; evaluate the expression into a matrix first");
    }

    // Element (i, j) of a matrix lives at data[i * row_stride + j * col_stride]
    // which covers both storage orders and any outer stride.
    const Eigen::Index m_a_row_stride;
    const Eigen::Index m_a_col_stride;
    const Eigen::Index m_b_row_stride;
    const Eigen::Index m_b_col_stride;
};

} // namespace ipc
