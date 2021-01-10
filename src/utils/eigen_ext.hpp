#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Eigen {

typedef Array<bool, Dynamic, 1> ArrayXb;
typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<bool, 3, 1> Vector3b;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

template <typename T, int dim, int max_dim = dim>
using Vector = Matrix<T, dim, 1, ColMajor, max_dim, 1>;
template <typename T, int dim, int max_dim = dim>
using RowVector = Matrix<T, 1, dim, RowMajor, 1, max_dim>;
template <typename T> using Vector1 = Vector<T, 1>;
template <typename T> using Vector2 = Vector<T, 2>;
template <typename T> using Vector3 = Vector<T, 3>;
template <typename T> using RowVector2 = RowVector<T, 2>;
template <typename T> using RowVector3 = RowVector<T, 3>;
template <typename T> using VectorX = Vector<T, Dynamic>;
template <typename T> using Matrix2 = Matrix<T, 2, 2>;
template <typename T> using Matrix3 = Matrix<T, 3, 3>;
template <typename T> using MatrixX = Matrix<T, Dynamic, Dynamic>;

using Vector1d = Vector1<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
template <typename T> using VectorX2 = Vector<T, Dynamic, 2>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
template <typename T> using VectorX3 = Vector<T, Dynamic, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 2 x 1
using VectorX2d = VectorX2<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
using VectorX3d = VectorX3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
using VectorX3i = VectorX3<int>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 1
template <typename T> using VectorX6 = Vector<T, Dynamic, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 1
using VectorX6d = VectorX6<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 1
using VectorX6b = VectorX6<bool>;
/// @brief A dynamic size matrix with a fixed maximum size of 9 x 1
template <typename T> using VectorX9 = Vector<T, Dynamic, 9>;
/// @brief A dynamic size matrix with a fixed maximum size of 9 x 1
using VectorX9d = VectorX9<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 12 x 1
template <typename T> using VectorX12 = Vector<T, Dynamic, 12>;
/// @brief A dynamic size matrix with a fixed maximum size of 12 x 1
using VectorX12d = VectorX12<double>;
template <typename T, int max_rows, int max_cols>
using MatrixXX = Matrix<T, Dynamic, Dynamic, ColMajor, max_rows, max_cols>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
template <typename T> using MatrixXX3 = MatrixXX<T, 3, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 6
template <typename T> using MatrixXX6 = MatrixXX<T, 6, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 9 x 9
template <typename T> using MatrixXX9 = MatrixXX<T, 9, 9>;
/// @brief A dynamic size matrix with a fixed maximum size of 12 x 12
template <typename T> using MatrixXX12 = MatrixXX<T, 12, 12>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
using MatrixXX3d = MatrixXX3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 6
using MatrixXX6d = MatrixXX6<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 12 x 12
using MatrixXX9d = MatrixXX9<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 12 x 12
using MatrixXX12d = MatrixXX12<double>;
/// @brief A dynamic size diagonal matrix
using DiagonalMatrixXd = DiagonalMatrix<double, Dynamic>;
/// @brief A dynamic size diagonal matrix with a fixed maximum size of 6 x 6
using DiagonalMatrixX6d = DiagonalMatrix<double, Dynamic, 6>;

/// Matrix Projection onto Positive Definite Cone
template <
    typename _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols>
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
project_to_pd(
    const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& A,
    double eps = 1e-8);
/// Matrix Projection onto Positive Semi-Definite Cone
template <
    typename _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols>
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
project_to_psd(
    const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>&
        A);

/// Cross product for dynamically sized matrices.
template <typename DerivedA, typename DerivedB>
auto cross(const MatrixBase<DerivedA>& a, const MatrixBase<DerivedB>& b);

} // namespace Eigen

#include "eigen_ext.tpp"
