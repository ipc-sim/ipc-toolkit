#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace ipc {

// Boolean scalar
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;
typedef Eigen::Matrix<bool, 3, 1> Vector3b;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

template <typename T, int dim, int max_dim = dim>
using Vector = Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1>;
template <typename T, int dim, int max_dim = dim>
using RowVector = Eigen::Matrix<T, 1, dim, Eigen::RowMajor, 1, max_dim>;
template <typename T> using Vector1 = Vector<T, 1>;
template <typename T> using Vector2 = Vector<T, 2>;
template <typename T> using Vector3 = Vector<T, 3>;
template <typename T> using VectorX = Vector<T, Eigen::Dynamic>;
template <typename T> using RowVector2 = RowVector<T, 2>;
template <typename T> using RowVector3 = RowVector<T, 3>;
template <typename T> using RowVectorX = RowVector<T, Eigen::Dynamic>;
template <typename T> using Matrix2 = Eigen::Matrix<T, 2, 2>;
template <typename T> using Matrix3 = Eigen::Matrix<T, 3, 3>;
template <typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

using Vector1d = Vector1<double>;
using Vector6d = Vector<double, 6>;
using Vector9d = Vector<double, 9>;
using Vector12d = Vector<double, 12>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Matrix9d = Eigen::Matrix<double, 9, 9>;
using Matrix12d = Eigen::Matrix<double, 12, 12>;

/// @brief A dynamic size matrix with a fixed maximum size of 3×1
template <typename T> using VectorMax2 = Vector<T, Eigen::Dynamic, 2>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×1
template <typename T> using VectorMax3 = Vector<T, Eigen::Dynamic, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 2×1
using VectorMax2d = VectorMax2<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×1
using VectorMax3d = VectorMax3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×1
using VectorMax3i = VectorMax3<int>;
/// @brief A dynamic size matrix with a fixed maximum size of 6×1
template <typename T> using VectorMax6 = Vector<T, Eigen::Dynamic, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 6×1
using VectorMax6d = VectorMax6<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6×1
using VectorMax6b = VectorMax6<bool>;
/// @brief A dynamic size matrix with a fixed maximum size of 9×1
template <typename T> using VectorMax9 = Vector<T, Eigen::Dynamic, 9>;
/// @brief A dynamic size matrix with a fixed maximum size of 9×1
using VectorMax9d = VectorMax9<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 12×1
template <typename T> using VectorMax12 = Vector<T, Eigen::Dynamic, 12>;
/// @brief A dynamic size matrix with a fixed maximum size of 12×1
using VectorMax12d = VectorMax12<double>;

/// @brief A dynamic size matrix with a fixed maximum size of 1×2
template <typename T> using RowVectorMax2 = RowVector<T, Eigen::Dynamic, 2>;
/// @brief A dynamic size matrix with a fixed maximum size of 1×3
template <typename T> using RowVectorMax3 = RowVector<T, Eigen::Dynamic, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 1×2
using RowVectorMax2d = RowVectorMax2<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 1×3
using RowVectorMax3d = RowVectorMax3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6×1
using RowVectorMax6d = RowVector<double, Eigen::Dynamic, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 9×1
using RowVectorMax9d = RowVector<double, Eigen::Dynamic, 9>;
/// @brief A dynamic size matrix with a fixed maximum size of 12×1
using RowVectorMax12d = RowVector<double, Eigen::Dynamic, 12>;

template <typename T, int max_rows, int max_cols>
using MatrixMax = Eigen::Matrix<
    T,
    Eigen::Dynamic,
    Eigen::Dynamic,
    Eigen::ColMajor,
    max_rows,
    max_cols>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×3
template <typename T> using MatrixMax2 = MatrixMax<T, 2, 2>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×3
template <typename T> using MatrixMax3 = MatrixMax<T, 3, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 6×6
template <typename T> using MatrixMax6 = MatrixMax<T, 6, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 9×9
template <typename T> using MatrixMax9 = MatrixMax<T, 9, 9>;
/// @brief A dynamic size matrix with a fixed maximum size of 12×12
template <typename T> using MatrixMax12 = MatrixMax<T, 12, 12>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×3
using MatrixMax2d = MatrixMax2<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3×3
using MatrixMax3d = MatrixMax3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6×6
using MatrixMax6d = MatrixMax6<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 12×12
using MatrixMax9d = MatrixMax9<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 12×12
using MatrixMax12d = MatrixMax12<double>;
/// @brief A dynamic size diagonal matrix
using DiagonalMatrixXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;
/// @brief A dynamic size diagonal matrix with a fixed maximum size of 6×6
using DiagonalMatrixMax6d = Eigen::DiagonalMatrix<double, Eigen::Dynamic, 6>;

/// @brief A dynamic size array with a fixed maximum size of 2×1
template <typename T>
using ArrayMax2 = Eigen::Array<T, Eigen::Dynamic, 1, Eigen::ColMajor, 2, 1>;
/// @brief A dynamic size array with a fixed maximum size of 2×1
template <typename T>
using ArrayMax3 = Eigen::Array<T, Eigen::Dynamic, 1, Eigen::ColMajor, 3, 1>;
/// @brief A dynamic size array with a fixed maximum size of 2×1
using ArrayMax2d = ArrayMax2<double>;
/// @brief A dynamic size array with a fixed maximum size of 3×1
using ArrayMax3d = ArrayMax3<double>;
/// @brief A dynamic size array with a fixed maximum size of 3×1
using ArrayMax3i = ArrayMax3<int>;

/// @brief Matrix projection onto positive definite cone
/// @param A Symmetric matrix to project
/// @return Projected matrix
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

/// @brief Matrix projection onto positive semi-definite cone
/// @param A Symmetric matrix to project
/// @return Projected matrix
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

inline Eigen::Vector3d to_3D(const VectorMax3d& v)
{
    assert(v.size() == 2 || v.size() == 3);
    return v.size() == 2 ? Eigen::Vector3d(v.x(), v.y(), 0) : v.head<3>();
}

/// Eigen IO Format to format vectors like vertex rows in an OBJ file.
static const Eigen::IOFormat OBJ_VERTEX_FORMAT = Eigen::IOFormat(
    Eigen::FullPrecision, Eigen::DontAlignCols, " ", "", "v ", "\n", "", "");

} // namespace ipc

#include "eigen_ext.tpp"
