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
template <typename T> using Vector1 = Vector<T, 1>;
template <typename T> using Vector2 = Vector<T, 2>;
template <typename T> using Vector3 = Vector<T, 3>;
template <typename T> using VectorX = Vector<T, Dynamic>;
template <typename T> using Matrix2 = Matrix<T, 2, 2>;
template <typename T> using Matrix3 = Matrix<T, 3, 3>;
template <typename T> using MatrixX = Matrix<T, Dynamic, Dynamic>;

using Vector1d = Vector1<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
template <typename T> using VectorX3 = Vector<T, Dynamic, 3>;
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
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
template <typename T>
using MatrixXX3 = Matrix<T, Dynamic, Dynamic, ColMajor, 3, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
template <typename T>
using MatrixXX6 = Matrix<T, Dynamic, Dynamic, ColMajor, 6, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
using MatrixXX3d = MatrixXX3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 6
using MatrixXX6d = MatrixXX6<double>;
/// @brief A dynamic size diagonal matrix
using DiagonalMatrixXd = DiagonalMatrix<double, Dynamic>;
/// @brief A dynamic size diagonal matrix with a fixed maximum size of 6 x 6
using DiagonalMatrixX6d = DiagonalMatrix<double, Dynamic, 6>;

/// Matrix Projection onto Positive Definite Cone
MatrixXd project_to_pd(const MatrixXd& A, double eps = 1e-8);
/// Matrix Projection onto Positive Semi-Definite Cone
MatrixXd project_to_psd(const MatrixXd& A);

/// Cross product for dynamically sized matrices.
template <typename DerivedA, typename DerivedB>
auto cross(const MatrixBase<DerivedA>& a, const MatrixBase<DerivedB>& b)
{
    assert(a.size() == 3 && b.size() == 3);
    Eigen::Matrix<
        typename DerivedA::Scalar, DerivedA::RowsAtCompileTime,
        DerivedA::ColsAtCompileTime>
        c(a.rows(), a.cols());
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

} // namespace Eigen
