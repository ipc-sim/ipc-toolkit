#include <ipc/utils/eigen_ext.hpp>

#include <iostream>

#include <Eigen/Eigenvalues>

#ifdef IPC_TOOLKIT_WITH_LOGGER
#include <ipc/utils/logger.hpp>
#endif

namespace Eigen {

// Matrix Projection onto Positive Definite Cone
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
    double eps)
{
    assert(eps > 0);

    // https://math.stackexchange.com/q/2776803
    SelfAdjointEigenSolver<
        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>
        eigensolver(A);
    if (eigensolver.info() != Success) {
#ifdef IPC_TOOLKIT_WITH_LOGGER
        ipc::logger().error(
            "unable to project matrix onto positive definite cone");
#else
        throw "unable to project matrix onto positive definite cone";
#endif
        return A;
    }
    // Check if all eigen values are positive.
    // The eigenvalues are sorted in increasing order.
    if (eigensolver.eigenvalues()[0] > 0.0) {
        return A;
    }
    DiagonalMatrix<double, Dynamic> D(eigensolver.eigenvalues());
    // Save a little time and only project the negative or zero values
    for (int i = 0; i < A.rows(); i++) {
        if (D.diagonal()[i] <= 0.0) {
            D.diagonal()[i] = eps;
        } else {
            break;
        }
    }
    return eigensolver.eigenvectors() * D
        * eigensolver.eigenvectors().transpose();
}

// Matrix Projection onto Positive Semi-Definite Cone
template <
    typename _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols>
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
project_to_psd(
    const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& A)
{
    // https://math.stackexchange.com/q/2776803
    SelfAdjointEigenSolver<
        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>
        eigensolver(A);
    if (eigensolver.info() != Success) {
#ifdef IPC_TOOLKIT_WITH_LOGGER
        ipc::logger().error(
            "unable to project matrix onto positive semi-definite cone");
#else
        throw "unable to project matrix onto positive semi-definite cone";
#endif
        return A;
    }
    // Check if all eigen values are zero or positive.
    // The eigenvalues are sorted in increasing order.
    if (eigensolver.eigenvalues()[0] >= 0.0) {
        return A;
    }
    DiagonalMatrix<double, Dynamic> D(eigensolver.eigenvalues());
    // Save a little time and only project the negative values
    for (int i = 0; i < A.rows(); i++) {
        if (D.diagonal()[i] < 0.0) {
            D.diagonal()[i] = 0.0;
        } else {
            break;
        }
    }
    return eigensolver.eigenvectors() * D
        * eigensolver.eigenvectors().transpose();
}

template <typename DerivedA, typename DerivedB, typename Result>
Result cross(const MatrixBase<DerivedA>& a, const MatrixBase<DerivedB>& b)
{
    assert(a.size() == 3 && b.size() == 3);
    Result c(a.rows(), a.cols());
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

} // namespace Eigen
