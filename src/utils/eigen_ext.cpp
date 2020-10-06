#include <ipc/utils/eigen_ext.hpp>

#include <iostream>

#include <Eigen/Eigenvalues>

#ifdef IPC_TOOLKIT_WITH_LOGGER
#include <ipc/utils/logger.hpp>
#endif

namespace Eigen {

// Matrix Projection onto Positive Definite Cone
MatrixXd project_to_pd(const MatrixXd& A, double eps)
{
    assert(eps > 0);

    // https://math.stackexchange.com/q/2776803
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
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
MatrixXd project_to_psd(const MatrixXd& A)
{
    // https://math.stackexchange.com/q/2776803
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
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

} // namespace Eigen
