#pragma once

#include "eigen_ext.hpp"

#include <ipc/utils/logger.hpp>

#include <Eigen/Eigenvalues>

#include <stdexcept> // std::runtime_error

namespace ipc {

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
    assert(A.isApprox(A.transpose()) && "A must be symmetric");

    // https://math.stackexchange.com/q/2776803
    Eigen::SelfAdjointEigenSolver<
        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>
        eigensolver(A);
    if (eigensolver.info() != Eigen::Success) {
        logger().error("unable to project matrix onto positive definite cone");
        throw std::runtime_error(
            "unable to project matrix onto positive definite cone");
    }
    // Check if all eigen values are positive.
    // The eigenvalues are sorted in increasing order.
    if (eigensolver.eigenvalues()[0] > 0.0) {
        return A;
    }
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(eigensolver.eigenvalues());
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
    const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& A,
    const PSDProjectionMethod method)
{
    assert(A.isApprox(A.transpose()) && "A must be symmetric");

    if (method == PSDProjectionMethod::NONE)
        return A;

    // https://math.stackexchange.com/q/2776803
    Eigen::SelfAdjointEigenSolver<
        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>
        eigensolver(A);
    if (eigensolver.info() != Eigen::Success) {
        logger().error(
            "unable to project matrix onto positive semi-definite cone");
        throw std::runtime_error(
            "unable to project matrix onto positive definite cone");
    }
    // Check if all eigen values are zero or positive.
    // The eigenvalues are sorted in increasing order.
    if (eigensolver.eigenvalues()[0] >= 0.0) {
        return A;
    }
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(eigensolver.eigenvalues());
    // Save a little time and only project the negative values
    for (int i = 0; i < A.rows(); i++) {
        if (D.diagonal()[i] < 0.0) {
            switch (method) {
            case PSDProjectionMethod::CLAMP:
                D.diagonal()[i] = 0.0;
                break;
            case PSDProjectionMethod::ABS:
                D.diagonal()[i] = std::abs(D.diagonal()[i]);
                break;
            default:
                throw std::runtime_error("Invalid type of PSD projection!");
            }
        } else {
            break;
        }
    }

    return eigensolver.eigenvectors() * D
        * eigensolver.eigenvectors().transpose();
}

} // namespace ipc
