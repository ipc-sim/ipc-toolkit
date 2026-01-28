#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc::affine {

struct AffineBody {
    MatrixMax3d A; // Affine matrix
    VectorMax3d p; // Translation vector
    double volume; // Volume of the body

    AffineBody() = default;

    AffineBody(
        Eigen::ConstRef<MatrixMax3d> A,
        Eigen::ConstRef<VectorMax3d> p,
        const double volume)
        : A(A)
        , p(p)
        , volume(volume)
    {
    }

    // NOLINTNEXTLINE(readability-identifier-naming)
    Eigen::SparseMatrix<double> J(const Eigen::MatrixXd& rest_positions) const
    {
        std::vector<Eigen::Triplet<double>> triplets;

        for (int i = 0; i < rest_positions.rows(); i++) {
            for (int j = 0; j < rest_positions.cols(); j++) {
                triplets.emplace_back(rest_positions.rows() * j + i, j, 1);
            }
        }

        // I ⊗ x̄
        for (int i = 0; i < rest_positions.rows(); i++) {
            for (int j = 0; j < rest_positions.cols(); j++) {
                for (int k = 0; k < A.cols(); k++) {
                    triplets.emplace_back(
                        i + k * rest_positions.rows(),
                        j + k * rest_positions.cols() + p.size(),
                        rest_positions(i, j));
                }
            }
        }

        Eigen::SparseMatrix<double> J(
            rest_positions.size(), A.size() + p.size());
        J.setFromTriplets(triplets.begin(), triplets.end());
        return J;
    }

    Eigen::MatrixXd
    transform_vertices(const Eigen::MatrixXd& rest_positions) const
    {
        // Compute: A x̄ + p
        // transpose because x is row-ordered
        return (rest_positions * A.transpose()).rowwise() + p.transpose();
    }

    Eigen::SparseMatrix<double>
    transform_vertices_gradient(const Eigen::MatrixXd& rest_positions) const
    {
        return J(rest_positions);
    }
};

} // namespace ipc::affine