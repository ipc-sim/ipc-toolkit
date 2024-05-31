#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc::affine {

struct AffineBody {
    MatrixMax3d A; // Affine matrix
    VectorMax3d p; // Translation vector
    double volume; // Volume of the body

    AffineBody() = default;

    AffineBody(const MatrixMax3d& A, const VectorMax3d& p, const double volume)
        : A(A)
        , p(p)
        , volume(volume)
    {
    }

    Eigen::MatrixXd
    transform_vertices(const Eigen::MatrixXd& rest_positions) const
    {
        // Compute: A x̄ + p
        // transpose because x is row-ordered
        return (rest_positions * A.transpose()).rowwise() + p.transpose();
    }
};

} // namespace ipc::affine