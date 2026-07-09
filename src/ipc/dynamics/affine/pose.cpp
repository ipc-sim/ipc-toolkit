#include "pose.hpp"

#include <ipc/dynamics/rigid/pose.hpp>

namespace ipc::affine {

VectorMax3d Pose::rotation_vector() const
{
    assert(rotation.rows() == rotation.cols());
    assert(rotation.rows() == 2 || rotation.rows() == 3);
    // A should be (approximately) a rotation. For rigid bodies (RBD) A = exp(θ)
    // is orthogonal to machine precision, but for affine bodies (ABD) A is only
    // *softly* held in SO(dim) by the orthogonality potential, so it drifts
    // from unitary by a small amount (~1e-5 in the near-rigid regime these
    // bodies operate in). Use a loose tolerance that still catches a genuine
    // non-rotation (real affine deformation drifts by orders of magnitude more).
    assert(rotation.isUnitary(1e-3));
    if (rotation.rows() == 2) {
        // For 2D, return the angle
        VectorMax3d angle(1);
        // rotation(1, 0) = sin(θ), rotation(0, 0) = cos(θ)
        // Thus, θ = atan2(sin(θ), cos(θ))
        angle(0) = std::atan2(rotation(1, 0), rotation(0, 0));
        return angle;
    } else {
        // For 3D, return the rotation vector
        return rigid::rotation_matrix_to_vector(rotation);
    }
}

void Pose::set_rotation_vector(Eigen::ConstRef<VectorMax3d> theta)
{
    assert(theta.size() == 1 || theta.size() == 3);
    if (theta.size() == 1) {
        // For 2D, set the rotation matrix directly
        rotation.resize(2, 2);
        rotation << std::cos(theta(0)), -std::sin(theta(0)), //
            std::sin(theta(0)), std::cos(theta(0));
    } else {
        // For 3D, convert the rotation vector to a rotation matrix
        rotation = rigid::rotation_vector_to_matrix(theta);
    }
}

Eigen::SparseMatrix<double>
Pose::J(Eigen::ConstRef<Eigen::MatrixXd> rest_positions)
{
    const int dim = rest_positions.cols();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(rest_positions.rows() * dim * (dim + 1));

    for (int i = 0; i < rest_positions.rows(); i++) {
        for (int k = 0; k < dim; k++) {
            // ∂(world_{i,k}) / ∂p_k = 1
            triplets.emplace_back(dim * i + k, k, 1);
            // world_{i,k} = p_k + Σ_j A(k, j) x̄(i, j)
            for (int j = 0; j < dim; j++) {
                triplets.emplace_back(
                    dim * i + k, dim + k + dim * j, rest_positions(i, j));
            }
        }
    }

    Eigen::SparseMatrix<double> J(rest_positions.size(), dim + dim * dim);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

} // namespace ipc::affine
