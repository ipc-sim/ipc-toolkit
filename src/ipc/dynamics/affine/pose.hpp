#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/SparseCore>

namespace ipc::affine {

/// @brief An affine body pose: world = A x̄ + p.
///
/// @p rotation is a dim×dim matrix A — a rotation matrix for rigid bodies (RBD)
/// or a general affine map for affine bodies (ABD). This is the shared
/// "affine-coordinate" representation of a body's configuration.
struct Pose {
    /// @brief Position (translation) of the body.
    VectorMax3d position;
    /// @brief The dim×dim affine/rotation matrix A.
    MatrixMax3d rotation;

    /// @brief Recover the rotation vector (3D) or angle (2D) from A.
    /// @note Assumes A ∈ SO(dim). This holds for rigid bodies but NOT for
    ///       affine bodies, whose A is only softly orthogonal — do not call
    ///       this on a raw affine A. Use the affine pose directly (e.g.
    ///       transform_vertices) instead of projecting to a rigid rotation.
    /// @return The rotation vector (3D) or 1-vector angle (2D).
    VectorMax3d rotation_vector() const;

    /// @brief Set A from a rotation vector (3D) or angle (2D).
    /// @param theta The rotation vector (3D) or angle (2D).
    void set_rotation_vector(Eigen::ConstRef<VectorMax3d> theta);

    /// @brief Transform vertices from local to world coordinates (A x̄ + p).
    /// @param V The vertices to transform (each row is a vertex).
    /// @return The transformed vertices.
    Eigen::MatrixXd transform_vertices(Eigen::ConstRef<Eigen::MatrixXd> V) const
    {
        // Compute: A x̄ + p
        // transpose because V is row-ordered
        return (V * rotation.transpose()).rowwise() + position.transpose();
    }

    /// @brief Compute the Jacobian of the transformed vertices w.r.t. the DOFs.
    /// @param V The vertices to transform (each row is a vertex).
    /// @return The Jacobian matrix of size (num_vertices * dim) x ndof.
    Eigen::SparseMatrix<double>
    transform_vertices_jacobian(Eigen::ConstRef<Eigen::MatrixXd> V) const
    {
        return J(V);
    }

    /// @brief Compute the Jacobian of the transformed vertices w.r.t. the DOFs.
    ///
    /// DOF layout: q = [p (dim); vec(A) column-major (dim²)], i.e.,
    /// q[dim + k + dim·j] = A(k, j) — matching the integrator's state layout
    /// and OrthogonalityPotential. Rows are vertex-major (dim·i + k), matching
    /// the collision-mesh gradient layout.
    ///
    /// @param rest_positions The rest vertices (each row is a vertex).
    /// @return The Jacobian matrix of size (num_vertices * dim) x (ndof).
    static Eigen::SparseMatrix<double>
    J(Eigen::ConstRef<Eigen::MatrixXd> rest_positions); // NOLINT
};

} // namespace ipc::affine
