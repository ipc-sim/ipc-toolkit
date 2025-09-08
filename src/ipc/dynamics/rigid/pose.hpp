#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc::rigid {

/// @brief Convert from a 3D rotation vector to a rotation matrix.
/// @param theta The rotation vector
/// @return The rotation matrix corresponding to the rotation vector
Eigen::Matrix3d
rotation_vector_to_matrix(Eigen::ConstRef<Eigen::Vector3d> theta);

/// @brief Compute the Jacobian of the rotation matrix with respect to the rotation vector.
/// @param theta The rotation vector
/// @return The Jacobian matrix
Eigen::Matrix<double, 9, 3>
rotation_vector_to_matrix_jacobian(Eigen::ConstRef<Eigen::Vector3d> theta);

/// @brief Compute the Hessian of the rotation matrix with respect to the rotation vector.
///
/// This is a 9x9 matrix where each row corresponds to the second derivative of
/// each element of the rotation matrix with respect to each component of the
/// rotation vector.
///
/// @param theta The rotation vector
/// @return The Hessian matrix
Eigen::Matrix<double, 9, 9>
rotation_vector_to_matrix_hessian(Eigen::ConstRef<Eigen::Vector3d> theta);

/// @brief Convert from a 3D rotation matrix to a rotation vector.
/// @param R The rotation matrix
/// @return The rotation vector corresponding to the rotation matrix
Eigen::Vector3d rotation_matrix_to_vector(Eigen::ConstRef<Eigen::Matrix3d> R);

// ----------------------------------------------------------------------------

struct Pose {
    // Position of the rigid body
    VectorMax3d position;
    // Rotation of the rigid body (rotation vector for 3D, angle for 1D)
    VectorMax3d rotation;

    Pose() = default;

    Pose(
        Eigen::ConstRef<VectorMax3d> _position,
        Eigen::ConstRef<VectorMax3d> _rotation)
        : position(_position)
        , rotation(_rotation)
    {
    }

    static Pose Zero(const int dim)
    {
        return Pose(VectorMax3d::Zero(dim), VectorMax3d::Zero(dim));
    }

    MatrixMax3d rotation_matrix() const
    {
        return rotation_vector_to_matrix(rotation);
    }

    Eigen::Quaternion<double> quaternion() const
    {
        assert(rotation.size() == 3);
        double angle = rotation.norm();
        Eigen::Vector3d axis = rotation / angle;
        return Eigen::Quaternion<double>(Eigen::AngleAxis<double>(angle, axis));
    }

    Eigen::MatrixXd transform_vertices(Eigen::ConstRef<Eigen::MatrixXd> V) const
    {
        // Compute: R(θ) V + p
        // transpose because x is row-ordered
        return (V * rotation_matrix().transpose()).rowwise()
            + position.transpose();
    }

    friend Pose operator*(const Pose& a, const Pose& b)
    {
        Pose c;

        c.position = a.position + b.position;

        if (c.position.size() == 2) {
            c.rotation = a.rotation + b.rotation;
        } else {
            MatrixMax3d R = a.rotation_matrix() * b.rotation_matrix();
            c.rotation = rotation_matrix_to_vector(R);
        }

        return c;
    }
};

struct AffinePose {
    // Position of the rigid body
    VectorMax3d position;
    // Rotation of the rigid body (rotation vector for 3D, angle for 1D)
    MatrixMax3d rotation;

    /// @brief Construct a rotation vector from the rotation matrix.
    /// @return The rotation vector corresponding to the rotation matrix
    VectorMax3d rotation_vector() const
    {
        assert(rotation.rows() == rotation.cols());
        assert(rotation.rows() == 2 || rotation.rows() == 3);
        assert(rotation.isUnitary(1e-9)); // Ensure it's a rotation matrix
        if (rotation.rows() == 2) {
            // For 2D, return the angle
            VectorMax3d angle(1);
            // rotation(1, 0) = sin(θ), rotation(0, 0) = cos(θ)
            // Thus, θ = atan2(sin(θ), cos(θ))
            angle(0) = std::atan2(rotation(1, 0), rotation(0, 0));
            return angle;
        } else {
            // For 3D, return the rotation vector
            return rotation_matrix_to_vector(rotation);
        }
    }

    /// @brief Set the rotation matrix from a rotation vector or angle.
    /// @param theta The rotation vector (3D) or angle (2D).
    void set_rotation_vector(Eigen::ConstRef<VectorMax3d> theta)
    {
        assert(theta.size() == 1 || theta.size() == 3);
        if (theta.size() == 1) {
            // For 2D, set the rotation matrix directly
            rotation << std::cos(theta(0)), -std::sin(theta(0)),
                std::sin(theta(0)), std::cos(theta(0));
        } else {
            // For 3D, convert the rotation vector to a rotation matrix
            rotation = rotation_vector_to_matrix(theta);
        }
    }
};

} // namespace ipc::rigid