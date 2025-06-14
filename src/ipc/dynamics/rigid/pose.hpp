#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc::rigid {

using RotationVector = VectorMax3d;
using RotationMatrix = MatrixMax3d;

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

template <typename RotationType = RotationVector> struct Pose {
    // Position of the rigid body
    VectorMax3d position;
    // Rotation of the rigid body (rotation vector for 3D, angle for 1D)
    RotationType rotation;

    RotationMatrix rotation_matrix() const
    {
        if constexpr (std::is_same_v<RotationType, RotationVector>) {
            return rotation_vector_to_matrix(rotation);
        } else {
            return rotation;
        }
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
            // Eigen::AngleAxisd r =
            // Eigen::AngleAxisd(Eigen::Matrix3d(rotation)); return r.angle() *
            // r.axis();
            return rotation_matrix_to_vector(rotation);
        }
    }

    /// @brief Set the rotation matrix from a rotation vector or angle.
    /// @param theta The rotation vector (3D) or angle (2D).
    void set_rotation_vector(Eigen::ConstRef<RotationVector> theta)
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