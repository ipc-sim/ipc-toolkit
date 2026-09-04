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

/// @brief Jacobian of the rotation matrix w.r.t. the rotation parameters,
/// dispatching on the dimension.
///
/// For a rotation parameterized by @p theta (a scalar angle in 2D or a rotation
/// vector in 3D), this returns d vec(R)/d theta as a (dim²)×rot_ndof matrix:
/// 2D → 4×1 (so(2) exp map, single generator), 3D → 9×3 (Rodrigues).
///
/// @param theta The rotation parameters (size 1 in 2D, size 3 in 3D)
/// @return The Jacobian, sized 4×1 (2D) or 9×3 (3D)
MatrixMax<double, 9, 3>
rotation_to_matrix_jacobian(Eigen::ConstRef<VectorMax3d> theta);

/// @brief Hessian of the rotation matrix w.r.t. the rotation parameters,
/// dispatching on the dimension.
///
/// Returns d² vec(R)/d theta² flattened as a (dim²)×(rot_ndof²) matrix:
/// 2D → 4×1 (d²R/dθ² = −R), 3D → 9×9.
///
/// @param theta The rotation parameters (size 1 in 2D, size 3 in 3D)
/// @return The Hessian, sized 4×1 (2D) or 9×9 (3D)
MatrixMax<double, 9, 9>
rotation_to_matrix_hessian(Eigen::ConstRef<VectorMax3d> theta);

// ----------------------------------------------------------------------------

struct Pose {
    /// @brief Position of the rigid body
    VectorMax3d position;
    /// @brief Rotation of the rigid body (rotation vector for 3D, angle for 1D)
    VectorMax3d rotation;

    /// @brief Default constructor.
    Pose() = default;

    /// @brief Construct a pose from position and rotation vectors.
    /// @param _position The position vector
    /// @param _rotation The rotation vector
    Pose(
        Eigen::ConstRef<VectorMax3d> _position,
        Eigen::ConstRef<VectorMax3d> _rotation)
        : position(_position)
        , rotation(_rotation)
    {
    }

    /// @brief Construct a pose from a concatenated vector.
    /// @param x The concatenated position and rotation vector
    explicit Pose(Eigen::ConstRef<VectorMax6d> x)
    {
        assert(x.size() == 6 || x.size() == 3);
        if (x.size() == 3) {
            position = x.head<2>();
            rotation = x.tail<1>();
        } else {
            position = x.head<3>();
            rotation = x.tail<3>();
        }
    }

    /// @brief Construct an identity pose with zero position and rotation.
    /// @param dim The dimension of the pose (2 or 3)
    /// @return Identity pose with zero position and rotation
    static Pose Identity(const int dim) // NOLINT(readability-identifier-naming)
    {
        assert(dim == 2 || dim == 3);
        // Rotation DOF: 2D uses a single angle, 3D uses a rotation vector.
        const int rot_ndof = dim == 2 ? 1 : 3;
        return Pose(VectorMax3d::Zero(dim), VectorMax3d::Zero(rot_ndof));
    }

    int ndof() const { return position.size() + rotation.size(); }

    /// @brief Construct a rotation matrix from the rotation vector.
    /// @return The rotation matrix corresponding to the rotation vector
    MatrixMax3d rotation_matrix() const
    {
        assert(rotation.size() == 1 || rotation.size() == 3);
        if (rotation.size() == 1) {
            // For 2D, set the rotation matrix directly
            MatrixMax3d R(2, 2);
            R << std::cos(rotation(0)), -std::sin(rotation(0)),
                std::sin(rotation(0)), std::cos(rotation(0));
            return R;
        } else {
            // For 3D, convert the rotation vector to a rotation matrix
            return rotation_vector_to_matrix(rotation);
        }
    }

    /// @brief Construct a quaternion from the rotation vector.
    /// @return The quaternion corresponding to the rotation vector
    Eigen::Quaternion<double> quaternion() const
    {
        assert(rotation.size() == 3);
        double angle = rotation.norm();
        if (angle == 0.0) {
            return Eigen::Quaternion<double>::Identity();
        }
        Eigen::Vector3d axis = rotation / angle;
        return Eigen::Quaternion<double>(Eigen::AngleAxis<double>(angle, axis));
    }

    /// @brief Clamp the rotation angle to the range [-π, π] to prevent numerical issues.
    /// This is important for optimization to ensure that the rotation does not
    /// grow unbounded and cause numerical instability.
    void clamp_rotation_to_pi()
    {
        double angle = rotation.norm();

        // 1. Small angle approximation / Zero check
        if (angle < std::numeric_limits<double>::epsilon()) {
            return;
        }

        // 2. Extract unit axis
        VectorMax3d axis = rotation / angle;

        // 3. Wrap angle to [-π, π] using remainder
        angle = std::remainder(angle, 2.0 * EIGEN_PI);

        // 4. Canonicalize PI: Ensure -π < angle <= π
        // This prevents the representation from jumping between PI and -PI
        if (angle <= -EIGEN_PI) {
            angle = EIGEN_PI;
            axis = -axis; // Flip axis to keep angle positive
        }

        rotation = angle * axis;
    }

    /// @brief Transform vertices from local to world coordinates using the pose.
    /// @param V The vertices to transform (each row is a vertex)
    /// @return The transformed vertices
    Eigen::MatrixXd transform_vertices(Eigen::ConstRef<Eigen::MatrixXd> V) const
    {
        // Compute: R(θ) V + p
        // transpose because x is row-ordered
        return (V * rotation_matrix().transpose()).rowwise()
            + position.transpose();
    }

    /// @brief Compute the Jacobian of the transformed vertices with respect to the pose.
    /// @param V The vertices to transform (each row is a vertex)
    /// @return The Jacobian matrix of size (num_vertices * dim) x ndof
    Eigen::MatrixXd
    transform_vertices_jacobian(Eigen::ConstRef<Eigen::MatrixXd> V) const;

    /// @brief Compute the Hessian of the transformed vertices with respect to the pose.
    /// @param V The vertices to transform (each row is a vertex)
    /// @return The Hessian matrix of size (num_vertices * dim) x (ndof * ndof)
    Eigen::MatrixXd
    transform_vertices_hessian(Eigen::ConstRef<Eigen::MatrixXd> V) const;

    /// @brief Compute the inverse of the pose.
    /// @return The inverse pose
    Pose inverse() const
    {
        Pose inv;
        const MatrixMax3d R_inv = rotation_matrix().transpose();
        inv.position = -(R_inv * position);
        if (position.size() == 2) {
            // Negate angle directly
            inv.rotation = -rotation;
        } else {
            // Convert inverse rotation to a vector
            inv.rotation = rotation_matrix_to_vector(R_inv);
        }
        return inv;
    }

    /// @brief Combine two poses.
    /// @param a The first pose
    /// @param b The second pose
    /// @return The combined pose
    friend Pose operator*(const Pose& a, const Pose& b)
    {
        Pose c;

        const MatrixMax3d Ra = a.rotation_matrix();

        c.position = Ra * b.position + a.position;

        if (c.position.size() == 2) {
            // Combine angles directly
            c.rotation = a.rotation + b.rotation;
        } else {
            // Combine rotation matrices multiplicatively
            c.rotation = rotation_matrix_to_vector(Ra * b.rotation_matrix());
        }

        return c;
    }

    /// @brief Check if two poses are equal.
    /// @param a The first pose
    /// @param b The second pose
    /// @return True if the poses are equal, false otherwise
    friend bool operator==(const Pose& a, const Pose& b)
    {
        return a.position == b.position && a.rotation == b.rotation;
    }

    static std::vector<Pose>
    to_poses(Eigen::ConstRef<Eigen::VectorXd> x, const int dim)
    {
        const int pose_ndof = dim == 2 ? 3 : 6;
        assert(x.size() % pose_ndof == 0);
        std::vector<Pose> poses(x.size() / pose_ndof);
        for (size_t i = 0; i < poses.size(); ++i) {
            poses[i] = Pose(x.segment(i * pose_ndof, pose_ndof));
        }
        return poses;
    }

    static Eigen::VectorXd from_poses(const std::vector<Pose>& poses)
    {
        const int pose_ndof = poses[0].ndof();
        Eigen::VectorXd x(poses.size() * pose_ndof);
        for (size_t i = 0; i < poses.size(); ++i) {
            assert(poses[i].ndof() == pose_ndof);
            x.segment(i * pose_ndof, pose_ndof) << poses[i].position,
                poses[i].rotation;
        }
        return x;
    }
};

} // namespace ipc::rigid