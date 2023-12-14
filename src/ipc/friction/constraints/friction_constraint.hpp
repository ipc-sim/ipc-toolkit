#pragma once

#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <ipc/config.hpp>

namespace ipc {

class FrictionConstraint : virtual public CollisionStencil {
protected:
    /// @brief Initialize the constraint.
    /// @param positions Vertex positions(rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param dmin Minimum distance
    void init(
        const Eigen::MatrixXd& positions,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double dmin);

public:
    virtual ~FrictionConstraint() { }

    /// @brief Get the dimension of the constraint.
    int dim() const { return tangent_basis.rows(); }

    /// @brief Get the number of degrees of freedom for the constraint.
    int ndof() const { return dim() * num_vertices(); };

    // -------------------------------------------------------------------------
    // Abstract methods
    // -------------------------------------------------------------------------

    /// @brief Compute the normal force magnitude.
    /// @param positions Constraint's vertex positions.
    /// @param dhat Barrier activiation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Normal force magnitude.
    double compute_normal_force_magnitude(
        const VectorMax12d& positions,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the gradient of the normal force magnitude.
    /// @param positions Constraint's vertex positions.
    /// @param dhat Barrier activiation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Gradient of the normal force magnitude wrt positions.
    VectorMax12d compute_normal_force_magnitude_gradient(
        const VectorMax12d& positions,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the tangent basis of the constraint.
    /// @param positions Constraint's vertex positions.
    /// @return Tangent basis of the constraint.
    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& positions) const = 0;

    /// @brief Compute the Jacobian of the tangent basis of the constraint.
    /// @param positions Constraint's vertex positions.
    /// @return Jacobian of the tangent basis of the constraint.
    virtual MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& positions) const = 0;

    /// @brief Compute the barycentric coordinates of the closest point.
    /// @param positions Constraint's vertex positions.
    /// @return Barycentric coordinates of the closest point.
    virtual VectorMax2d
    compute_closest_point(const VectorMax12d& positions) const = 0;

    /// @brief Compute the Jacobian of the barycentric coordinates of the closest point.
    /// @param positions Constraint's vertex positions.
    /// @return Jacobian of the barycentric coordinates of the closest point.
    virtual MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& positions) const = 0;

    /// @brief Compute the relative velocity of the constraint.
    /// @param velocities Constraint's vertex velocities.
    /// @return Relative velocity of the constraint.
    virtual VectorMax3d
    relative_velocity(const VectorMax12d& velocities) const = 0;

    /// @brief Construct the premultiplier matrix for the relative velocity.
    /// @note Uses the cached closest point.
    /// @return A matrix M such that `relative_velocity = M * velocities`.
    virtual MatrixMax<double, 3, 12> relative_velocity_matrix() const
    {
        return relative_velocity_matrix(closest_point);
    }

    /// @brief Construct the premultiplier matrix for the relative velocity.
    /// @param closest_point Barycentric coordinates of the closest point.
    /// @return A matrix M such that `relative_velocity = M * velocities`.
    virtual MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const = 0;

    /// @brief Construct the Jacobian of the relative velocity premultiplier wrt the closest points.
    /// @param closest_point Barycentric coordinates of the closest point.
    /// @return Jacobian of the relative velocity premultiplier wrt the closest points.
    virtual MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const = 0;

public:
    /// @brief Contact force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
    double mu;

    /// @brief Weight
    double weight = 1;

    /// @brief Gradient of weight with respect to all DOF
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Barycentric coordinates of the closest point(s)
    VectorMax2d closest_point;

    /// @brief Tangent basis of the contact (max size 3×2)
    MatrixMax<double, 3, 2> tangent_basis;
};

} // namespace ipc
