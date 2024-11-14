#pragma once

#include <ipc/potentials/normal_potential.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class TangentialCollision : virtual public CollisionStencil {
protected:
    /// @brief Initialize the collision.
    /// @param collision NormalCollision stencil.
    /// @param positions Collision stencil's vertex positions.
    /// @param normal_potential Normal potential used for normal force.
    /// @param normal_stiffness Normal potential stiffness.
    void init(
        const NormalCollision& collision,
        const VectorMax12d& positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);

public:
    virtual ~TangentialCollision() = default;

    /// @brief Get the dimension of the collision.
    int dim() const { return tangent_basis.rows(); }

    /// @brief Get the number of degrees of freedom for the collision.
    int ndof() const { return dim() * num_vertices(); };

    // -- Abstract methods -----------------------------------------------------

    /// @brief Compute the tangent basis of the collision.
    /// @param positions Collision stencil's vertex positions.
    /// @return Tangent basis of the collision.
    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& positions) const = 0;

    /// @brief Compute the Jacobian of the tangent basis of the collision.
    /// @param positions Collision stencil's vertex positions.
    /// @return Jacobian of the tangent basis of the collision.
    virtual MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& positions) const = 0;

    /// @brief Compute the barycentric coordinates of the closest point.
    /// @param positions Collision stencil's vertex positions.
    /// @return Barycentric coordinates of the closest point.
    virtual VectorMax2d
    compute_closest_point(const VectorMax12d& positions) const = 0;

    /// @brief Compute the Jacobian of the barycentric coordinates of the closest point.
    /// @param positions Collision stencil's vertex positions.
    /// @return Jacobian of the barycentric coordinates of the closest point.
    virtual MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& positions) const = 0;

    /// @brief Compute the relative velocity of the collision.
    /// @param positions Collision stencil's vertex velocities.
    /// @return Relative velocity of the collision.
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
    /// @brief Normal force magnitude
    double normal_force_magnitude;

    /// @brief Ratio between normal and tangential forces (e.g., friction coefficient)
    double mu;

    /// @brief Weight
    double weight = 1;

    /// @brief Gradient of weight with respect to all DOF
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Barycentric coordinates of the closest point(s)
    VectorMax2d closest_point;

    /// @brief Tangent basis of the collision (max size 3×2)
    MatrixMax<double, 3, 2> tangent_basis;
};

} // namespace ipc
