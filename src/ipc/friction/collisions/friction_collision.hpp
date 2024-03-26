#pragma once

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/barrier/barrier.hpp>

#include <ipc/config.hpp>

namespace ipc {

class FrictionCollision : virtual public CollisionStencil {
protected:
    /// @brief Initialize the collision.
    /// @param collision Collision stencil.
    /// @param positions Collision stencil's vertex positions.
    /// @param barrier_potential Barrier potential used for normal force.
    /// @param barrier_stiffness Barrier potential stiffness.
    void init(
        const Collision& collision,
        const VectorMax12d& positions,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness);

public:
    virtual ~FrictionCollision() = default;

    /// @brief Get the dimension of the collision.
    int dim() const { return tangent_basis.rows(); }

    /// @brief Get the number of degrees of freedom for the collision.
    int ndof() const { return dim() * num_vertices(); };

    // -- Abstract methods -----------------------------------------------------

    /// @brief Compute the normal force magnitude.
    /// @param positions Collision stencil's vertex positions.
    /// @param dhat Barrier activation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Normal force magnitude.
    double compute_normal_force_magnitude(
        const VectorMax12d& positions,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the gradient of the normal force magnitude.
    /// @param positions Collision stencil's vertex positions.
    /// @param dhat Barrier activation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Gradient of the normal force magnitude wrt positions.
    VectorMax12d compute_normal_force_magnitude_gradient(
        const VectorMax12d& positions,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const double dmin = 0) const;

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
    /// @brief Collision force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
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
