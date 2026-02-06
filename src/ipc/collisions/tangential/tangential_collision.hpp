#pragma once

#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/potentials/normal_potential.hpp>
#include <ipc/smooth_contact/collisions/smooth_collision.hpp>
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
        Eigen::ConstRef<VectorMax12d> positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);

    /// @brief Initialize the collision.
    /// @param collision NormalCollision stencil.
    /// @param positions Collision stencil's vertex positions.
    /// @param normal_force Magnitude of the normal force.
    void init(
        const NormalCollision& collision,
        Eigen::ConstRef<VectorMax12d> positions,
        const double normal_force);

public:
    virtual ~TangentialCollision() = default;

    /// @brief Get the dimension of the collision.
    int dim() const { return tangent_basis.rows(); }

    /// @brief Get the number of degrees of freedom for the collision.
    int ndof() const { return dim() * num_vertices(); }

    // -- Abstract methods -----------------------------------------------------

    /// @brief Compute the tangent basis of the collision.
    /// @param positions Collision stencil's vertex positions.
    /// @return Tangent basis of the collision.
    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the Jacobian of the tangent basis of the collision.
    /// @param positions Collision stencil's vertex positions.
    /// @return Jacobian of the tangent basis of the collision.
    virtual MatrixMax<double, 36, 2> compute_tangent_basis_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the barycentric coordinates of the closest point.
    /// @param positions Collision stencil's vertex positions.
    /// @return Barycentric coordinates of the closest point.
    virtual VectorMax2d
    compute_closest_point(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the Jacobian of the barycentric coordinates of the closest point.
    /// @param positions Collision stencil's vertex positions.
    /// @return Jacobian of the barycentric coordinates of the closest point.
    virtual MatrixMax<double, 2, 12> compute_closest_point_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the relative velocity of the collision.
    /// @param velocities Collision stencil's vertex velocities.
    /// @return Relative velocity of the collision.
    virtual VectorMax3d
    relative_velocity(Eigen::ConstRef<VectorMax12d> velocities) const = 0;

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
    virtual MatrixMax<double, 3, 12> relative_velocity_matrix(
        Eigen::ConstRef<VectorMax2d> closest_point) const = 0;

    /// @brief Construct the Jacobian of the relative velocity premultiplier wrt the closest points.
    /// @param closest_point Barycentric coordinates of the closest point.
    /// @return Jacobian of the relative velocity premultiplier wrt the closest points.
    virtual MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        Eigen::ConstRef<VectorMax2d> closest_point) const = 0;

public:
    /// @brief Normal force magnitude
    double normal_force_magnitude = 0;

    /// @brief SmoothCollision instance to compute normal force magnitude and its derivatives
    std::shared_ptr<SmoothCollision> smooth_collision;

    /// @brief Ratio between normal and static tangential forces (e.g., friction coefficient)
    double mu_s = 0;

    /// @brief Ratio between normal and kinetic tangential forces (e.g., friction coefficient)
    double mu_k = 0;

    /// @brief Anisotropic static friction coefficients (2D, one per tangent direction).
    /// @note Zero vector → scalar mu_s (backward compatible). Elliptical model;
    ///       see ipc::smooth_mu and Erleben et al., CGF 2019,
    ///       DOI 10.1111/cgf.13885.
    Eigen::Vector2d mu_s_aniso = Eigen::Vector2d::Zero();

    /// @brief Anisotropic kinetic friction coefficients (2D, one per tangent direction).
    /// @note Zero vector → scalar mu_k (backward compatible). Elliptical model;
    ///       see ipc::smooth_mu and Erleben et al., CGF 2019,
    ///       DOI 10.1111/cgf.13885.
    Eigen::Vector2d mu_k_aniso = Eigen::Vector2d::Zero();

    /// @brief Weight
    double weight = 1;

    /// @brief Tangential anisotropy scaling in the collision's tangent basis.
    /// @note Default (1,1) preserves current isotropic behavior.
    ///       Requires a_i > 0. Values scale tau before friction evaluation.
    ///       Used with mu_s_aniso/mu_k_aniso by the elliptical model in
    ///       ipc::smooth_mu.
    Eigen::Vector2d mu_aniso = Eigen::Vector2d::Ones();

    /// @brief Gradient of weight with respect to all DOF
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Barycentric coordinates of the closest point(s)
    VectorMax2d closest_point;

    /// @brief Tangent basis of the collision (max size 3×2)
    MatrixMax<double, 3, 2> tangent_basis;
};

} // namespace ipc
