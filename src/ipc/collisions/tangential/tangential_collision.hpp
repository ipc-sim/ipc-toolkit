#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/potentials/normal_potential.hpp>

namespace ipc {

class TangentialCollision {
public:
    TangentialCollision() = default;

    void init(
        const NormalCollision& collision,
        Eigen::ConstRef<VectorMax12d> positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);

    virtual ~TangentialCollision() = default;

    // =========================================================================
    // Methods that must be implemented by each collision type

    /// @brief Get the dimension of the space.
    virtual int dim() const = 0;

    /// @brief Get the number of degrees of freedom.
    virtual int ndof() const = 0;

    /// @brief Get the number of vertices in the stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the vertex IDs of the contact stencil.
    /// @param edges The edge matrix of mesh (#E × 2).
    /// @param faces The face matrix of mesh (#F × 3).
    /// @return The vertex IDs involved in the contact stencil.
    virtual std::array<long, 4>
    vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges, 
        Eigen::ConstRef<Eigen::MatrixXi> faces) const = 0;

    /// @brief Extract the stencil's DOF from a global DOF vector
    /// @param dof The global DOF vector (#V × dim).
    /// @param edges The edge matrix of mesh (#E × 2).
    /// @param faces The face matrix of mesh (#F × 3).
    /// @return The stencil's DOF (≤ 4 × dim).
    virtual VectorMax12d dof(
        Eigen::ConstRef<Eigen::MatrixXd> dof,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const = 0;

    /// @brief Compute the distance.
    /// @param positions The stencil's vertex positions.
    /// @return The distance.
    virtual double compute_distance(Eigen::ConstRef<VectorMax12d> positions)
        const = 0;

    /// @brief Compute the distance gradient wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The distance gradient.
    virtual VectorMax12d
    compute_distance_gradient(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    // =========================================================================
    // Main methods for tangential collision

    /// @brief Compute the tangent basis.
    /// @param positions The stencil's vertex positions.
    /// @return The tangent basis matrix ∈ ℝ^{3×2}.
    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the Jacobian of the tangent basis wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The Jacobian of the tangent basis matrix ∈ ℝ^{(3×dim)×2}.
    virtual MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(Eigen::ConstRef<VectorMax12d> positions)
        const = 0;

    /// @brief Compute the closest points.
    /// @param positions The stencil's vertex positions.
    /// @return The closest points coordinates ∈ ℝ^{0…2}.
    virtual VectorMax2d
    compute_closest_point(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the Jacobian of the closest points wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The Jacobian of the closest point.
    virtual MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(Eigen::ConstRef<VectorMax12d> positions)
        const = 0;

    /// @brief Compute the relative velocity of the two objects in contact.
    /// @param velocities The stencil's vertex velocities.
    /// @return The relative velocity ∈ ℝ^3.
    virtual VectorMax3d
    relative_velocity(Eigen::ConstRef<VectorMax12d> velocities) const = 0;

    /// @brief Compute the relative velocity matrix.
    /// @return The relative velocity matrix ∈ ℝ^{3×(≤4×dim)}.
    MatrixMax<double, 3, 12> relative_velocity_matrix() const;

    /// @brief Compute the relative velocity matrix.
    /// @param closest_point The closest point coordinates ∈ ℝ^{0…2}.
    /// @return The relative velocity matrix ∈ ℝ^{3×(≤4×dim)}.
    virtual MatrixMax<double, 3, 12> relative_velocity_matrix(
        Eigen::ConstRef<VectorMax2d> closest_point) const = 0;

    /// @brief Compute the Jacobian of the relative velocity matrix wrt the closest point.
    /// @param closest_point The closest point coordinates ∈ ℝ^{0…2}.
    /// @return The Jacobian of the relative velocity matrix ∈ ℝ^{6×(≤4×dim)}.
    virtual MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        Eigen::ConstRef<VectorMax2d> closest_point) const = 0;

    // =========================================================================
    // Variables

    /// @brief The contact basis expressed in the world frame.
    MatrixMax<double, 3, 2> tangent_basis;

    /// @brief The coordinates of the closest points.
    VectorMax2d closest_point;

    /// @brief The magnitude of the contact force.
    double normal_force_magnitude = 0;

    /// @brief The term's weight (e.g., collision area)
    double weight = 1;

    /// @brief The gradient of the term's weight wrt the rest positions.
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Global friction coefficient.
    double mu = 0;

    /// @brief Static friction coefficient.
    double s_mu = -1;

    /// @brief Kinetic friction coefficient.
    double k_mu = -1;

    /// @brief The material ID for the first object in this collision
    int material_id1 = NO_MATERIAL_ID;

    /// @brief The material ID for the second object in this collision
    int material_id2 = NO_MATERIAL_ID;

    /// @brief Check if material IDs are being used for this collision
    /// @return true if both material IDs are valid, false otherwise
    bool has_material_ids() const { return material_id1 != NO_MATERIAL_ID && material_id2 != NO_MATERIAL_ID; }
};

} // namespace ipc
