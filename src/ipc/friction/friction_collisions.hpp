#pragma once

#include <ipc/friction/collisions/friction_collision.hpp>
#include <ipc/friction/collisions/vertex_vertex.hpp>
#include <ipc/friction/collisions/edge_vertex.hpp>
#include <ipc/friction/collisions/edge_edge.hpp>
#include <ipc/friction/collisions/face_vertex.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/smooth_contact/smooth_collisions.hpp>

namespace ipc {

class FrictionCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = FrictionCollision;

public:
    FrictionCollisions() = default;

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Collisions& collisions,
        const BarrierPotential& barrier_potential,
        double barrier_stiffness,
        double mu)
    {
        this->build(
            mesh, vertices, collisions, barrier_potential, barrier_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu));
    }

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Collisions& collisions,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        const Eigen::VectorXd& mus,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    template <int dim>
    void build_for_smooth_contact(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const SmoothCollisions<dim>& collisions,
        const ParameterType &params,
        const double barrier_stiffness,
        double mu)
    {
        this->build_for_smooth_contact(
            mesh, vertices, collisions, params, barrier_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu));
    }

    template <int dim>
    void build_for_smooth_contact(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const SmoothCollisions<dim>& collisions,
        const ParameterType &params,
        const double barrier_stiffness,
        const Eigen::VectorXd& mus,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    // ------------------------------------------------------------------------

    /// @brief Get the number of friction collisions.
    size_t size() const;

    /// @brief Get if the friction collisions are empty.
    bool empty() const;

    /// @brief Clear the friction collisions.
    void clear();

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    FrictionCollision& operator[](const size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const FrictionCollision& operator[](const size_t i) const;

    static double default_blend_mu(double mu0, double mu1)
    {
        // return mu0 * mu1;
        // return std::min(mu0, mu1);
        // return std::max(mu0, mu1);
        return (mu0 + mu1) / 2;
    }

public:
    std::vector<VertexVertexFrictionCollision> vv_collisions;
    std::vector<EdgeVertexFrictionCollision> ev_collisions;
    std::vector<EdgeEdgeFrictionCollision> ee_collisions;
    std::vector<FaceVertexFrictionCollision> fv_collisions;

    double barrier_stiffness_;
};

} // namespace ipc
