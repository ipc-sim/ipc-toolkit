#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/edge_edge.hpp>
#include <ipc/collisions/tangential/edge_vertex.hpp>
#include <ipc/collisions/tangential/face_vertex.hpp>
#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/collisions/tangential/vertex_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

class TangentialCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = TangentialCollision;

public:
    TangentialCollisions() = default;

    /// @brief Build the tangential collisions.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param normal_potential The normal potential.
    /// @param normal_stiffness Stiffness of the normal potential.
    /// @param mu The coefficient of friction.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double normal_stiffness,
        double mu)
    {
        this->build(
            mesh, vertices, collisions, normal_potential, normal_stiffness, mu,
            mu);
    }

    /// @brief Build the tangential collisions.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param normal_potential The normal potential.
    /// @param normal_stiffness Stiffness of the normal potential.
    /// @param mu_s The static friction coefficient.
    /// @param mu_k The kinetic friction coefficient.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double normal_stiffness,
        double mu_s,
        double mu_k)
    {
        this->build(
            mesh, vertices, collisions, normal_potential, normal_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu_s),
            Eigen::VectorXd::Constant(vertices.rows(), mu_k));
    }

    /// @brief Build the tangential collisions.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param normal_potential The normal potential.
    /// @param normal_stiffness Stiffness of the normal potential.
    /// @param mu_k The kinetic friction coefficient per vertex.
    /// @param mu_s The static friction coefficient per vertex.
    /// @param blend_mu Function to blend vertex-based coefficients of friction.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
        Eigen::ConstRef<Eigen::VectorXd> mu_k,
        Eigen::ConstRef<Eigen::VectorXd> mu_s,
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
    TangentialCollision& operator[](const size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const TangentialCollision& operator[](const size_t i) const;

    static double default_blend_mu(double mu0, double mu1)
    {
        // return mu0 * mu1;
        // return std::min(mu0, mu1);
        // return std::max(mu0, mu1);
        return (mu0 + mu1) / 2;
    }

public:
    /// @brief Vertex-vertex tangential collisions.
    std::vector<VertexVertexTangentialCollision> vv_collisions;
    /// @brief Edge-vertex tangential collisions.
    std::vector<EdgeVertexTangentialCollision> ev_collisions;
    /// @brief Edge-edge tangential collisions.
    std::vector<EdgeEdgeTangentialCollision> ee_collisions;
    /// @brief Face-vertex tangential collisions.
    std::vector<FaceVertexTangentialCollision> fv_collisions;
};

} // namespace ipc
