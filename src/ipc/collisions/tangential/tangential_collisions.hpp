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

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double normal_stiffness,
        double mu)
    {
        this->build(
            mesh, vertices, collisions, normal_potential, normal_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu));
    }

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
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
