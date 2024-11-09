#pragma once

#include <ipc/friction/collisions/friction_collision.hpp>
#include <ipc/friction/collisions/vertex_vertex.hpp>
#include <ipc/friction/collisions/edge_vertex.hpp>
#include <ipc/friction/collisions/edge_edge.hpp>
#include <ipc/friction/collisions/face_vertex.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <optional>
#include <functional>
#include <map>
#include <tuple>

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
        const std::function<double(double, double, std::optional<BlendType>)>&
            blend_mu = default_blend_mu);

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Collisions& collisions,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness,
        double mu,
        const double static_mu,
        const double kinetic_mu,
        const std::map<std::tuple<int, int>, std::pair<double, double>>& pairwise_friction,
        const std::function<double(double, double, std::optional<BlendType>)>& blend_mu = default_blend_mu);

    // ------------------------------------------------------------------------

    /// @brief Get the number of friction collisions.
    size_t size() const;

    /// @brief Get if the friction collisions are empty.
    bool empty() const;

    /// @brief Clear the friction collisions.
    void clear();

    /// @brief Get a reference to the collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    FrictionCollision& operator[](const size_t i);

    /// @brief Get a const reference to the collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const FrictionCollision& operator[](const size_t i) const;

    /// @brief Default blend function for friction coefficients.
    /// Uses the specified blending type or defaults to average blending.
    static double default_blend_mu(
        double mu0, double mu1, std::optional<BlendType> type = BlendType::AVG)
    {
        return blend_mu(mu0, mu1, type);
    }

    /// @brief Retrieve the friction coefficients for a pair of vertices.
    /// @param id1 The ID of the first vertex.
    /// @param id2 The ID of the second vertex.
    /// @param pairwise_friction The map containing pairwise friction coefficients.
    /// @param static_mu The default static friction coefficient.
    /// @param kinetic_mu The default kinetic friction coefficient.
    /// @return The static and kinetic friction coefficients.
    std::pair<double, double> retrieve_friction_coefficients(
        int id1, int id2,
        const std::map<std::tuple<int, int>, std::pair<double, double>>& pairwise_friction_map,
        std::optional<double> static_mu = std::nullopt,
        std::optional<double> kinetic_mu = std::nullopt) const;

    /// @brief Retrieves the static and kinetic friction coefficients for a pair of vertices.
    /// @param mesh Reference to the collision mesh to access material IDs.
    /// @param vertex1 Index of the first vertex.
    /// @param vertex2 Index of the second vertex.
    /// @param pairwise_friction Map containing pairwise friction coefficients.
    /// @param default_static_mu Default static friction coefficient.
    /// @param default_kinetic_mu Default kinetic friction coefficient.
    /// @return std::pair<double, double> Pair containing {static_mu, kinetic_mu}.
    std::pair<double, double> get_pairwise_friction_coefficients(
        const CollisionMesh& mesh,
        int vertex1,
        int vertex2,
        const std::map<std::tuple<int, int>, std::pair<double, double>>&
            pairwise_friction,
        double default_static_mu,
        double default_kinetic_mu) const;

public:
    std::vector<VertexVertexFrictionCollision> vv_collisions;
    std::vector<EdgeVertexFrictionCollision> ev_collisions;
    std::vector<EdgeEdgeFrictionCollision> ee_collisions;
    std::vector<FaceVertexFrictionCollision> fv_collisions;
};
} // namespace ipc
