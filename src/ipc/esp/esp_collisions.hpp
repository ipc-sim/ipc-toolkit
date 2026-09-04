#pragma once

#include "adaptive_support.hpp"
#include "collisions/esp_collision.hpp"
#include "collisions/esp_collision_dict.hpp"

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>

#include <map>

namespace ipc {
class ESPCollisions {
public:
    ESPCollisions() = default;
    virtual ~ESPCollisions() = default;

    /// @brief Compute per-vertex adaptive dhat values. The returned object can
    ///        be passed to build() to avoid recomputing it on every rebuild.
    static std::unique_ptr<AdaptiveSupport> compute_adaptive_dhat(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters& params);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters params,
        BroadPhase* broad_phase = nullptr);

    /// @brief Build using a pre-computed AdaptiveSupport (copied internally;
    ///        pass nullptr to build without adaptive dhat).
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters params,
        const AdaptiveSupport* adaptive,
        BroadPhase* broad_phase = nullptr);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    void build(
        const Candidates& _candidates,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters params);

    /// @brief Build from candidates using a pre-computed AdaptiveSupport (copied internally).
    void build(
        const Candidates& _candidates,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters params,
        const AdaptiveSupport* adaptive);

    // ------------------------------------------------------------------------

    /// @brief Get the number of collisions.
    size_t size() const;

    /// @brief Get if the collision set are empty.
    bool empty() const;

    /// @brief Clear the collision set.
    void clear();

    /// @brief Compute minimum distance between all contact candidates
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @return Squared minimum distance
    double compute_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    /// @brief Convert contact pairs to string
    std::string to_string(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters& params) const;

    /// @brief Number of contact candidates
    int n_candidates() const { return m_candidates.size(); }

    /// @brief Count occurrences of each edge id across all edge_edge_collisions keys.
    /// @return A map from occurrence count to the number of edge ids with that count.
    std::map<size_t, size_t> edge_id_count_distribution() const;

    /// @brief Get per-edge collision counts as a vector.
    /// @param num_edges Total number of edges in the collision mesh.
    /// @return A vector of size num_edges where entry i counts how many
    ///         edge-edge collision pairs involve edge i (as first element).
    Eigen::VectorXd edge_collision_counts(size_t num_edges) const;

public:
    /// @brief per-vertex adaptive dhat interpolated on edges/faces
    std::unique_ptr<AdaptiveSupport> adaptive_dhat = nullptr;

    /// @brief Collision candidates
    Candidates m_candidates;

    /// @brief collision sets for 3D quadrature

    // vertex_collisions[vi] provides the contact set for vertex vi
    unordered_map<index_t, std::unique_ptr<ESPCollisionDict<PointType::VERTEX>>>
        vertex_collisions;
    // edge_edge_collisions[(ei, ej)] provides the contact set for the closest
    // point on ei, between edge ei and ej.
    unordered_map<
        std::pair<index_t, index_t>,
        std::unique_ptr<ESPCollisionDict<PointType::EDGE>>>
        edge_edge_collisions;
    // face_collisions[fi][qi] provides the contact set for quadrature point qi
    // of face fi
    unordered_map<
        index_t,
        std::vector<std::unique_ptr<ESPCollisionDict<PointType::FACE>>>>
        face_collisions;

    /// @brief collision sets for 2D quadrature
    // edge_collisions_2d[ei][qi] provides the contact set for Gauss-Lobatto QP
    // qi on edge ei
    unordered_map<
        index_t,
        std::vector<std::unique_ptr<ESPCollisionDict<PointType::EDGE, 2>>>>
        edge_collisions_2d;

    /// @brief Total number of collision pairs counted across all quadrature build functions
    size_t num_quadrature_collision_pairs = 0;
};
} // namespace ipc