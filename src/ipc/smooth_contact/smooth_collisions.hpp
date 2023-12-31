#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision.hpp>
#include "edge_vertex.hpp"
#include "edge_edge.hpp"
// #include <ipc/collisions/face_vertex.hpp>
// #include <ipc/collisions/plane_vertex.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/candidates/candidates.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

class SmoothCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = Collision;

public:
    SmoothCollisions() { }

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param dmin Minimum distance.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0,
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param  dmin  Minimum distance.
    void build(
        const Candidates& candidates,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0);

    // ------------------------------------------------------------------------

    /// @brief Computes the minimum distance between any non-adjacent elements.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @returns The minimum distance between any non-adjacent elements.
    double compute_minimum_distance(
        const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

    // ------------------------------------------------------------------------

    /// @brief Get the number of collisions.
    size_t size() const;

    /// @brief Get if the collision set are empty.
    bool empty() const;

    /// @brief Clear the collision set.
    void clear();

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    Collision& operator[](size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const Collision& operator[](size_t i) const;

    /// @brief Get if the collision at i is a vertex-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is a vertex-vertex collision.
    bool is_vertex_vertex(size_t i) const;

    /// @brief Get if the collision at i is an edge-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is an edge-vertex collision.
    bool is_edge_vertex(size_t i) const;

    /// @brief Get if the collision at i is an edge-edge collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is an edge-edge collision.
    bool is_edge_edge(size_t i) const;

    /// @brief Get if the collision at i is an face-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is an face-vertex collision.
    bool is_face_vertex(size_t i) const;

    /// @brief Get if the collision at i is an plane-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is an plane-vertex collision.
    bool is_plane_vertex(size_t i) const;

    std::string
    to_string(const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

public:
    // std::vector<SmoothVertexVertexCollision> vv_collisions;
    std::vector<SmoothEdgeVertexCollision> ev_collisions;
    std::vector<SmoothEdgeEdgeCollision> ee_collisions;
    // std::vector<SmoothFaceVertexCollision> fv_collisions;
    // std::vector<SmoothPlaneVertexCollision> pv_collisions;
};

} // namespace ipc
