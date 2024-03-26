#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision.hpp>
#include <ipc/collisions/vertex_vertex.hpp>
#include <ipc/collisions/edge_vertex.hpp>
#include <ipc/collisions/edge_edge.hpp>
#include <ipc/collisions/face_vertex.hpp>
#include <ipc/collisions/plane_vertex.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/candidates/candidates.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

class Collisions {
public:
    /// @brief The type of the collisions.
    using value_type = Collision;

public:
    Collisions() = default;

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

    /// @brief Get if the collision set should use the convergent formulation.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return If the collision set should use the convergent formulation.
    bool use_convergent_formulation() const
    {
        return m_use_convergent_formulation;
    }

    /// @brief Set if the collision set should use the convergent formulation.
    /// @warning This must be set before the collisions are built.
    /// @param use_convergent_formulation If the collision set should use the convergent formulation.
    void set_use_convergent_formulation(const bool use_convergent_formulation);

    /// @brief Get if the collision set are using the convergent formulation.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return If the collision set are using the convergent formulation.
    bool are_shape_derivatives_enabled() const
    {
        return m_are_shape_derivatives_enabled;
    }

    /// @brief Set if the collision set should enable shape derivative computation.
    /// @warning This must be set before the collisions are built.
    /// @param are_shape_derivatives_enabled If the collision set should enable shape derivative computation.
    void
    set_are_shape_derivatives_enabled(const bool are_shape_derivatives_enabled);

    std::string
    to_string(const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

public:
    std::vector<VertexVertexCollision> vv_collisions;
    std::vector<EdgeVertexCollision> ev_collisions;
    std::vector<EdgeEdgeCollision> ee_collisions;
    std::vector<FaceVertexCollision> fv_collisions;
    std::vector<PlaneVertexCollision> pv_collisions;

protected:
    bool m_use_convergent_formulation = false;
    bool m_are_shape_derivatives_enabled = false;
};

} // namespace ipc
