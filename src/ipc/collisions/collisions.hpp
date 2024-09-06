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

    /// @brief Get if the collision set should use area weighting.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return If the collision set should use area weighting.
    bool use_area_weighting() const { return m_use_area_weighting; }

    /// @brief Set if the collision set should use area weighting.
    /// @warning This must be set before the collisions are built.
    /// @param use_area_weighting If the collision set should use area weighting.
    void set_use_area_weighting(const bool use_area_weighting);

    /// @brief Get if the collision set should use the improved max approximator.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return If the collision set should use the improved max approximator.
    bool use_improved_max_approximator() const
    {
        return m_use_improved_max_approximator;
    }

    /// @brief Set if the collision set should use the improved max approximator.
    /// @warning This must be set before the collisions are built.
    /// @param use_improved_max_approximator If the collision set should use the improved max approximator.
    void
    set_use_improved_max_approximator(const bool use_improved_max_approximator);

    /// @brief Get if the collision set are using the convergent formulation.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return If the collision set are using the convergent formulation.
    bool enable_shape_derivatives() const { return m_enable_shape_derivatives; }

    /// @brief Set if the collision set should enable shape derivative computation.
    /// @warning This must be set before the collisions are built.
    /// @param enable_shape_derivatives If the collision set should enable shape derivative computation.
    void set_enable_shape_derivatives(const bool enable_shape_derivatives);

    std::string
    to_string(const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

public:
    std::vector<VertexVertexCollision> vv_collisions;
    std::vector<EdgeVertexCollision> ev_collisions;
    std::vector<EdgeEdgeCollision> ee_collisions;
    std::vector<FaceVertexCollision> fv_collisions;
    std::vector<PlaneVertexCollision> pv_collisions;

protected:
    bool m_use_area_weighting = false;
    bool m_use_improved_max_approximator = false;
    bool m_enable_shape_derivatives = false;
};

} // namespace ipc
