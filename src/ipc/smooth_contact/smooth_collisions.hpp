#pragma once

#include <ipc/collisions/collisions.hpp>
#include "edge_vertex.hpp"
#include "edge_edge.hpp"
#include "face_vertex.hpp"

namespace ipc {

template <int dim>
class SmoothCollisions : public VirtualCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = Collision;

public:
    SmoothCollisions() = default;
    SmoothCollisions(bool _use_adaptive_eps)
    : use_adaptive_eps(_use_adaptive_eps)
    {
    }

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
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD) override;

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
        const double dmin = 0) override;

    // ------------------------------------------------------------------------

    /// @brief Get the number of collisions.
    size_t size() const override;

    /// @brief Get if the collision set are empty.
    bool empty() const override;

    /// @brief Clear the collision set.
    void clear() override;

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    Collision& operator[](size_t i) override;

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const Collision& operator[](size_t i) const override;

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

    void set_use_convergent_formulation(
        const bool use_convergent_formulation) override
    {
        logger().error("Smooth contact formulation doesn't have convergent version!");
    }

    void set_are_shape_derivatives_enabled(const bool are_shape_derivatives_enabled) override
    {
        logger().error("Smooth contact formulation doesn't have shape derivatives implemented!");
    }

    std::vector<CandidateType> get_candidate_types(const int &_dim) const override;

public:
    // std::vector<SmoothVertexVertexCollision> vv_collisions;
    std::vector<SmoothEdgeVertexCollision> ev_collisions;
    std::vector<SmoothEdgeEdgeCollision<dim>> ee_collisions;
    std::vector<SmoothFaceVertexCollision> fv_collisions;
    // std::vector<SmoothPlaneVertexCollision> pv_collisions;

    const bool use_adaptive_eps = false;
};

} // namespace ipc
