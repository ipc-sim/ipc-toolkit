#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>
#include <ipc/collisions/normal/edge_vertex.hpp>
#include <ipc/collisions/normal/face_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/collisions/normal/plane_vertex.hpp>
#include <ipc/collisions/normal/vertex_vertex.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

class NormalCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = NormalCollision;

    /// @brief How to construct the collision set.
    /// This choice is described in the "Offset Geometric Contact" [Chen et al.
    /// 2025] using the notation: \f$\mathcal{F}\f$ and \f$\mathcal{E}\f$ for
    /// the face and edge sets, respectively.
    enum class CollisionSetType : uint8_t {
        /// @brief IPC set type, which uses the original formulation described in [Li et al. 2020].
        IPC,
        /// @brief Improved max approximation set type, which uses the improved max approximation formulation described in [Li et al. 2023].
        IMPROVED_MAX_APPROX,
        /// @brief Offset Geometric Contact set type, which uses the formulation described in [Chen et al. 2025].
        OGC,
    };

public:
    NormalCollisions() = default;

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param dmin Minimum distance.
    /// @param broad_phase Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const double dhat,
        const double dmin = 0,
        const std::shared_ptr<BroadPhase>& broad_phase =
            make_default_broad_phase());

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param  dmin  Minimum distance.
    void build(
        const Candidates& candidates,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const double dhat,
        const double dmin = 0);

    // ------------------------------------------------------------------------

    /// @brief Computes the minimum distance between any non-adjacent elements.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @returns The minimum distance between any non-adjacent elements.
    double compute_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

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
    NormalCollision& operator[](size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const NormalCollision& operator[](size_t i) const;

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
    [[deprecated(
        "use_improved_max_approximator() is deprecated. Use collision_set_type() instead.")]]
    bool use_improved_max_approximator() const
    {
        return m_collision_set_type == CollisionSetType::IMPROVED_MAX_APPROX;
    }

    /// @brief Set if the collision set should use the improved max approximator.
    /// @warning This must be set before the collisions are built.
    /// @param use_improved_max_approximator If the collision set should use the improved max approximator.
    [[deprecated(
        "set_use_improved_max_approximator() is deprecated. Use set_collision_set_type() instead.")]]
    void
    set_use_improved_max_approximator(const bool use_improved_max_approximator)
    {
        // Default to IPC if not using improved max approximator
        set_collision_set_type(
            use_improved_max_approximator
                ? CollisionSetType::IMPROVED_MAX_APPROX
                : CollisionSetType::IPC);
    }

    /// @brief Get the set type of the collision set.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return The set type of the collision set.
    CollisionSetType collision_set_type() const { return m_collision_set_type; }

    /// @brief Set the set type of the collision set.
    /// @param set_type The set type of the collision set.
    void set_collision_set_type(const CollisionSetType type);

    /// @brief Get if the collision set are using the convergent formulation.
    /// @note If not empty, this is the current value not necessarily the value used to build the collisions.
    /// @return If the collision set are using the convergent formulation.
    bool enable_shape_derivatives() const { return m_enable_shape_derivatives; }

    /// @brief Set if the collision set should enable shape derivative computation.
    /// @warning This must be set before the collisions are built.
    /// @param enable_shape_derivatives If the collision set should enable shape derivative computation.
    void set_enable_shape_derivatives(const bool enable_shape_derivatives);

    std::string to_string(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

public:
    /// @brief Vertex-vertex normal collisions.
    std::vector<VertexVertexNormalCollision> vv_collisions;
    /// @brief Edge-vertex normal collisions.
    std::vector<EdgeVertexNormalCollision> ev_collisions;
    /// @brief Edge-edge normal collisions.
    std::vector<EdgeEdgeNormalCollision> ee_collisions;
    /// @brief Face-vertex normal collisions.
    std::vector<FaceVertexNormalCollision> fv_collisions;
    /// @brief Plane-vertex normal collisions.
    std::vector<PlaneVertexNormalCollision> pv_collisions;

protected:
    CollisionSetType m_collision_set_type = CollisionSetType::IPC;
    bool m_use_area_weighting = false;
    bool m_enable_shape_derivatives = false;
};

} // namespace ipc
