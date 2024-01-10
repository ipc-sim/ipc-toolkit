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

class VirtualCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = Collision;
public:
    VirtualCollisions() = default;
    virtual ~VirtualCollisions() = default;

    virtual void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0,
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD);

    virtual void build(
        const Candidates& candidates,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0) = 0;

    /// @brief Computes the minimum distance between any non-adjacent elements.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @returns The minimum distance between any non-adjacent elements.
    double compute_minimum_distance(
        const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

    /// @brief Get the number of collisions.
    virtual size_t size() const = 0;

    /// @brief Get if the collision set are empty.
    virtual bool empty() const = 0;

    /// @brief Clear the collision set.
    virtual void clear() = 0;

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    virtual Collision& operator[](size_t i) = 0;

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    virtual const Collision& operator[](size_t i) const = 0;

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
    virtual void set_use_convergent_formulation(
        const bool use_convergent_formulation)
    {
        if (!empty()
            && use_convergent_formulation != m_use_convergent_formulation) {
            logger().warn(
                "Setting use_convergent_formulation after building collisions. "
                "Re-build collisions for this to have an effect.");
        }

        m_use_convergent_formulation = use_convergent_formulation;
    }

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
    virtual void set_are_shape_derivatives_enabled(const bool are_shape_derivatives_enabled) 
    {
        if (!empty()
            && are_shape_derivatives_enabled != m_are_shape_derivatives_enabled) {
            logger().warn(
                "Setting enable_shape_derivatives after building collisions. "
                "Re-build collisions for this to have an effect.");
        }

        m_are_shape_derivatives_enabled = are_shape_derivatives_enabled;
    }

    void set_edge_quadrature_type(const SurfaceQuadratureType type) { quad_type = type; }
    SurfaceQuadratureType get_edge_quadrature_type() const { return quad_type; }

    virtual std::vector<CandidateType> get_candidate_types(const int &dim) const = 0;

    // add edge-edge, edge-face, face-face neighbors to candidates
    virtual bool include_neighbor() const = 0;

protected:
    bool m_use_convergent_formulation = false;
    bool m_are_shape_derivatives_enabled = false;
    SurfaceQuadratureType quad_type = SurfaceQuadratureType::SinglePoint;
};

class Collisions : public VirtualCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = Collision;

public:
    Collisions() { }

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

    std::vector<CandidateType> get_candidate_types(const int &dim) const override;

    bool include_neighbor() const override { return false; }

    std::vector<VertexVertexCollision> vv_collisions;
    std::vector<EdgeVertexCollision> ev_collisions;
    std::vector<EdgeEdgeCollision> ee_collisions;
    std::vector<FaceVertexCollision> fv_collisions;
    std::vector<PlaneVertexCollision> pv_collisions;
};

} // namespace ipc
