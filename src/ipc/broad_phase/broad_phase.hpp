#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/broad_phase/aabb.hpp>
#include <ipc/broad_phase/collision_candidate.hpp>

#include <Eigen/Core>

namespace ipc {

/// Enumeration of implemented broad phase methods.
enum class BroadPhaseMethod {
    BRUTE_FORCE = 0,
    HASH_GRID,
    SPATIAL_HASH,
    SWEEP_AND_TINIEST_QUEUE,
    SWEEP_AND_TINIEST_QUEUE_GPU, // Requires CUDA
    NUM_METHODS
};

class BroadPhase {
public:
    virtual ~BroadPhase() { clear(); }

    /// @brief Construct a registered broad phase object.
    /// @param method The broad phase method to use.
    /// @return The constructed broad phase object.
    static std::unique_ptr<BroadPhase>
    make_broad_phase(const BroadPhaseMethod method);

    /// @brief Build the broad phase for static collision detection.
    /// @param V0 Positions of the vertices.
    /// @param E Edges of the mesh.
    /// @param F Faces of the mesh.
    /// @param inflation_radius Radius of inflation around all elements.
    virtual void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0);

    /// @brief Build the broad phase for continuous collision detection.
    /// @param V0 Starting positions of the vertices.
    /// @param V1 Ending positions of the vertices.
    /// @param E Edges of the mesh.
    /// @param F Faces of the mesh.
    /// @param inflation_radius Radius of inflation around all elements.
    virtual void build(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0);

    /// @brief Clear any built data.
    virtual void clear();

    /// @brief Find the candidate edge-vertex collisisons.
    /// @param[out] candidates The candidate edge-vertex collisisons.
    virtual void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const = 0;

    /// @brief Find the candidate edge-edge collisions.
    /// @param[out] candidates The candidate edge-edge collisisons.
    virtual void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const = 0;

    /// @brief Find the candidate face-vertex collisions.
    /// @param[out] candidates The candidate face-vertex collisisons.
    virtual void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const = 0;

    /// @brief Find the candidate edge-face intersections.
    /// @param[out] candidates The candidate edge-face intersections.
    virtual void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const = 0;

    /// @brief Detect all collision candidates needed for a given dimensional simulation.
    /// @param dim The dimension of the simulation (i.e., 2 or 3).
    /// @param candidates The detected collision candidates.
    virtual void
    detect_collision_candidates(int dim, Candidates& candidates) const;

    /// @brief Function for determining if two vertices can collide.
    std::function<bool(size_t, size_t)> can_vertices_collide =
        [](size_t, size_t) { return true; };

protected:
    virtual bool can_edge_vertex_collide(size_t ei, size_t vi) const;
    virtual bool can_edges_collide(size_t eai, size_t ebi) const;
    virtual bool can_face_vertex_collide(size_t fi, size_t vi) const;
    virtual bool can_edge_face_collide(size_t ei, size_t fi) const;

    std::vector<AABB> vertex_boxes;
    std::vector<AABB> edge_boxes;
    std::vector<AABB> face_boxes;
};

/// @brief Construct a set of discrete collision detection candidates.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface Vertex positions at start as rows of a matrix.
/// @param[out] candidates The constructed candidate set as output.
/// @param[in] inflation_radius Amount to inflate the bounding boxes.
/// @param[in] method Broad phase method to use.
void construct_collision_candidates(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    Candidates& candidates,
    double inflation_radius = 0,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID);

/// @brief Construct a set of continuous collision detection candidates.
/// @note Assumes the trajectory is linear.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V0 Surface vertex positions at start as rows of a matrix.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @param[out] candidates The constructed candidate set as output.
/// @param[in] inflation_radius Amount to inflate the bounding boxes.
/// @param[in] method Broad phase method to use.
void construct_collision_candidates(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    Candidates& candidates,
    double inflation_radius = 0,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID);

} // namespace ipc
