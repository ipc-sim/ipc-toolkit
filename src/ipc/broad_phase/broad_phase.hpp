#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/broad_phase/aabb.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/edge_face.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_face.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/candidates/vertex_vertex.hpp>

#include <Eigen/Core>

namespace ipc {

class Candidates; // Forward declaration

/// @brief Base class for broad phase collision detection methods.
class BroadPhase {
public:
    virtual ~BroadPhase() { clear(); }

    /// @brief Get the name of the broad phase method.
    /// @return The name of the broad phase method.
    virtual std::string name() const = 0;

    /// @brief Build the broad phase for static collision detection.
    /// @param vertices Vertex positions
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    virtual void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0);

    /// @brief Build the broad phase for continuous collision detection.
    /// @param vertices_t0 Starting vertices of the vertices.
    /// @param vertices_t1 Ending vertices of the vertices.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    virtual void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0);

    /// @brief Clear any built data.
    virtual void clear();

    /// @brief Detect all collision candidates needed for a given dimensional simulation.
    /// @param dim The dimension of the simulation (i.e., 2 or 3).
    /// @param candidates The detected collision candidates.
    void detect_collision_candidates(int dim, Candidates& candidates) const;

    /// @brief Find the candidate vertex-vertex collisions.
    /// @param[out] candidates The candidate vertex-vertex collisions.
    virtual void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const = 0;

    /// @brief Find the candidate edge-vertex collisions.
    /// @param[out] candidates The candidate edge-vertex collisions.
    virtual void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const = 0;

    /// @brief Find the candidate edge-edge collisions.
    /// @param[out] candidates The candidate edge-edge collisions.
    virtual void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const = 0;

    /// @brief Find the candidate face-vertex collisions.
    /// @param[out] candidates The candidate face-vertex collisions.
    virtual void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const = 0;

    /// @brief Find the candidate edge-face intersections.
    /// @param[out] candidates The candidate edge-face intersections.
    virtual void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const = 0;

    /// @brief Find the candidate face-face collisions.
    /// @param[out] candidates The candidate face-face collisions.
    virtual void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const = 0;

    /// @brief Function for determining if two vertices can collide.
    std::function<bool(size_t, size_t)> can_vertices_collide =
        default_can_vertices_collide;

protected:
    virtual bool can_edge_vertex_collide(size_t ei, size_t vi) const;
    virtual bool can_edges_collide(size_t eai, size_t ebi) const;
    virtual bool can_face_vertex_collide(size_t fi, size_t vi) const;
    virtual bool can_edge_face_collide(size_t ei, size_t fi) const;
    virtual bool can_faces_collide(size_t fai, size_t fbi) const;

    static bool default_can_vertices_collide(size_t, size_t) { return true; }

    std::vector<AABB> vertex_boxes;
    std::vector<AABB> edge_boxes;
    std::vector<AABB> face_boxes;
};

} // namespace ipc
