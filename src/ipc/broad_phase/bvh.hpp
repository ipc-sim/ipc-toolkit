#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

#include <SimpleBVH/BVH.hpp>

namespace ipc {

class BVH : public BroadPhase {
public:
    BVH() = default;

    /// @brief Build the broad phase for static collision detection.
    /// @param vertices Vertex positions
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double inflation_radius = 0) override;

    /// @brief Build the broad phase for continuous collision detection.
    /// @param vertices_t0 Starting vertices of the vertices.
    /// @param vertices_t1 Ending vertices of the vertices.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double inflation_radius = 0) override;

    /// @brief Clear any built data.
    void clear() override;

    /// @brief Find the candidate vertex-vertex collisions.
    /// @param[out] candidates The candidate vertex-vertex collisions.
    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-vertex collisions.
    /// @param[out] candidates The candidate edge-vertex collisions.
    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-edge collisions.
    /// @param[out] candidates The candidate edge-edge collisions.
    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override;

    /// @brief Find the candidate face-vertex collisions.
    /// @param[out] candidates The candidate face-vertex collisions.
    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-face intersections.
    /// @param[out] candidates The candidate edge-face intersections.
    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override;

    /// @brief Find the candidate face-face collisions.
    /// @param[out] candidates The candidate face-face collisions.
    void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const override;

protected:
    /// @brief Initialize a BVH from a set of boxes.
    /// @param[in] boxes Set of boxes to initialize the BVH with.
    /// @param[out] bvh The BVH to initialize.
    static void init_bvh(const std::vector<AABB>& boxes, SimpleBVH::BVH& bvh);

    /// @brief Detect candidate collisions between a BVH and a sets of boxes.
    /// @tparam Candidate Type of candidate collision.
    /// @tparam swap_order Whether to swap the order of box id with the BVH id when adding to the candidates.
    /// @tparam triangular Whether to consider (i, j) and (j, i) as the same.
    /// @param[in] boxes The boxes to detect collisions with.
    /// @param[in] bvh The BVH to detect collisions with.
    /// @param[in] can_collide Function to determine if two primitives can collide given their ids.
    /// @param[out] candidates The candidate collisions.
    template <
        typename Candidate,
        bool swap_order = false,
        bool triangular = false>
    static void detect_candidates(
        const std::vector<AABB>& boxes,
        const SimpleBVH::BVH& bvh,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates);

    /// @brief BVH containing the vertices.
    SimpleBVH::BVH vertex_bvh;
    /// @brief BVH containing the edges.
    SimpleBVH::BVH edge_bvh;
    /// @brief BVH containing the faces.
    SimpleBVH::BVH face_bvh;
};

} // namespace ipc