#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

#include <SimpleBVH/BVH.hpp>

namespace ipc {

class BVH : public BroadPhase {
public:
    /// @brief Build the broad phase for static collision detection.
    /// @param V0 Positions of the vertices.
    /// @param E Edges of the mesh.
    /// @param F Faces of the mesh.
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double inflation_radius = 0) override;

    /// @brief Build the broad phase for continuous collision detection.
    /// @param V0 Starting positions of the vertices.
    /// @param V1 Ending positions of the vertices.
    /// @param E Edges of the mesh.
    /// @param F Faces of the mesh.
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
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

protected:
    static void init_bvh(const std::vector<AABB>& boxes, SimpleBVH::BVH& bvh);

    template <
        typename Candidate,
        bool swap_order = false,
        bool triangular = false>
    static void detect_candidates(
        const std::vector<AABB>& boxes,
        const SimpleBVH::BVH& bvh,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates);

    SimpleBVH::BVH vertex_bvh;
    SimpleBVH::BVH edge_bvh;
    SimpleBVH::BVH face_bvh;
};

} // namespace ipc