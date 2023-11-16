#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

namespace ipc {

class BruteForce : public BroadPhase {
public:
    /// @brief Find the candidate vertex-vertex collisions.
    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-vertex collisions.
    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-edge collisions.
    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override;

    /// @brief Find the candidate face-vertex collisions.
    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-face intersections.
    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override;

private:
    template <typename Candidate, bool triangular = false>
    void detect_candidates(
        const std::vector<AABB>& boxes0,
        const std::vector<AABB>& boxes1,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates) const;
};

} // namespace ipc
