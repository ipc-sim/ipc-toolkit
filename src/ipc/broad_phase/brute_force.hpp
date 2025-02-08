#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

namespace ipc {

class BruteForce : public BroadPhase {
public:
    BruteForce() = default;

    /// @brief Get the name of the broad phase method.
    /// @return The name of the broad phase method.
    std::string name() const override { return "BruteForce"; }

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

private:
    /// @brief Detect candidates for collisions between two sets of boxes.
    /// @tparam Candidate Type of the candidate.
    /// @tparam triangular Whether to consider (i, j) and (j, i) as the same.
    /// @param[in] boxes0 First set of boxes.
    /// @param[in] boxes1 Second set of boxes.
    /// @param[in] can_collide Function to determine if two primitives can collide given their ids.
    /// @param[out] candidates The candidate collisions.
    template <typename Candidate, bool triangular = false>
    void detect_candidates(
        const std::vector<AABB>& boxes0,
        const std::vector<AABB>& boxes1,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates) const;
};

} // namespace ipc
