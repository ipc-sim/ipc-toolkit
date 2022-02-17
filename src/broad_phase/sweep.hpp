#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

#include <stq/cpu/aabb.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <stq/gpu/aabb.cuh>
#endif

namespace ipc {

class SweepAndTiniestQueue : public BroadPhase {
public:
    void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0) override;

    void build(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0) override;

    void clear() override;

    /// @brief Find the candidate edge-vertex collisisons.
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

protected:
    std::vector<stq::cpu::Aabb> boxes;
    std::vector<std::pair<int, int>> overlaps;
    size_t num_vertices;
    size_t num_edges;
    size_t num_faces;
};

#ifdef IPC_TOOLKIT_WITH_CUDA
class SweepAndTiniestQueueGPU : public BroadPhase {
public:
    void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0) override;

    void build(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0) override;

    void clear() override;

    /// @brief Find the candidate edge-vertex collisisons.
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
    std::vector<stq::gpu::Aabb> boxes;
    std::vector<std::pair<int, int>> overlaps;
};
#endif

} // namespace ipc