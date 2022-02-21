#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

#include <stq/cpu/aabb.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <stq/gpu/aabb.cuh>
#endif

namespace ipc {

// A version of the BP that copies the meshes into the class rather than making
// the AABBs.
class CopyMeshBroadPhase : public BroadPhase {
protected:
    void copy_mesh(const Eigen::MatrixXi& E, const Eigen::MatrixXi& F);

    bool can_edge_vertex_collide(size_t ei, size_t vi) const override;
    bool can_edges_collide(size_t eai, size_t ebi) const override;
    bool can_face_vertex_collide(size_t fi, size_t vi) const override;
    bool can_edge_face_collide(size_t ei, size_t fi) const override;

    Eigen::MatrixXi edges;
    Eigen::MatrixXi faces;
};

class SweepAndTiniestQueue : public CopyMeshBroadPhase {
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
    long to_edge_id(long id) const;
    long to_face_id(long id) const;

    bool is_vertex(long id) const;
    bool is_edge(long id) const;
    bool is_face(long id) const;

    std::vector<stq::cpu::Aabb> boxes;
    std::vector<std::pair<int, int>> overlaps;
    long num_vertices;
};

#ifdef IPC_TOOLKIT_WITH_CUDA
class SweepAndTiniestQueueGPU : public CopyMeshBroadPhase {
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