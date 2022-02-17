#include <ipc/broad_phase/sweep.hpp>

#include <tbb/parallel_sort.h>

#include <stq/cpu/io.hpp>
#include <stq/cpu/sweep.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <ccdgpu/helper.cuh>
#endif

namespace ipc {

void SweepAndTiniestQueue::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double inflation_radius)
{
    build(V, V, E, F, inflation_radius);
}

void SweepAndTiniestQueue::build(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double inflation_radius)
{
    num_vertices = V0.rows();
    num_edges = E.rows();
    num_faces = F.rows();
    stq::cpu::constructBoxes(V0, V1, E, F, boxes, inflation_radius);
    stq::cpu::run_sweep_cpu(boxes, overlaps);
}

void SweepAndTiniestQueue::clear()
{
    BroadPhase::clear();
    num_vertices = 0;
    num_edges = 0;
    num_faces = 0;
    boxes.clear();
    overlaps.clear();
}

/// @brief Find the candidate edge-vertex collisisons.
void SweepAndTiniestQueue::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    throw "Not implemented!";
}

/// @brief Find the candidate edge-edge collisions.
void SweepAndTiniestQueue::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    using namespace stq::cpu;
    for (const std::pair<int, int>& overlap : overlaps) {
        const Aabb& boxA = boxes[overlap.first];
        const Aabb& boxB = boxes[overlap.second];
        if (is_edge(boxA.vertexIds) && is_edge(boxB.vertexIds)) {
            // && can_edges_collide(boxA.id, boxB.id)) { // EE
            candidates.emplace_back(
                boxA.id - num_vertices, boxB.id - num_vertices);
        }
    }
}

/// @brief Find the candidate face-vertex collisions.
void SweepAndTiniestQueue::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    using namespace stq::cpu;
    for (const std::pair<int, int>& overlap : overlaps) {
        const Aabb& boxA = boxes[overlap.first];
        const Aabb& boxB = boxes[overlap.second];
        if (is_face(boxA.vertexIds) && is_vertex(boxB.vertexIds)) {
            // && can_face_vertex_collide(boxA.id, boxB.id)) { // FV
            candidates.emplace_back(
                boxA.id - num_vertices - num_edges, boxB.id);
        } else if (is_face(boxB.vertexIds) && is_vertex(boxA.vertexIds)) {
            // && can_face_vertex_collide(boxB.id, boxA.id)) { // VF
            candidates.emplace_back(
                boxB.id - num_vertices - num_edges, boxA.id);
        }
    }
}

/// @brief Find the candidate edge-face intersections.
void SweepAndTiniestQueue::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    throw "Not implemented!";
}

////////////////////////////////////////////////////////////////////////////////

#ifdef IPC_TOOLKIT_WITH_CUDA
void SweepAndTiniestQueueGPU::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double inflation_radius)
{
    ccd::gpu::construct_static_collision_candidates(
        V, E, F, overlaps, boxes, inflation_radius);
}

void SweepAndTiniestQueueGPU::build(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double inflation_radius)
{
    ccd::gpu::construct_continuous_collision_candidates(
        V0, V1, E, F, overlaps, boxes, inflation_radius);
}

void SweepAndTiniestQueueGPU::clear()
{
    BroadPhase::clear();
    overlaps.clear();
    boxes.clear();
}

// Find the candidate edge-vertex collisisons.
void SweepAndTiniestQueueGPU::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    // 2D SQT is not implemented!
    throw "Not implemented!";
    // using namespace stq::gpu;
    // for (const std::pair<int, int>& overlap : overlaps) {
    //     const Aabb& boxA = boxes[overlap.first];
    //     const Aabb& boxB = boxes[overlap.second];
    //     if (is_edge(boxA) && is_vertex(boxB)
    //         && can_edge_vertex_collide(boxA.ref_id, boxB.ref_id)) { // EV
    //         candidates.emplace_back(boxA.ref_id, boxB.ref_id);
    //     } else if (
    //         is_edge(boxB) && is_vertex(boxA)
    //         && can_edge_vertex_collide(boxB.ref_id, boxA.ref_id)) { // VE
    //         candidates.emplace_back(boxB.ref_id, boxA.ref_id);
    //     }
    // }
}

// Find the candidate edge-edge collisions.
void SweepAndTiniestQueueGPU::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    using namespace stq::gpu;
    for (const std::pair<int, int>& overlap : overlaps) {
        const Aabb& boxA = boxes[overlap.first];
        const Aabb& boxB = boxes[overlap.second];
        if (is_edge(boxA) && is_edge(boxB)) {
            // && can_edges_collide(boxA.ref_id, boxB.ref_id)) { // EE
            candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        }
    }
}

// Find the candidate face-vertex collisions.
void SweepAndTiniestQueueGPU::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    using namespace stq::gpu;
    for (const std::pair<int, int>& overlap : overlaps) {
        const Aabb& boxA = boxes[overlap.first];
        const Aabb& boxB = boxes[overlap.second];
        if (is_face(boxA) && is_vertex(boxB)) {
            // && can_face_vertex_collide(boxA.ref_id, boxB.ref_id)) { // FV
            candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        } else if (is_face(boxB) && is_vertex(boxA)) {
            // && can_face_vertex_collide(boxB.ref_id, boxA.ref_id)) { // VF
            candidates.emplace_back(boxB.ref_id, boxA.ref_id);
        }
    }
}

// Find the candidate edge-face intersections.
void SweepAndTiniestQueueGPU::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    throw "Not implemented!";
}
#endif

} // namespace ipc