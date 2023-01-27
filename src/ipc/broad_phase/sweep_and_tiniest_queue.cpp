#include "sweep_and_tiniest_queue.hpp"

#include <ipc/config.hpp>

#include <stq/cpu/io.hpp>
#include <stq/cpu/sweep.hpp>
#ifdef IPC_TOOLKIT_WITH_CUDA
#include <ccdgpu/helper.cuh>
#endif

namespace ipc {

void SweepAndTiniestQueue::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& _edges,
    const Eigen::MatrixXi& _faces,
    const double inflation_radius)
{
    build(vertices, vertices, _edges, _faces, inflation_radius);
}

void SweepAndTiniestQueue::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& _edges,
    const Eigen::MatrixXi& _faces,
    const double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(_edges, _faces);
    num_vertices = vertices_t0.rows();
    stq::cpu::constructBoxes(
        vertices_t0, vertices_t1, edges, faces, boxes, inflation_radius);
    int n = boxes.size();
    stq::cpu::sort_along_xaxis(boxes);
    stq::cpu::run_sweep_cpu(boxes, n, overlaps);
}

void SweepAndTiniestQueue::clear()
{
    BroadPhase::clear();
    num_vertices = 0;
    boxes.clear();
    overlaps.clear();
}

void SweepAndTiniestQueue::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

void SweepAndTiniestQueue::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

void SweepAndTiniestQueue::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    for (const auto& [id1, id2] : overlaps) {
        if (is_edge(id1) && is_edge(id2)
            && can_edges_collide(to_edge_id(id1), to_edge_id(id2))) { // EE
            candidates.emplace_back(to_edge_id(id1), to_edge_id(id2));
        }
    }
}

void SweepAndTiniestQueue::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    for (const auto& [id1, id2] : overlaps) {
        if (is_face(id1) && is_vertex(id2)
            && can_face_vertex_collide(to_face_id(id1), id2)) { // FV
            candidates.emplace_back(to_face_id(id1), id2);
        } else if (
            is_face(id2) && is_vertex(id1)
            && can_face_vertex_collide(to_face_id(id2), id1)) { // VF
            candidates.emplace_back(to_face_id(id2), id1);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

long SweepAndTiniestQueue::to_edge_id(long id) const
{
    assert(id >= num_vertices);
    assert(id < num_vertices + this->edges.rows());
    return id - num_vertices;
}

long SweepAndTiniestQueue::to_face_id(long id) const
{
    assert(id >= num_vertices + this->edges.rows());
    assert(id < num_vertices + this->edges.rows() + this->faces.rows());
    return id - num_vertices - this->edges.rows();
}

bool SweepAndTiniestQueue::is_vertex(long id) const
{
    return id >= 0 && id < num_vertices;
}

bool SweepAndTiniestQueue::is_edge(long id) const
{
    return id >= num_vertices && id < num_vertices + this->edges.rows();
}

bool SweepAndTiniestQueue::is_face(long id) const
{
    return id >= num_vertices + this->edges.rows()
        && id < num_vertices + this->edges.rows() + this->faces.rows();
}

// ============================================================================

#ifdef IPC_TOOLKIT_WITH_CUDA
void SweepAndTiniestQueueGPU::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& _edges,
    const Eigen::MatrixXi& _faces,
    const double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(_edges, _faces);
    ccd::gpu::construct_static_collision_candidates(
        vertices, edges, faces, overlaps, boxes, inflation_radius);
}

void SweepAndTiniestQueueGPU::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& _edges,
    const Eigen::MatrixXi& _faces,
    const double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(_edges, _faces);
    ccd::gpu::construct_continuous_collision_candidates(
        vertices_t0, vertices_t1, edges, faces, overlaps, boxes,
        inflation_radius);
}

void SweepAndTiniestQueueGPU::clear()
{
    BroadPhase::clear();
    overlaps.clear();
    boxes.clear();
}

void SweepAndTiniestQueueGPU::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

void SweepAndTiniestQueueGPU::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    // 2D STQ is not implemented!
    throw std::runtime_error("Not implemented!");
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

void SweepAndTiniestQueueGPU::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    using namespace stq::gpu;
    for (const std::pair<int, int>& overlap : overlaps) {
        const Aabb& boxA = boxes[overlap.first];
        const Aabb& boxB = boxes[overlap.second];
        if (is_edge(boxA) && is_edge(boxB)
            && can_edges_collide(boxA.ref_id, boxB.ref_id)) { // EE
            candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        }
    }
}

void SweepAndTiniestQueueGPU::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    using namespace stq::gpu;
    for (const std::pair<int, int>& overlap : overlaps) {
        const Aabb& boxA = boxes[overlap.first];
        const Aabb& boxB = boxes[overlap.second];
        if (is_face(boxA) && is_vertex(boxB)
            && can_face_vertex_collide(boxA.ref_id, boxB.ref_id)) { // FV
            candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        } else if (
            is_face(boxB) && is_vertex(boxA)
            && can_face_vertex_collide(boxB.ref_id, boxA.ref_id)) { // VF
            candidates.emplace_back(boxB.ref_id, boxA.ref_id);
        }
    }
}

void SweepAndTiniestQueueGPU::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}
#endif

// ============================================================================

void CopyMeshBroadPhase::copy_mesh(
    const Eigen::MatrixXi& p_edges, const Eigen::MatrixXi& p_faces)
{
    edges = p_edges;
    faces = p_faces;
}

bool CopyMeshBroadPhase::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const long e0i = edges(ei, 0), e1i = edges(ei, 1);

    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool CopyMeshBroadPhase::can_edges_collide(size_t eai, size_t ebi) const
{
    const long ea0i = edges(eai, 0), ea1i = edges(eai, 1);
    const long eb0i = edges(ebi, 0), eb1i = edges(ebi, 1);

    const bool share_endpoint =
        ea0i == eb0i || ea0i == eb1i || ea1i == eb0i || ea1i == eb1i;

    return !share_endpoint
        && (can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
            || can_vertices_collide(ea1i, eb0i)
            || can_vertices_collide(ea1i, eb1i));
}

bool CopyMeshBroadPhase::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const long f0i = faces(fi, 0), f1i = faces(fi, 1), f2i = faces(fi, 2);

    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool CopyMeshBroadPhase::can_edge_face_collide(size_t ei, size_t fi) const
{
    const long e0i = edges(ei, 0), e1i = edges(ei, 1);
    const long f0i = faces(fi, 0), f1i = faces(fi, 1), f2i = faces(fi, 2);

    const bool share_endpoint = e0i == f0i || e0i == f1i || e0i == f2i
        || e1i == f0i || e1i == f1i || e1i == f2i;

    return !share_endpoint
        && (can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
            || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
            || can_vertices_collide(e1i, f1i)
            || can_vertices_collide(e1i, f2i));
}

} // namespace ipc