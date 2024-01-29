#include "sweep_and_tiniest_queue.hpp"

#include <ipc/config.hpp>

#include <scalable_ccd/stq/sweep.hpp>
#ifdef IPC_TOOLKIT_WITH_CUDA
#include <scalable_ccd/cuda/tight_inclusion/helper.cuh>
#endif

namespace ipc {

void SweepAndTiniestQueue::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    build(vertices, vertices, edges, faces, inflation_radius);
}

void SweepAndTiniestQueue::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    clear();

    scalable_ccd::build_vertex_boxes(
        vertices_t0, vertices_t1, vertex_boxes, inflation_radius);
    scalable_ccd::build_edge_boxes(vertex_boxes, edges, edge_boxes);
    scalable_ccd::build_face_boxes(vertex_boxes, faces, face_boxes);
}

void SweepAndTiniestQueue::clear()
{
    BroadPhase::clear();
    vertex_boxes.clear();
    edge_boxes.clear();
    face_boxes.clear();
}

void SweepAndTiniestQueue::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    std::vector<std::pair<int, int>> overlaps;
    scalable_ccd::sort_and_sweep(vertex_boxes, vv_sort_axis, overlaps);

    for (const auto& [vai, vbi] : overlaps) {
        if (can_vertices_collide(vai, vbi)) {
            candidates.emplace_back(vai, vbi);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    std::vector<std::pair<int, int>> overlaps;
    scalable_ccd::sort_and_sweep(
        edge_boxes, vertex_boxes, ev_sort_axis, overlaps);

    for (const auto& [ei, vi] : overlaps) {
        if (can_edge_vertex_collide(ei, vi)) {
            candidates.emplace_back(ei, vi);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    std::vector<std::pair<int, int>> overlaps;
    scalable_ccd::sort_and_sweep(edge_boxes, ee_sort_axis, overlaps);

    for (const auto& [eai, ebi] : overlaps) {
        if (can_edges_collide(eai, ebi)) {
            candidates.emplace_back(eai, ebi);
        }
    }
}

void SweepAndTiniestQueue::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    std::vector<std::pair<int, int>> overlaps;
    scalable_ccd::sort_and_sweep(
        face_boxes, vertex_boxes, fv_sort_axis, overlaps);

    for (const auto& [fi, vi] : overlaps) {
        if (can_face_vertex_collide(fi, vi)) {
            candidates.emplace_back(fi, vi);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    std::vector<std::pair<int, int>> overlaps;
    scalable_ccd::sort_and_sweep(
        edge_boxes, face_boxes, ef_sort_axis, overlaps);

    for (const auto& [ei, fi] : overlaps) {
        if (can_edge_face_collide(ei, fi)) {
            candidates.emplace_back(ei, fi);
        }
    }
}

void SweepAndTiniestQueue::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    std::vector<std::pair<int, int>> overlaps;
    scalable_ccd::sort_and_sweep(face_boxes, ff_sort_axis, overlaps);

    for (const auto& [fai, fbi] : overlaps) {
        if (can_faces_collide(fai, fbi)) {
            candidates.emplace_back(fai, fbi);
        }
    }
}

// ----------------------------------------------------------------------------

bool SweepAndTiniestQueue::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const auto& [e0i, e1i, _] = edge_boxes[ei].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(vi != e0i && vi != e1i);

    return can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i);
}

bool SweepAndTiniestQueue::can_edges_collide(size_t eai, size_t ebi) const
{
    const auto& [ea0i, ea1i, _] = edge_boxes[eai].vertex_ids;
    const auto& [eb0i, eb1i, __] = edge_boxes[ebi].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(ea0i != eb0i && ea0i != eb1i && ea1i != eb0i && ea1i != eb1i);

    return can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
        || can_vertices_collide(ea1i, eb0i) || can_vertices_collide(ea1i, eb1i);
}

bool SweepAndTiniestQueue::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const auto& [f0i, f1i, f2i] = face_boxes[fi].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(vi != f0i && vi != f1i && vi != f2i);

    return can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
        || can_vertices_collide(vi, f2i);
}

bool SweepAndTiniestQueue::can_edge_face_collide(size_t ei, size_t fi) const
{
    const auto& [e0i, e1i, _] = edge_boxes[ei].vertex_ids;
    const auto& [f0i, f1i, f2i] = face_boxes[fi].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(
        e0i != f0i && e0i != f1i && e0i != f2i && e1i != f0i && e1i != f1i
        && e1i != f2i);

    return can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
        || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
        || can_vertices_collide(e1i, f1i) || can_vertices_collide(e1i, f2i);
}

bool SweepAndTiniestQueue::can_faces_collide(size_t fai, size_t fbi) const
{
    const auto& [fa0i, fa1i, fa2i] = face_boxes[fai].vertex_ids;
    const auto& [fb0i, fb1i, fb2i] = face_boxes[fbi].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(
        fa0i != fb0i && fa0i != fb1i && fa0i != fb2i && fa1i != fb0i
        && fa1i != fb1i && fa1i != fb2i && fa2i != fb0i && fa2i != fb1i
        && fa2i != fb2i);

    return can_vertices_collide(fa0i, fb0i) || can_vertices_collide(fa0i, fb1i)
        || can_vertices_collide(fa0i, fb2i) || can_vertices_collide(fa1i, fb0i)
        || can_vertices_collide(fa1i, fb1i) || can_vertices_collide(fa1i, fb2i)
        || can_vertices_collide(fa2i, fb0i) || can_vertices_collide(fa2i, fb1i)
        || can_vertices_collide(fa2i, fb2i);
}

// ============================================================================

#ifdef IPC_TOOLKIT_WITH_CUDA
void SweepAndTiniestQueueGPU::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& _edges,
    const Eigen::MatrixXi& _faces,
    const double inflation_radius)
{
    assert(_edges.size() == 0 || _edges.cols() == 2);
    assert(_faces.size() == 0 || _faces.cols() == 3);

    clear();

    edges = _edges;
    faces = _faces;

    scalable_ccd::cuda::construct_static_collision_candidates(
        vertices, edges, faces, overlaps, boxes, inflation_radius);
}

void SweepAndTiniestQueueGPU::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& _edges,
    const Eigen::MatrixXi& _faces,
    const double inflation_radius)
{
    assert(_edges.size() == 0 || _edges.cols() == 2);
    assert(_faces.size() == 0 || _faces.cols() == 3);

    clear();

    edges = _edges;
    faces = _faces;

    scalable_ccd::cuda::construct_continuous_collision_candidates(
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
    throw std::runtime_error(
        "SweepAndTiniestQueueGPU::detect_vertex_vertex_candidates not implemented!");
}

void SweepAndTiniestQueueGPU::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    // 2D STQ is not implemented!
    throw std::runtime_error(
        "SweepAndTiniestQueueGPU::detect_edge_vertex_candidates not implemented!");
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
    using namespace scalable_ccd::cuda::stq;
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
    using namespace scalable_ccd::cuda::stq;
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
    throw std::runtime_error(
        "SweepAndTiniestQueueGPU::detect_edge_face_candidates not implemented!");
}

void SweepAndTiniestQueueGPU::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    throw std::runtime_error(
        "SweepAndTiniestQueueGPU::detect_face_face_candidates not implemented!");
}

// ----------------------------------------------------------------------------

bool SweepAndTiniestQueueGPU::can_edge_vertex_collide(
    size_t ei, size_t vi) const
{
    const long e0i = edges(ei, 0), e1i = edges(ei, 1);

    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool SweepAndTiniestQueueGPU::can_edges_collide(size_t eai, size_t ebi) const
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

bool SweepAndTiniestQueueGPU::can_face_vertex_collide(
    size_t fi, size_t vi) const
{
    const long f0i = faces(fi, 0), f1i = faces(fi, 1), f2i = faces(fi, 2);

    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool SweepAndTiniestQueueGPU::can_edge_face_collide(size_t ei, size_t fi) const
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

bool SweepAndTiniestQueueGPU::can_faces_collide(size_t fai, size_t fbi) const
{
    const long fa0i = faces(fai, 0), fa1i = faces(fai, 1), fa2i = faces(fai, 2);
    const long fb0i = faces(fbi, 0), fb1i = faces(fbi, 1), fb2i = faces(fbi, 2);

    const bool share_endpoint = fa0i == fb0i || fa0i == fb1i || fa0i == fb2i
        || fa1i == fb0i || fa1i == fb1i || fa1i == fb2i || fa2i == fb0i
        || fa2i == fb1i || fa2i == fb2i;

    return !share_endpoint
        && (can_vertices_collide(fa0i, fb0i) //
            || can_vertices_collide(fa0i, fb1i)
            || can_vertices_collide(fa0i, fb2i)
            || can_vertices_collide(fa1i, fb0i)
            || can_vertices_collide(fa1i, fb1i)
            || can_vertices_collide(fa1i, fb2i)
            || can_vertices_collide(fa2i, fb0i)
            || can_vertices_collide(fa2i, fb1i)
            || can_vertices_collide(fa2i, fb2i));
}

#endif

} // namespace ipc