#include "sweep_and_tiniest_queue.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <scalable_ccd/cuda/broad_phase/broad_phase.cuh>

namespace ipc {

void SweepAndTiniestQueue::build(
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

    stq.build(vertices, edges, faces, boxes, inflation_radius);
    overlaps = stq.detect_overlaps();
}

void SweepAndTiniestQueue::build(
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

    stq.build(vertices_t0, vertices_t1, edges, faces, boxes, inflation_radius);
    overlaps = stq.detect_overlaps();
}

void SweepAndTiniestQueue::clear()
{
    BroadPhase::clear();
    stq.clear();
    overlaps.clear();
    boxes.clear();
}

void SweepAndTiniestQueue::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    throw std::runtime_error(
        "SweepAndTiniestQueue::detect_vertex_vertex_candidates not implemented!");
}

void SweepAndTiniestQueue::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    throw std::runtime_error(
        "SweepAndTiniestQueue::detect_edge_vertex_candidates not implemented!");
}

void SweepAndTiniestQueue::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    for (const std::pair<int, int>& overlap : overlaps) {
        const scalable_ccd::cuda::AABB& boxA = boxes[overlap.first];
        const scalable_ccd::cuda::AABB& boxB = boxes[overlap.second];
        if (boxA.is_edge() && boxB.is_edge()
            && can_edges_collide(boxA.element_id, boxB.element_id)) {
            candidates.emplace_back(boxA.element_id, boxB.element_id);
        }
    }
}

void SweepAndTiniestQueue::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    for (const std::pair<int, int>& overlap : overlaps) {
        const scalable_ccd::cuda::AABB& boxA = boxes[overlap.first];
        const scalable_ccd::cuda::AABB& boxB = boxes[overlap.second];
        if (boxA.is_face() && boxB.is_vertex()
            && can_face_vertex_collide(boxA.element_id, boxB.element_id)) {
            candidates.emplace_back(boxA.element_id, boxB.element_id);
        } else if (
            boxB.is_face() && boxA.is_vertex()
            && can_face_vertex_collide(boxB.element_id, boxA.element_id)) {
            candidates.emplace_back(boxB.element_id, boxA.element_id);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    throw std::runtime_error(
        "SweepAndTiniestQueue::detect_edge_face_candidates not implemented!");
}

void SweepAndTiniestQueue::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    throw std::runtime_error(
        "SweepAndTiniestQueue::detect_face_face_candidates not implemented!");
}

// ----------------------------------------------------------------------------

bool SweepAndTiniestQueue::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const long e0i = edges(ei, 0), e1i = edges(ei, 1);

    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool SweepAndTiniestQueue::can_edges_collide(size_t eai, size_t ebi) const
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

bool SweepAndTiniestQueue::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const long f0i = faces(fi, 0), f1i = faces(fi, 1), f2i = faces(fi, 2);

    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool SweepAndTiniestQueue::can_edge_face_collide(size_t ei, size_t fi) const
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

bool SweepAndTiniestQueue::can_faces_collide(size_t fai, size_t fbi) const
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

} // namespace ipc

#endif