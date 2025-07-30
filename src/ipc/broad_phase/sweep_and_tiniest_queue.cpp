#include "sweep_and_tiniest_queue.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <scalable_ccd/cuda/broad_phase/broad_phase.cuh>

namespace ipc {

struct SweepAndTiniestQueue::Boxes {
    ~Boxes() = default;

    std::vector<scalable_ccd::cuda::AABB> vertices;
    std::vector<scalable_ccd::cuda::AABB> edges;
    std::vector<scalable_ccd::cuda::AABB> faces;
};

SweepAndTiniestQueue::SweepAndTiniestQueue()
    : BroadPhase()
    , boxes(std::make_unique<Boxes>())
{
}

SweepAndTiniestQueue::~SweepAndTiniestQueue() = default;

void SweepAndTiniestQueue::build(
    Eigen::ConstRef<Eigen::MatrixXd> _vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    clear();

    // Make sure the vertices are 3D
    const Eigen::MatrixXd vertices = to_X3d(_vertices);

    scalable_ccd::cuda::build_vertex_boxes(
        vertices, boxes->vertices, inflation_radius);
    scalable_ccd::cuda::build_edge_boxes(boxes->vertices, edges, boxes->edges);
    scalable_ccd::cuda::build_face_boxes(boxes->vertices, faces, boxes->faces);
}

void SweepAndTiniestQueue::build(
    Eigen::ConstRef<Eigen::MatrixXd> _vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> _vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    assert(_vertices_t0.rows() == _vertices_t1.rows());
    assert(_vertices_t0.cols() == _vertices_t1.cols());
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    clear();

    // Mutable copies of the vertices
    const Eigen::MatrixXd vertices_t0 = to_X3d(_vertices_t0);
    const Eigen::MatrixXd vertices_t1 = to_X3d(_vertices_t1);

    scalable_ccd::cuda::build_vertex_boxes(
        vertices_t0, vertices_t1, boxes->vertices, inflation_radius);
    scalable_ccd::cuda::build_edge_boxes(boxes->vertices, edges, boxes->edges);
    scalable_ccd::cuda::build_face_boxes(boxes->vertices, faces, boxes->faces);
}

void SweepAndTiniestQueue::clear()
{
    BroadPhase::clear();
    boxes->vertices.clear();
    boxes->edges.clear();
    boxes->faces.clear();
}

void SweepAndTiniestQueue::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    scalable_ccd::cuda::BroadPhase broad_phase;
    // TODO: Precompute d_vertex_boxes
    broad_phase.build(
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->vertices));

    for (const auto& [vai, vbi] : broad_phase.detect_overlaps()) {
        if (can_vertices_collide(vai, vbi)) {
            candidates.emplace_back(vai, vbi);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    scalable_ccd::cuda::BroadPhase broad_phase;
    // TODO: Precompute d_vertex_boxes and d_edge_boxes
    broad_phase.build(
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->edges),
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->vertices));

    for (const auto& [ei, vi] : broad_phase.detect_overlaps()) {
        if (can_edge_vertex_collide(ei, vi)) {
            candidates.emplace_back(ei, vi);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    scalable_ccd::cuda::BroadPhase broad_phase;
    // TODO: Precompute d_edge_boxes
    broad_phase.build(
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->edges));

    for (const auto& [eai, ebi] : broad_phase.detect_overlaps()) {
        if (can_edges_collide(eai, ebi)) {
            candidates.emplace_back(eai, ebi);
        }
    }
}

void SweepAndTiniestQueue::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    scalable_ccd::cuda::BroadPhase broad_phase;
    // TODO: Precompute d_vertex_boxes and d_face_boxes
    broad_phase.build(
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->faces),
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->vertices));

    for (const auto& [fi, vi] : broad_phase.detect_overlaps()) {
        if (can_face_vertex_collide(fi, vi)) {
            candidates.emplace_back(fi, vi);
        }
    }
}

void SweepAndTiniestQueue::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    scalable_ccd::cuda::BroadPhase broad_phase;
    // TODO: Precompute d_face_boxes and d_edge_boxes
    broad_phase.build(
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->edges),
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->faces));

    for (const auto& [ei, fi] : broad_phase.detect_overlaps()) {
        if (can_edge_face_collide(ei, fi)) {
            candidates.emplace_back(ei, fi);
        }
    }
}

void SweepAndTiniestQueue::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    scalable_ccd::cuda::BroadPhase broad_phase;
    // TODO: Precompute d_face_boxes
    broad_phase.build(
        std::make_shared<scalable_ccd::cuda::DeviceAABBs>(boxes->faces));

    for (const auto& [fai, fbi] : broad_phase.detect_overlaps()) {
        if (can_faces_collide(fai, fbi)) {
            candidates.emplace_back(fai, fbi);
        }
    }
}

// ----------------------------------------------------------------------------

bool SweepAndTiniestQueue::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const auto& [e0i, e1i, _] = boxes->edges[ei].vertex_ids;

    // Checked by scalable_ccd
    assert(vi != e0i && vi != e1i);

    return can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i);
}

bool SweepAndTiniestQueue::can_edges_collide(size_t eai, size_t ebi) const
{
    const auto& [ea0i, ea1i, _] = boxes->edges[eai].vertex_ids;
    const auto& [eb0i, eb1i, __] = boxes->edges[ebi].vertex_ids;

    // Checked by scalable_ccd
    assert(ea0i != eb0i && ea0i != eb1i && ea1i != eb0i && ea1i != eb1i);

    return can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
        || can_vertices_collide(ea1i, eb0i) || can_vertices_collide(ea1i, eb1i);
}

bool SweepAndTiniestQueue::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const auto& [f0i, f1i, f2i] = boxes->faces[fi].vertex_ids;

    // Checked by scalable_ccd
    assert(vi != f0i && vi != f1i && vi != f2i);

    return can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
        || can_vertices_collide(vi, f2i);
}

bool SweepAndTiniestQueue::can_edge_face_collide(size_t ei, size_t fi) const
{
    const auto& [e0i, e1i, _] = boxes->edges[ei].vertex_ids;
    const auto& [f0i, f1i, f2i] = boxes->faces[fi].vertex_ids;

    // Checked by scalable_ccd
    assert(
        e0i != f0i && e0i != f1i && e0i != f2i && e1i != f0i && e1i != f1i
        && e1i != f2i);

    return can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
        || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
        || can_vertices_collide(e1i, f1i) || can_vertices_collide(e1i, f2i);
}

bool SweepAndTiniestQueue::can_faces_collide(size_t fai, size_t fbi) const
{
    const auto& [fa0i, fa1i, fa2i] = boxes->faces[fai].vertex_ids;
    const auto& [fb0i, fb1i, fb2i] = boxes->faces[fbi].vertex_ids;

    // Checked by scalable_ccd
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

} // namespace ipc

#endif