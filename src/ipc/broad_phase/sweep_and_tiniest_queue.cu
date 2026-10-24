#include "sweep_and_tiniest_queue.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/broad_phase/details/connectivity_filters.hpp>

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

    dim = _vertices.cols();

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

    dim = _vertices_t0.cols();

    // Mutable copies of the vertices
    const Eigen::MatrixXd vertices_t0 = to_X3d(_vertices_t0);
    const Eigen::MatrixXd vertices_t1 = to_X3d(_vertices_t1);

    scalable_ccd::cuda::build_vertex_boxes(
        vertices_t0, vertices_t1, boxes->vertices, inflation_radius);
    scalable_ccd::cuda::build_edge_boxes(boxes->vertices, edges, boxes->edges);
    scalable_ccd::cuda::build_face_boxes(boxes->vertices, faces, boxes->faces);
}

void SweepAndTiniestQueue::build(
    const AABBs& _vertex_boxes, // not the (unused) inherited member
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const uint8_t _dim)
{
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    clear();

    dim = _dim;

    // Convert from ipc::AABB to scalable_ccd::cuda::AABB
    boxes->vertices.resize(_vertex_boxes.size());
    for (int i = 0; i < _vertex_boxes.size(); ++i) {
        boxes->vertices[i].min.x = _vertex_boxes[i].min.x();
        boxes->vertices[i].min.y = _vertex_boxes[i].min.y();
        boxes->vertices[i].min.z = _vertex_boxes[i].min.z();

        boxes->vertices[i].max.x = _vertex_boxes[i].max.x();
        boxes->vertices[i].max.y = _vertex_boxes[i].max.y();
        boxes->vertices[i].max.z = _vertex_boxes[i].max.z();

        // If vertex id == -1 it means this slot is not used.
        // But Scalable CCD does not have this kind of special value so we map
        // it to unique negative id.
        const auto [vi, vj, vk] = _vertex_boxes[i].vertex_ids;
        assert(vi >= 0);
        boxes->vertices[i].vertex_ids.x = vi;
        boxes->vertices[i].vertex_ids.y = vj >= 0 ? vj : (-vi - 1);
        boxes->vertices[i].vertex_ids.z = vk >= 0 ? vk : (-vi - 1);

        boxes->vertices[i].element_id = i;
    }

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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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

// Scalable CCD already excludes pairs that share a vertex, so these apply only
// the user-filter half of ipc::details' connectivity rule (and assert the
// other half held).

bool SweepAndTiniestQueue::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const auto& [e0i, e1i, _] = boxes->edges[ei].vertex_ids;
    const index_t e[2] = { e0i, e1i };
    const index_t v[1] = { static_cast<index_t>(vi) };

    assert(!details::share_vertex(e, v)); // Checked by scalable_ccd

    return details::any_vertex_pair_can_collide(v, e, can_vertices_collide);
}

bool SweepAndTiniestQueue::can_edges_collide(size_t eai, size_t ebi) const
{
    const auto& [ea0i, ea1i, _] = boxes->edges[eai].vertex_ids;
    const auto& [eb0i, eb1i, __] = boxes->edges[ebi].vertex_ids;
    const index_t ea[2] = { ea0i, ea1i };
    const index_t eb[2] = { eb0i, eb1i };

    assert(!details::share_vertex(ea, eb)); // Checked by scalable_ccd

    return details::any_vertex_pair_can_collide(ea, eb, can_vertices_collide);
}

bool SweepAndTiniestQueue::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const auto& [f0i, f1i, f2i] = boxes->faces[fi].vertex_ids;
    const index_t f[3] = { f0i, f1i, f2i };
    const index_t v[1] = { static_cast<index_t>(vi) };

    assert(!details::share_vertex(f, v)); // Checked by scalable_ccd

    return details::any_vertex_pair_can_collide(v, f, can_vertices_collide);
}

bool SweepAndTiniestQueue::can_edge_face_collide(size_t ei, size_t fi) const
{
    const auto& [e0i, e1i, _] = boxes->edges[ei].vertex_ids;
    const auto& [f0i, f1i, f2i] = boxes->faces[fi].vertex_ids;
    const index_t e[2] = { e0i, e1i };
    const index_t f[3] = { f0i, f1i, f2i };

    assert(!details::share_vertex(e, f)); // Checked by scalable_ccd

    return details::any_vertex_pair_can_collide(e, f, can_vertices_collide);
}

bool SweepAndTiniestQueue::can_faces_collide(size_t fai, size_t fbi) const
{
    const auto& [fa0i, fa1i, fa2i] = boxes->faces[fai].vertex_ids;
    const auto& [fb0i, fb1i, fb2i] = boxes->faces[fbi].vertex_ids;
    const index_t fa[3] = { fa0i, fa1i, fa2i };
    const index_t fb[3] = { fb0i, fb1i, fb2i };

    assert(!details::share_vertex(fa, fb)); // Checked by scalable_ccd

    return details::any_vertex_pair_can_collide(fa, fb, can_vertices_collide);
}

} // namespace ipc

#endif
