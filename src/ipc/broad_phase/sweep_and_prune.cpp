#include "sweep_and_prune.hpp"

#include <scalable_ccd/broad_phase/sort_and_sweep.hpp>

namespace ipc {

void SweepAndPrune::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    clear();

    scalable_ccd::build_vertex_boxes(vertices, vertex_boxes, inflation_radius);
    scalable_ccd::build_edge_boxes(vertex_boxes, edges, edge_boxes);
    scalable_ccd::build_face_boxes(vertex_boxes, faces, face_boxes);
}

void SweepAndPrune::build(
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

void SweepAndPrune::clear()
{
    BroadPhase::clear();
    vertex_boxes.clear();
    edge_boxes.clear();
    face_boxes.clear();
}

void SweepAndPrune::detect_vertex_vertex_candidates(
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

void SweepAndPrune::detect_edge_vertex_candidates(
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

void SweepAndPrune::detect_edge_edge_candidates(
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

void SweepAndPrune::detect_face_vertex_candidates(
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

void SweepAndPrune::detect_edge_face_candidates(
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

void SweepAndPrune::detect_face_face_candidates(
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

bool SweepAndPrune::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const auto& [e0i, e1i, _] = edge_boxes[ei].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(vi != e0i && vi != e1i);

    return can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i);
}

bool SweepAndPrune::can_edges_collide(size_t eai, size_t ebi) const
{
    const auto& [ea0i, ea1i, _] = edge_boxes[eai].vertex_ids;
    const auto& [eb0i, eb1i, __] = edge_boxes[ebi].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(ea0i != eb0i && ea0i != eb1i && ea1i != eb0i && ea1i != eb1i);

    return can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
        || can_vertices_collide(ea1i, eb0i) || can_vertices_collide(ea1i, eb1i);
}

bool SweepAndPrune::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const auto& [f0i, f1i, f2i] = face_boxes[fi].vertex_ids;

    // Checked by scalable_ccd::sort_and_sweep
    assert(vi != f0i && vi != f1i && vi != f2i);

    return can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
        || can_vertices_collide(vi, f2i);
}

bool SweepAndPrune::can_edge_face_collide(size_t ei, size_t fi) const
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

bool SweepAndPrune::can_faces_collide(size_t fai, size_t fbi) const
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

} // namespace ipc