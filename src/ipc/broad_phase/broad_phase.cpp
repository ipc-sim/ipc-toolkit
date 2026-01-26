#include "broad_phase.hpp"

#include <ipc/config.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

namespace ipc {

void BroadPhase::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    IPC_TOOLKIT_PROFILE_BLOCK("BroadPhase::build(static)");
    clear();
    dim = static_cast<uint8_t>(vertices.cols());
    build_vertex_boxes(vertices, vertex_boxes, inflation_radius);
    build(edges, faces);
}

void BroadPhase::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    IPC_TOOLKIT_PROFILE_BLOCK("BroadPhase::build(dynamic)");
    assert(vertices_t0.rows() == vertices_t1.rows());
    assert(vertices_t0.cols() == vertices_t1.cols());
    clear();
    dim = static_cast<uint8_t>(vertices_t0.cols());
    build_vertex_boxes(
        vertices_t0, vertices_t1, vertex_boxes, inflation_radius);
    build(edges, faces);
}

void BroadPhase::build(
    const AABBs& _vertex_boxes,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const uint8_t _dim)
{
    IPC_TOOLKIT_PROFILE_BLOCK("BroadPhase::build(boxes)");

    clear();

    assert(&(this->vertex_boxes) != &_vertex_boxes);
    this->vertex_boxes = _vertex_boxes;
    this->dim = _dim;

    build(edges, faces);
}

void BroadPhase::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    IPC_TOOLKIT_PROFILE_BLOCK("BroadPhase::build(edges_faces)");
    assert(!vertex_boxes.empty());
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);
    build_edge_boxes(vertex_boxes, edges, edge_boxes);
    build_face_boxes(vertex_boxes, faces, face_boxes);
}

void BroadPhase::clear()
{
    vertex_boxes.clear();
    edge_boxes.clear();
    face_boxes.clear();
    dim = 0; // reset dimension
}

void BroadPhase::detect_collision_candidates(Candidates& candidates) const
{
    candidates.clear();
    assert(dim == 2 || dim == 3);
    if (dim == 2) {
        // This is not needed for 3D
        detect_edge_vertex_candidates(candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        detect_edge_edge_candidates(candidates.ee_candidates);
        detect_face_vertex_candidates(candidates.fv_candidates);
    }
}

void BroadPhase::compute_mesh_aabb(
    Eigen::Ref<Eigen::Array3d> mesh_min,
    Eigen::Ref<Eigen::Array3d> mesh_max) const
{
    IPC_TOOLKIT_PROFILE_BLOCK("BroadPhase::compute_mesh_aabb");
    assert(!vertex_boxes.empty());

    struct MeshDomain {
        Eigen::Array3d min;
        Eigen::Array3d max;
    };

    MeshDomain domain {
        Eigen::Array3d::Constant(std::numeric_limits<double>::max()),
        Eigen::Array3d::Constant(std::numeric_limits<double>::lowest())
    };

    domain = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, vertex_boxes.size()), domain,
        [&](const tbb::blocked_range<size_t>& r, MeshDomain local) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                local.min = local.min.min(vertex_boxes[i].min);
                local.max = local.max.max(vertex_boxes[i].max);
            }
            return local;
        },
        [](const MeshDomain& a, const MeshDomain& b) {
            return MeshDomain { a.min.min(b.min), a.max.max(b.max) };
        });

    mesh_min = domain.min;
    mesh_max = domain.max;
}

// ============================================================================

bool BroadPhase::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    assert(ei < edge_boxes.size());
    const auto& [e0i, e1i, _] = edge_boxes[ei].vertex_ids;

    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool BroadPhase::can_edges_collide(size_t eai, size_t ebi) const
{
    assert(eai < edge_boxes.size());
    const auto& [ea0i, ea1i, _] = edge_boxes[eai].vertex_ids;
    assert(ebi < edge_boxes.size());
    const auto& [eb0i, eb1i, __] = edge_boxes[ebi].vertex_ids;

    const bool share_endpoint =
        ea0i == eb0i || ea0i == eb1i || ea1i == eb0i || ea1i == eb1i;

    return !share_endpoint
        && (can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
            || can_vertices_collide(ea1i, eb0i)
            || can_vertices_collide(ea1i, eb1i));
}

bool BroadPhase::can_face_vertex_collide(size_t fi, size_t vi) const
{
    assert(fi < face_boxes.size());
    const auto& [f0i, f1i, f2i] = face_boxes[fi].vertex_ids;

    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool BroadPhase::can_edge_face_collide(size_t ei, size_t fi) const
{
    assert(ei < edge_boxes.size());
    const auto& [e0i, e1i, _] = edge_boxes[ei].vertex_ids;
    assert(fi < face_boxes.size());
    const auto& [f0i, f1i, f2i] = face_boxes[fi].vertex_ids;

    const bool share_endpoint = e0i == f0i || e0i == f1i || e0i == f2i
        || e1i == f0i || e1i == f1i || e1i == f2i;

    return !share_endpoint
        && (can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
            || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
            || can_vertices_collide(e1i, f1i)
            || can_vertices_collide(e1i, f2i));
}

bool BroadPhase::can_faces_collide(size_t fai, size_t fbi) const
{
    assert(fai < face_boxes.size());
    const auto& [fa0i, fa1i, fa2i] = face_boxes[fai].vertex_ids;
    assert(fbi < face_boxes.size());
    const auto& [fb0i, fb1i, fb2i] = face_boxes[fbi].vertex_ids;

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
