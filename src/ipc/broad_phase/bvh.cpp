#include "bvh.hpp"

#include <ipc/utils/merge_thread_local.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {
void BVH::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double inflation_radius)
{
    BroadPhase::build(V, E, F, inflation_radius);
    init_bvh(vertex_boxes, vertex_bvh);
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);
}

void BVH::build(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double inflation_radius)
{
    BroadPhase::build(V0, V1, E, F, inflation_radius);
    init_bvh(vertex_boxes, vertex_bvh);
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);
}

void BVH::init_bvh(const std::vector<AABB>& boxes, SimpleBVH::BVH& bvh)
{
    if (boxes.size() == 0)
        return;

    std::vector<std::array<Eigen::Vector3d, 2>> vector_boxes(boxes.size());
    for (int i = 0; i < boxes.size(); i++) {
        vector_boxes[i] = { { to_3D(boxes[i].min), to_3D(boxes[i].max) } };
    }

    bvh.init(vector_boxes);
}

void BVH::clear()
{
    vertex_bvh.clear();
    edge_bvh.clear();
    face_bvh.clear();
}

template <typename Candidate, bool swap_order, bool triangular>
void BVH::detect_candidates(
    const std::vector<AABB>& boxes,
    const SimpleBVH::BVH& bvh,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates)
{
    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0ul, boxes.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_candidates = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {

                std::vector<unsigned int> js;
                bvh.intersect_box(to_3D(boxes[i].min), to_3D(boxes[i].max), js);

                for (const unsigned int j : js) {
                    int ai = i, bi = j;
                    if constexpr (swap_order) {
                        std::swap(ai, bi);
                    }

                    if constexpr (triangular) {
                        if (ai >= bi) {
                            continue;
                        }
                    }

                    if (!can_collide(ai, bi)) {
                        continue;
                    }

                    local_candidates.emplace_back(ai, bi);
                }
            }
        });

    merge_thread_local_vectors(storage, candidates);
}

void BVH::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (vertex_boxes.size() == 0) {
        return;
    }

    detect_candidates<
        VertexVertexCandidate, /*swap_order=*/false, /*triangular=*/true>(
        vertex_boxes, vertex_bvh, can_vertices_collide, candidates);
}

void BVH::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    // O(n^2) or O(n^3) to build
    // O(klog(n)) to do a single look up
    // O(knlog(n)) to do all look ups

    if (edge_boxes.size() == 0 || vertex_boxes.size() == 0) {
        return;
    }

    if (edge_boxes.size() > vertex_boxes.size()) {
        detect_candidates<EdgeVertexCandidate, true>(
            vertex_boxes, edge_bvh,
            [&](size_t ei, size_t vi) {
                return can_edge_vertex_collide(ei, vi);
            },
            candidates);
    } else {
        detect_candidates<EdgeVertexCandidate, false>(
            edge_boxes, vertex_bvh,
            [&](size_t ei, size_t vi) {
                return can_edge_vertex_collide(ei, vi);
            },
            candidates);
    }
}

void BVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    if (edge_boxes.size() == 0) {
        return;
    }

    detect_candidates<
        EdgeEdgeCandidate, /*swap_order=*/false, /*triangular=*/true>(
        edge_boxes, edge_bvh,
        [&](size_t eai, size_t ebi) { return can_edges_collide(eai, ebi); },
        candidates);
}

void BVH::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    if (face_boxes.size() == 0 || vertex_boxes.size() == 0) {
        return;
    }

    if (face_boxes.size() > vertex_boxes.size()) {
        detect_candidates<FaceVertexCandidate, true>(
            vertex_boxes, face_bvh,
            [&](size_t fi, size_t vi) {
                return can_face_vertex_collide(fi, vi);
            },
            candidates);
    } else {
        detect_candidates<FaceVertexCandidate, false>(
            face_boxes, vertex_bvh,
            [&](size_t fi, size_t vi) {
                return can_face_vertex_collide(fi, vi);
            },
            candidates);
    }
}

void BVH::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    if (edge_boxes.size() == 0 || face_boxes.size() == 0) {
        return;
    }

    if (edge_boxes.size() > face_boxes.size()) {
        detect_candidates<EdgeFaceCandidate, true>(
            face_boxes, edge_bvh,
            [&](size_t ei, size_t fi) { return can_edge_face_collide(ei, fi); },
            candidates);
    } else {
        detect_candidates<EdgeFaceCandidate, false>(
            edge_boxes, face_bvh,
            [&](size_t ei, size_t fi) { return can_edge_face_collide(ei, fi); },
            candidates);
    }
}
} // namespace ipc