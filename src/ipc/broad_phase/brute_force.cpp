#include "brute_force.hpp"

#include <ipc/utils/merge_thread_local.hpp>

#include <tbb/blocked_range2d.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm> // std::min/max

using namespace std::placeholders;

namespace ipc {

template <typename Candidate, bool triangular>
void BruteForce::detect_candidates(
    const std::vector<AABB>& boxes0,
    const std::vector<AABB>& boxes1,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    tbb::parallel_for(
        tbb::blocked_range2d<size_t>(0ul, boxes0.size(), 0ul, boxes1.size()),
        [&](const tbb::blocked_range2d<size_t>& r) {
            auto& local_candidates = storage.local();

            size_t i_end;
            if constexpr (triangular) {
                i_end = std::min(r.rows().end(), r.cols().end()); // i < j
            } else {
                i_end = r.rows().end();
            }

            for (size_t i = r.rows().begin(); i < i_end; i++) {
                const AABB& box0 = boxes0[i];

                size_t j_begin;
                if constexpr (triangular) {
                    // i < r.cols().end() → i + 1 <= r.cols().end()
                    j_begin = std::max(i + 1, r.cols().begin());
                } else {
                    j_begin = r.cols().begin();
                }

                for (size_t j = j_begin; j < r.cols().end(); j++) {
                    if (!can_collide(i, j)) {
                        continue;
                    }

                    const AABB& box1 = boxes1[j];
                    if (box0.intersects(box1)) {
                        local_candidates.emplace_back(i, j);
                    }
                }
            }
        });

    merge_thread_local_vectors(storage, candidates);
}

void BruteForce::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    detect_candidates<VertexVertexCandidate, true>(
        vertex_boxes, vertex_boxes, can_vertices_collide, candidates);
}

void BruteForce::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    detect_candidates(
        edge_boxes, vertex_boxes,
        std::bind(&BruteForce::can_edge_vertex_collide, this, _1, _2),
        candidates);
}

void BruteForce::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    detect_candidates<EdgeEdgeCandidate, true>(
        edge_boxes, edge_boxes,
        std::bind(&BruteForce::can_edges_collide, this, _1, _2), candidates);
}

void BruteForce::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    detect_candidates(
        face_boxes, vertex_boxes,
        std::bind(&BruteForce::can_face_vertex_collide, this, _1, _2),
        candidates);
}

void BruteForce::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    detect_candidates(
        edge_boxes, face_boxes,
        std::bind(&BruteForce::can_edge_face_collide, this, _1, _2),
        candidates);
}

void BruteForce::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    detect_candidates<FaceFaceCandidate, true>(
        face_boxes, face_boxes,
        std::bind(&BruteForce::can_faces_collide, this, _1, _2), candidates);
}

} // namespace ipc
