#include "bvh.hpp"

#include <ipc/utils/merge_thread_local.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

using namespace std::placeholders;

namespace ipc {

void BVH::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    BroadPhase::build(vertices, edges, faces, inflation_radius);
    init_bvh(vertex_boxes, vertex_bvh);
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);
}

void BVH::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    BroadPhase::build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
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
    // O(n^2) or O(n^3) to build
    // O(klog(n)) to do a single look up
    // O(knlog(n)) to do all look ups

    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), boxes.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_candidates = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                std::vector<unsigned int> js;
                bvh.intersect_box(boxes[i].min, boxes[i].max, js);

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
    if (edge_boxes.size() == 0 || vertex_boxes.size() == 0) {
        return;
    }

    // In 2D and for codimensional edge-vertex collisions, there are more
    // vertices than edges, so we want to iterate over the edges.
    detect_candidates(
        edge_boxes, vertex_bvh,
        std::bind(&BVH::can_edge_vertex_collide, this, _1, _2), candidates);
}

void BVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    if (edge_boxes.size() == 0) {
        return;
    }

    detect_candidates<
        EdgeEdgeCandidate, /*swap_order=*/false, /*triangular=*/true>(
        edge_boxes, edge_bvh, std::bind(&BVH::can_edges_collide, this, _1, _2),
        candidates);
}

void BVH::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    if (face_boxes.size() == 0 || vertex_boxes.size() == 0) {
        return;
    }

    // The ratio vertices:faces is 1:2, so we want to iterate over the vertices.
    detect_candidates<FaceVertexCandidate, /*swap_order=*/true>(
        vertex_boxes, face_bvh,
        std::bind(&BVH::can_face_vertex_collide, this, _1, _2), candidates);
}

void BVH::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    if (edge_boxes.size() == 0 || face_boxes.size() == 0) {
        return;
    }

    // The ratio edges:faces is 3:2, so we want to iterate over the faces.
    detect_candidates<EdgeFaceCandidate, /*swap_order=*/true>(
        face_boxes, edge_bvh,
        std::bind(&BVH::can_edge_face_collide, this, _1, _2), candidates);
}

void BVH::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    if (face_boxes.size() == 0) {
        return;
    }

    detect_candidates<
        FaceFaceCandidate, /*swap_order=*/false, /*triangular=*/true>(
        face_boxes, face_bvh, std::bind(&BVH::can_faces_collide, this, _1, _2),
        candidates);
}
} // namespace ipc