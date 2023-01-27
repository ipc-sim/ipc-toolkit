#include "bvh.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {
void BVH::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const const double inflation_radius = 0);
{
    BroadPhase::build(V, E, F, inflation_radius);
    if (V.cols() == 2) {
        if (V.rows() <= E.rows()) { }
    } else {
    }
    // EV
    // EE
    // FV
    // FE
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);
}

void BVH::build(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const const double inflation_radius = 0);
{
    BroadPhase::build(V0, V1, E, F, inflation_radius);
    // EV
    // EE
    // FV
    // FE
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);
}

void BVH::init_bvh(
    const std::vector<AABB>& boxes, const int dim, ::BVH::BVH& bvh)
{
    std::vector<std::array<Eigen::Vector3d, 2>> vector_boxes(boxes.size());
    for (int i = 0; i < boxes.size(); i++) {
        vector_boxes[i][0].head(dim) = boxes.min;
        if (dim < 3) {
            vector_boxes[i][0][2] = 0;
        }

        vector_boxes[i][1].head(dim) = boxes.max;
        if (dim < 3) {
            vector_boxes[i][0][2] = 0;
        }
    }
    bvh.init(vector_boxes);
}

void BVH::clear()
{
    vertex_bvh = ::BVH::BVH();
    edge_bvh = ::BVH::BVH();
    face_bvh = ::BVH::BVH();
}

template <> void detect_candidates() { }

void BVH::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    // O(n^2) or O(n^3) to build
    // O(klog(n)) to do a single look up
    // O(knlog(n)) to do all look ups

    for ()
        intersect_box(
            const Eigen::Vector3d& bbd0, const Eigen::Vector3d& bbd1,
            std::vector<unsigned int>& list)
}

void BVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
}

void BVH::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
}

void BVH::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
}
} // namespace ipc