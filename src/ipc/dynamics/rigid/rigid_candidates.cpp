#include "rigid_candidates.hpp"

#include <ipc/config.hpp>
#include <ipc/ipc.hpp>
#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/io/write_candidates_obj.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <igl/remove_unreferenced.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace ipc::rigid {

void RigidCandidates::build(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses,
    const double inflation_radius,
    BroadPhase* broad_phase)
{
    build(bodies, poses, poses, inflation_radius, broad_phase);
}

void RigidCandidates::build(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double inflation_radius,
    BroadPhase* broad_phase)
{
    assert(poses_t0.size() == bodies.num_bodies());
    assert(poses_t1.size() == bodies.num_bodies());
    if (bodies.num_bodies() < 2) {
        return; // No collisions possible
    }

    std::unique_ptr<BroadPhase> default_broad_phase;
    if (broad_phase == nullptr) {
        default_broad_phase = make_default_broad_phase();
        broad_phase = default_broad_phase.get();
    }

    const int dim = poses_t0[0].position.size();

    clear();

    /* TODO: Set can_vertices_collide to use rigid body grouping and
     * bodies->can_vertices_collide */

    // 1. Broad phase between bodies

    // a. Build body AABBs
    AABBs body_boxes(bodies.num_bodies());
    for (int i = 0; i < bodies.num_bodies(); i++) {
        const double r = bodies[i].bounding_radius();
        body_boxes[i] = AABB(
            AABB::from_point(poses_t0[i].position, r),
            AABB::from_point(poses_t1[i].position, r));
    }

    // b. Build the broad phase
    broad_phase->build(
        body_boxes, /*edges=*/Eigen::MatrixXi(), /*faces=*/Eigen::MatrixXi(),
        dim);

    // c. Detect body-body candidates: bodies are stored as "vertices," so we
    // can reuse vertex-vertex candidate detection.
    std::vector<VertexVertexCandidate> body_candidates;
    broad_phase->detect_vertex_vertex_candidates(body_candidates);

    // 2. Broad phase between colliding bodies
    for (const VertexVertexCandidate& body_candidate : body_candidates) {
        auto [body_i, body_j] = body_candidate; // Body indices

        // Ensure bi has more vertices than bj to take advantage of log(n)
        // traversal.
        if (bodies.num_body_vertices(body_i)
            < bodies.num_body_vertices(body_j)) {
            std::swap(body_i, body_j);
        }

        // a. Build the boxes for body bj inside the body space of body bi
        //    TODO: This should use the nonlinear trajectory of the body
        const Pose body_j_to_i_t0 =
            poses_t0[body_i].inverse() * poses_t0[body_j];
        const Eigen::MatrixXd j_vertices_t0 =
            bodies.body_vertices(body_i, body_j_to_i_t0);

        const Pose body_j_to_i_t1 =
            poses_t1[body_i].inverse() * poses_t1[body_j];
        const Eigen::MatrixXd j_vertices_t1 =
            bodies.body_vertices(body_j, body_j_to_i_t1);

        AABBs body_j_vertex_boxes(bodies.num_body_vertices(body_j));
        build_vertex_boxes(
            j_vertices_t0, j_vertices_t1, body_j_vertex_boxes,
            inflation_radius);

        AABBs body_j_edge_boxes, body_j_face_boxes;
        build_edge_boxes(
            body_j_vertex_boxes, bodies.body_edges(body_j), body_j_edge_boxes);
        build_face_boxes(
            body_j_vertex_boxes, bodies.body_faces(body_j), body_j_face_boxes);

        // b. Detect candidates between body bi and bj using the BVH of body bi
        if (dim == 2) {
            // TODO: Edge-vertex candidates
        } else {
            // TODO: Edge-edge candidates
            // TODO: Face-vertex candidates
        }
    }
}

} // namespace ipc::rigid
