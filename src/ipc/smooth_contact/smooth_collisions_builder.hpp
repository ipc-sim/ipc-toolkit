#pragma once

#include "smooth_collisions.hpp"

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/maybe_parallel_for.hpp>

#include <Eigen/Core>

namespace ipc {

template <int dim> class SmoothCollisionsBuilder;

template <> class SmoothCollisionsBuilder<2> {
public:
    SmoothCollisionsBuilder() { }

    void add_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const SmoothContactParameters& params,
        const std::function<double(const index_t)>& vert_dhat,
        const std::function<double(const index_t)>& edge_dhat,
        const size_t start_i,
        const size_t end_i);

    // -------------------------------------------------------------------------

    static void merge(
        const ParallelCacheType<SmoothCollisionsBuilder<2>>& local_storage,
        SmoothCollisions& merged_collisions);

    // Constructed collisions
    std::vector<std::shared_ptr<SmoothCollision>> collisions;

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<
        std::pair<index_t, index_t>,
        std::shared_ptr<SmoothCollisionTemplate<Point2, Point2>>>
        vert_vert_2_to_id;
    unordered_map<
        std::pair<index_t, index_t>,
        std::shared_ptr<SmoothCollisionTemplate<Edge2, Point2>>>
        vert_edge_2_to_id;
};

template <> class SmoothCollisionsBuilder<3> {
public:
    SmoothCollisionsBuilder() { }

    void add_edge_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const SmoothContactParameters& params,
        const std::function<double(const index_t)>& vert_dhat,
        const std::function<double(const index_t)>& edge_dhat,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const SmoothContactParameters& params,
        const std::function<double(const index_t)>& vert_dhat,
        const std::function<double(const index_t)>& edge_dhat,
        const std::function<double(const index_t)>& face_dhat,
        const size_t start_i,
        const size_t end_i);

    // -------------------------------------------------------------------------

    static void merge(
        const ParallelCacheType<SmoothCollisionsBuilder<3>>& local_storage,
        SmoothCollisions& merged_collisions);

    // Constructed collisions
    std::vector<std::shared_ptr<SmoothCollision>> collisions;

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates, no need for Face-Vertex
    // and Edge-Edge
    unordered_map<
        std::pair<index_t, index_t>,
        std::shared_ptr<SmoothCollisionTemplate<Point3, Point3>>>
        vert_vert_3_to_id;
    unordered_map<
        std::pair<index_t, index_t>,
        std::shared_ptr<SmoothCollisionTemplate<Edge3, Point3>>>
        edge_vert_3_to_id;
};

} // namespace ipc