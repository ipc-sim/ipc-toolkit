#pragma once

#include <ipc/collision_mesh.hpp>
#include "smooth_collisions.hpp"
#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

template <int dim> class SmoothCollisionsBuilder {
public:
    SmoothCollisionsBuilder() {}

    // only for 2D, transform edge-vertex to edge-edge
    void add_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const ParameterType& param,
        const std::function<double(const long&)>& vert_dhat,
        const std::function<double(const long&)>& edge_dhat,
        const size_t start_i,
        const size_t end_i);

    // ------------------------------------------------------------------------

    static void merge(
#if defined(IPC_TOOLKIT_WITH_TBB)
        const tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>>&
#elif defined(IPC_TOOLKIT_WITH_CPP_THREADS)
        const std::vector<SmoothCollisionsBuilder<dim>>&
#else
        const std::array<SmoothCollisionsBuilder<dim>, 1>&
#endif
            local_storage,
        SmoothCollisions<dim>& merged_collisions);

    // -------------------------------------------------------------------------

    // only for 3D, transform edge-edge to edge-edge-face
    void add_edge_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const ParameterType& param,
        const std::function<double(const long&)>& vert_dhat,
        const std::function<double(const long&)>& edge_dhat,
        const size_t start_i,
        const size_t end_i);

    // only for 3D, transform face-vertex to face-face
    void add_face_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const ParameterType& param,
        const std::function<double(const long&)>& vert_dhat,
        const std::function<double(const long&)>& edge_dhat,
        const std::function<double(const long&)>& face_dhat,
        const size_t start_i,
        const size_t end_i);

    template <typename TCollision>
    static void add_collision(
        const std::shared_ptr<TCollision>& pair,
        unordered_map<std::pair<long, long>, std::tuple<TCollision, long>>&
            cc_to_id_,
        std::vector<
            std::shared_ptr<typename SmoothCollisions<dim>::value_type>>&
            collisions_);

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<
        std::pair<long, long>,
        std::tuple<SmoothCollisionTemplate<max_vert_2d, Point2, Point2>, long>>
        vert_vert_2_to_id;
    unordered_map<
        std::pair<long, long>,
        std::tuple<SmoothCollisionTemplate<max_vert_2d, Edge2, Point2>, long>>
        vert_edge_2_to_id;

    unordered_map<
        std::pair<long, long>,
        std::tuple<SmoothCollisionTemplate<max_vert_3d, Face, Point3>, long>>
        face_vert_to_id;
    unordered_map<
        std::pair<long, long>,
        std::tuple<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>, long>>
        vert_vert_3_to_id;
    unordered_map<
        std::pair<long, long>,
        std::tuple<SmoothCollisionTemplate<max_vert_3d, Edge3, Point3>, long>>
        edge_vert_3_to_id;
    unordered_map<
        std::pair<long, long>,
        std::tuple<SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3>, long>>
        edge_edge_3_to_id;

    // Constructed collisions
    std::vector<std::shared_ptr<typename SmoothCollisions<dim>::value_type>>
        collisions;
};

} // namespace ipc