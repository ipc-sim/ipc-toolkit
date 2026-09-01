#pragma once

#include <array>
#include <ipc/collision_mesh.hpp>
#include <ipc/high_order_contact/high_order_collisions.hpp>

#include <Eigen/Core>
#include <tbb/enumerable_thread_specific.h>

#include <memory>

namespace ipc {

template <int dim> class HighOrderCollisionsBuilder;
class PointPotential;
class QuadratureCollisionsBuilder;

template <> class HighOrderCollisionsBuilder<2> {
public:
    HighOrderCollisionsBuilder() = default;
    // Copy creates an empty builder (used by tbb::enumerable_thread_specific).
    HighOrderCollisionsBuilder(const HighOrderCollisionsBuilder&)
        : HighOrderCollisionsBuilder()
    {
    }

    /// @brief Build per-edge, per-QP collision dicts for the 2D quadrature path.
    /// For each edge ei in [start, end), places Gauss-Lobatto QPs on ei and
    /// finds nearby vertices/edges from candidates.ev_set(ei) and
    /// candidates.ee_set(ei). Results are stored in edge_collisions_2d.
    void build_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Candidates& candidates,
        const HighOrderContactParameters& params,
        size_t start,
        size_t end);

    // -------------------------------------------------------------------------

    static void merge(
        tbb::enumerable_thread_specific<HighOrderCollisionsBuilder<2>>&
            local_storage,
        HighOrderCollisions& merged_collisions);

    // Per-edge QP collision dicts: each entry is {edge_id, [dict_qp0, ...]}.
    // Stored as a vector of pairs (not a map) so structured-binding iteration
    // gives mutable references, enabling std::move in merge().
    std::vector<std::pair<
        index_t,
        std::vector<
            std::unique_ptr<HighOrderCollisionDict<PointType::EDGE, 2>>>>>
        edge_collisions_2d;
};

template <> class HighOrderCollisionsBuilder<3> {
public:
    HighOrderCollisionsBuilder() { }

    static std::shared_ptr<HighOrderCollision> reduce_point_triangle_collision(
        const FaceVertexCandidate& candidate,
        const HighOrderContactParameters& params,
        const CollisionMesh& mesh,
        const VertexMatrixView<3>& vertices,
        PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

    static std::shared_ptr<HighOrderCollision> reduce_point_edge_collision(
        const EdgeVertexCandidate& candidate,
        const HighOrderContactParameters& params,
        const CollisionMesh& mesh,
        const VertexMatrixView<3>& vertices,
        PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    void add_face_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const HighOrderContactParameters& params,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_negative_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const HighOrderContactParameters& params,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_positive_vertex_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const HighOrderContactParameters& params,
        const size_t start_i,
        const size_t end_i);

    /*/
    -------------------------------------------------------------------------

    void add_negative_edge_edge_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<std::array<index_t, 3>>& candidates,
        const HighOrderContactParameters& params,
        const double dhat,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_face_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<std::array<index_t, 3>>& candidates,
        const HighOrderContactParameters& params,
        const double dhat,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<std::array<index_t, 3>>& candidates,
        const HighOrderContactParameters& params,
        const double dhat,
        const size_t start_i,
        const size_t end_i);

    /*/// -------------------------------------------------------------------------

    static void merge(
        const tbb::enumerable_thread_specific<HighOrderCollisionsBuilder<3>>&
            local_storage,
        HighOrderCollisions& merged_collisions);

    // Constructed collisions
    std::vector<std::shared_ptr<HighOrderCollision>> collisions;

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<std::pair<index_t, index_t>, index_t> vert_vert_3_to_id;
    unordered_map<std::pair<index_t, index_t>, index_t> vert_edge_3_to_id;
    unordered_map<std::pair<index_t, index_t>, index_t> vert_face_3_to_id;

    unordered_map<std::array<index_t, 3>, index_t> eev_3_to_id;
    unordered_map<std::array<index_t, 3>, index_t> eef_3_to_id;
    unordered_map<std::array<index_t, 3>, index_t> eee_3_to_id;
};

class QuadratureCollisionsBuilder {
public:
    QuadratureCollisionsBuilder(
        const CollisionMesh& mesh,
        const Candidates& candidates,
        const HighOrderContactParameters& params);
    QuadratureCollisionsBuilder(QuadratureCollisionsBuilder&&) = default;
    QuadratureCollisionsBuilder&
    operator=(QuadratureCollisionsBuilder&&) = default;
    QuadratureCollisionsBuilder(const QuadratureCollisionsBuilder& other);
    QuadratureCollisionsBuilder&
    operator=(const QuadratureCollisionsBuilder& other);
    ~QuadratureCollisionsBuilder();

    void build_vertex_collisions(
        const Eigen::MatrixXd& vertices,
        const std::vector<index_t>& vertex_indices,
        size_t start,
        size_t end);

    void build_face_collisions(
        const Eigen::MatrixXd& vertices,
        const std::vector<index_t>& face_indices,
        size_t start,
        size_t end);

    void build_edge_edge_collisions(
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& ee_candidates,
        const size_t start_i,
        const size_t end_i);

    static void merge(
        tbb::enumerable_thread_specific<QuadratureCollisionsBuilder>&
            local_storage,
        HighOrderCollisions& merged_collisions);

    // Local storage
    std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::VERTEX>>>
        vertex_collisions;
    std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::EDGE>>>
        edge_edge_collisions;
    // face_collisions[i] = {fid, [dict_for_qp0, dict_for_qp1, ...]}
    std::vector<std::pair<
        index_t,
        std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::FACE>>>>>
        face_collisions;

    size_t num_collision_pairs = 0;

    std::shared_ptr<PointPotential> point_potential;
};
} // namespace ipc
