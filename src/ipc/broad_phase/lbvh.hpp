#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

#include <memory>

namespace ipc {

/// @brief Linear Bounding Volume Hierarchy (LBVH) broad phase collision detection.
class LBVH : public BroadPhase {
public:
    static constexpr index_t INVALID_ID = 0xFFFFFFFF;

    struct Node {
        static constexpr int32_t INVALID_POINTER = 0x0; // do not change

        /// @brief The min corner of the AABB
        Eigen::Array3f aabb_min;
        /// @brief The max corner of the AABB
        Eigen::Array3f aabb_max;

        // Union to overlap Leaf data and Internal Node data.
        // This compresses the Node size to 64 bytes (1 cache line),
        // reducing cache misses during traversal.
        union {
            struct {
                /// @brief Pointer to the left child or INVALID_POINTER in case of leaf
                int left;
                /// @brief Pointer to the right child or INVALID_POINTER in case of leaf
                int right;
            };

            struct {
                /// @brief The primitive id (INVALID_ID <=> inner node)
                /// @note We use this to distinguish leaves from internal nodes to safely access the union.
                int primitive_id;

                /// @brief Marker to indicate this is an inner node
                /// If is_inner == 0 then right = 0 which is INVALID_POINTER
                /// If is_inner != 0 then right = actual right pointer
                int is_inner_marker;
            };
        };

        /// @brief Default constructor
        Node();

        /// @brief Check if this node is an inner node.
        bool is_inner() const { return is_inner_marker; }

        /// @brief Check if this node is a leaf node.
        bool is_leaf() const { return !is_inner(); }

        /// @brief Check if this node is valid.
        bool is_valid() const
        {
            return is_inner()
                ? (left != INVALID_POINTER && right != INVALID_POINTER)
                : primitive_id >= 0;
        }

        /// @brief Check if this node's AABB intersects with another node's AABB.
        bool intersects(const Node& other) const
        {
            return (aabb_min <= other.aabb_max).all()
                && (other.aabb_min <= aabb_max).all();
        }
    };

    struct MortonCodeElement {
        uint64_t morton_code; ///< Key for sorting
        size_t box_id;        ///< Pointer into boxes buffer
    };

    struct ConstructionInfo {
        /// @brief Parent to the parent
        int parent = -1;
        /// @brief Number of threads that arrived
        std::atomic<int> visitation_count = 0;
    };

public:
    LBVH();
    ~LBVH();

    /// @brief Get the name of the broad phase method.
    /// @return The name of the broad phase method.
    std::string name() const override { return "LBVH"; }

    using BroadPhase::build;

    /// @brief Clear any built data.
    void clear() override;

    /// @brief Find the candidate vertex-vertex collisions.
    /// @param[out] candidates The candidate vertex-vertex collisions.
    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-vertex collisions.
    /// @param[out] candidates The candidate edge-vertex collisions.
    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-edge collisions.
    /// @param[out] candidates The candidate edge-edge collisions.
    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override;

    /// @brief Find the candidate face-vertex collisions.
    /// @param[out] candidates The candidate face-vertex collisions.
    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-face intersections.
    /// @param[out] candidates The candidate edge-face intersections.
    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override;

    /// @brief Find the candidate face-face collisions.
    /// @param[out] candidates The candidate face-face collisions.
    void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const override;

    const std::vector<Node>& vertex_nodes() const { return vertex_bvh; }
    const std::vector<Node>& edge_nodes() const { return edge_bvh; }
    const std::vector<Node>& face_nodes() const { return face_bvh; }

protected:
    /// @brief Build the broad phase for collision detection.
    /// @note Assumes the vertex_boxes have been built.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    void build(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) override;

    /// @brief Initialize a LBVH from a set of boxes.
    /// @param[in] boxes Set of boxes to initialize the LBVH with.
    /// @param[out] bvh The LBVH to initialize.
    void
    init_bvh(const std::vector<AABB>& boxes, std::vector<Node>& lbvh) const;

    /// @brief Detect candidate collisions between a LBVH and a sets of boxes.
    /// @tparam Candidate Type of candidate collision.
    /// @tparam swap_order Whether to swap the order of box id with the LBVH id when adding to the candidates.
    /// @tparam triangular Whether to consider (i, j) and (j, i) as the same.
    /// @param[in] boxes The boxes to detect collisions with.
    /// @param[in] bvh The LBVH to detect collisions with.
    /// @param[in] can_collide Function to determine if two primitives can collide given their ids.
    /// @param[out] candidates The candidate collisions.
    template <
        typename Candidate,
        bool swap_order = false,
        bool triangular = false>
    static void detect_candidates(
        const std::vector<Node>& source,
        const std::vector<Node>& target,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates);

    template <typename Candidate>
    static void detect_candidates(
        const std::vector<Node>& source_and_target,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        detect_candidates<Candidate, false, true>(
            source_and_target, source_and_target, can_collide, candidates);
    }

    // -------------------------------------------------------------------------
    // BroadPhase::*_boxes are cleared, so we have to override these methods to
    // not depend on them.
    bool has_vertices() const { return !vertex_bvh.empty(); }
    bool has_edges() const { return !edge_bvh.empty(); }
    bool has_faces() const { return !face_bvh.empty(); }
    bool can_edge_vertex_collide(size_t ei, size_t vi) const override;
    bool can_edges_collide(size_t eai, size_t ebi) const override;
    bool can_face_vertex_collide(size_t fi, size_t vi) const override;
    bool can_edge_face_collide(size_t ei, size_t fi) const override;
    bool can_faces_collide(size_t fai, size_t fbi) const override;
    // -------------------------------------------------------------------------

    /// @brief BVH containing the vertices.
    std::vector<Node> vertex_bvh;
    /// @brief BVH containing the edges.
    std::vector<Node> edge_bvh;
    /// @brief BVH containing the faces.
    std::vector<Node> face_bvh;

    /// @brief Edge vertices in the original mesh order.
    std::vector<std::array<index_t, 2>> edge_vertex_ids;
    /// @brief Face vertices in the original mesh order.
    std::vector<std::array<index_t, 3>> face_vertex_ids;

    /// @brief Dimension of the simulation for which the broad phase was built.
    int dim;

    /// @brief The axis-aligned bounding box of the entire mesh.
    struct {
        Eigen::Array3d min;
        Eigen::Array3d max;
    } mesh_aabb;
};

} // namespace ipc