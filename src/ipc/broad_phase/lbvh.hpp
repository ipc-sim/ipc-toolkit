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

        // true to use absolute pointers (left/right child pointer is the
        // absolute index of the child in the buffer/array) or false for
        // relative pointers (left/right child pointer is the relative pointer
        // from the parent index to the child index in the buffer, i.e. absolute
        // child pointer = absolute parent pointer + relative child pointer)
        static constexpr bool ABSOLUTE_POINTERS = true;

        // helper function to handle relative pointers on
        // CPU side, i.e. convert them to absolute
        // pointers for array indexing
        static constexpr uint32_t POINTER(uint32_t index, uint32_t pointer)
        {
            return ABSOLUTE_POINTERS ? pointer : index + pointer;
        }

        /// @brief The min corner of the AABB
        ArrayMax3d aabb_min;
        /// @brief The max corner of the AABB
        ArrayMax3d aabb_max;
        /// @brief The vertex ids (v0, v1, v2)
        std::array<index_t, 3> vertex_ids;
        /// @brief Pointer to the left child or INVALID_POINTER in case of leaf
        int left;
        /// @brief Pointer to the right child or INVALID_POINTER in case of leaf
        int right;
        /// @brief The primitive id (INVALID_ID <=> inner node)
        index_t primitive_id;

        bool is_leaf() const
        {
            assert(is_valid());
            return left == INVALID_POINTER && right == INVALID_POINTER;
        }

        bool is_inner() const
        {
            return left != INVALID_POINTER && right != INVALID_POINTER;
        }

        bool is_valid() const
        {
            return !((left == INVALID_POINTER) ^ (right == INVALID_POINTER));
        }

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
        int parent;
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

    /// @brief BVH containing the vertices.
    std::vector<Node> vertex_bvh;
    /// @brief BVH containing the edges.
    std::vector<Node> edge_bvh;
    /// @brief BVH containing the faces.
    std::vector<Node> face_bvh;
    /// @brief The axis-aligned bounding box of the entire mesh.
    struct {
        ArrayMax3d min;
        ArrayMax3d max;
    } mesh_aabb;
};

} // namespace ipc