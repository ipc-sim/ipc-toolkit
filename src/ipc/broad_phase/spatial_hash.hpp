// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#pragma once

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <vector>

namespace ipc {

class SpatialHash : public BroadPhase {
public: // data
    /// @brief The left bottom corner of the world bounding box.
    ArrayMax3d left_bottom_corner;

    /// @brief The right top corner of the world bounding box.
    ArrayMax3d right_top_corner;

    /// @brief The number of voxels in each dimension.
    ArrayMax3i voxel_count;

    /// @brief 1.0 / voxel_size
    double one_div_voxelSize;

    /// @brief The number of voxels in the first two dimensions.
    int voxel_count_0x1;

    // // The index of the first edge in voxel_occupancies
    int edge_start_ind;
    // // The index of the first triangle in voxel_occupancies
    int tri_start_ind;

    /// @brief Map from voxel index to the primitive indices it contains.
    unordered_map<int, std::vector<int>> voxel_to_primitives;

    /// @brief Map from point index to the voxel indices it occupies.
    std::vector<std::vector<int>> point_to_voxels;

    /// @brief Map from edge index to the voxel indices it occupies.
    std::vector<std::vector<int>> edge_to_voxels;

    /// @brief Map from face index to the voxel indices it occupies.
    std::vector<std::vector<int>> face_to_voxels;

protected:
    int dim;
    double built_in_radius;

public: // constructor
    SpatialHash() = default;

    SpatialHash(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0,
        double voxel_size = -1)
    {
        build(vertices, edges, faces, inflation_radius, voxel_size);
    }

    SpatialHash(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0,
        double voxel_size = -1)
    {
        build(
            vertices_t0, vertices_t1, edges, faces, inflation_radius,
            voxel_size);
    }

public: // API
    void build(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0) override
    {
        build(vertices, edges, faces, inflation_radius, /*voxel_size=*/-1);
    }

    void build(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0) override
    {
        build(
            vertices_t0, vertices_t1, edges, faces, inflation_radius,
            /*voxel_size=*/-1);
    }

    void build(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius,
        double voxel_size);

    void build(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius,
        double voxel_size);

    void clear() override
    {
        BroadPhase::clear();
        voxel_to_primitives.clear();
        point_to_voxels.clear();
        edge_to_voxels.clear();
        face_to_voxels.clear();
    }

    /// @brief Check if primitive index refers to a vertex.
    inline bool is_vertex_index(int idx) const { return idx < edge_start_ind; }

    /// @brief Check if primitive index refers to an edge.
    inline bool is_edge_index(int idx) const
    {
        return idx >= edge_start_ind && idx < tri_start_ind;
    }

    /// @brief Check if primitive index refers to a triangle.
    inline bool is_triangle_index(int idx) const
    {
        return idx >= tri_start_ind;
    }

    /// @brief Convert a primitive index to an edge index.
    inline int to_edge_index(int idx) const
    {
        assert(is_edge_index(idx));
        return idx - edge_start_ind;
    }

    /// @brief Convert a primitive index to a triangle index.
    inline int to_triangle_index(int idx) const
    {
        assert(is_triangle_index(idx));
        return idx - tri_start_ind;
    }

    // ========================================================================
    // BroadPhase API

    /// @brief Find the candidate vertex-vertex collisions.
    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-vertex collisions.
    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-edge collisions.
    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override;

    /// @brief Find the candidate face-vertex collisions.
    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-face intersections.
    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override;

    /// @brief Find the candidate face-face collisions.
    /// @param[out] candidates The candidate face-face collisions.
    void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const override;

protected: // helper functions
    void query_point_for_points(int vi, unordered_set<int>& vert_inds) const;

    void query_point_for_edges(int vi, unordered_set<int>& edge_inds) const;

    void query_point_for_triangles(int vi, unordered_set<int>& tri_inds) const;

    // will only put edges with larger than ei index into edge_inds
    void query_edge_for_edges(int eai, unordered_set<int>& edge_inds) const;

    void query_edge_for_triangles(int ei, unordered_set<int>& tri_inds) const;

    // will only put triangles with larger than ti index into tri_inds
    void
    query_triangle_for_triangles(int ti, unordered_set<int>& tri_inds) const;

    int locate_voxel_index(const VectorMax3d& p) const;

    ArrayMax3i locate_voxel_axis_index(const VectorMax3d& p) const;

    void locate_box_voxel_axis_index(
        ArrayMax3d min_corner,
        ArrayMax3d max_corner,
        ArrayMax3i& min_index,
        ArrayMax3i& max_index,
        const double inflation_radius = 0) const;

    int
    voxel_axis_index_to_voxel_index(const ArrayMax3i& voxel_axis_index) const;

    int voxel_axis_index_to_voxel_index(int ix, int iy, int iz) const;

private:
    /// @brief Detect candidate collisions between type A and type B.
    /// @tparam Candidate Type of candidate collision.
    /// @tparam swap_order Whether to swap the order of A and B when adding to the candidates.
    /// @tparam triangular Whether to consider (i, j) and (j, i) as the same.
    /// @param[in] boxesA The boxes of type A to detect collisions with.
    /// @param[in] boxesB The boxes of type B to detect collisions with.
    /// @param[in] query_A_for_Bs Function to query boxes of type B for boxes of type A.
    /// @param[in] can_collide Function to determine if two primitives can collide given their ids.
    /// @param[out] candidates The candidate collisions.
    template <typename Candidate, bool swap_order, bool triangular = false>
    void detect_candidates(
        const std::vector<AABB>& boxesA,
        const std::vector<AABB>& boxesB,
        const std::function<void(int, unordered_set<int>&)>& query_A_for_Bs,
        const std::function<bool(int, int)>& can_collide,
        std::vector<Candidate>& candidates) const;

    /// @brief Detect candidate collisions between type A and type A.
    /// @tparam Candidate Type of candidate collision.
    /// @param[in] boxesA The boxes of type A to detect collisions with.
    /// @param[in] query_A_for_As Function to query boxes of type A for boxes of type A.
    /// @param[in] can_collide Function to determine if two primitives can collide given their ids.
    /// @param[out] candidates The candidate collisions.
    template <typename Candidate>
    void detect_candidates(
        const std::vector<AABB>& boxesA,
        const std::function<void(int, unordered_set<int>&)>& query_A_for_As,
        const std::function<bool(int, int)>& can_collide,
        std::vector<Candidate>& candidates) const;
};

} // namespace ipc
