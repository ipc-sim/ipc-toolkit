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
    ArrayMax3d left_bottom_corner, right_top_corner;
    ArrayMax3i voxel_count;
    double one_div_voxelSize;
    int voxel_count_0x1;

    int edge_start_ind, tri_start_ind;

    unordered_map<int, std::vector<int>> voxel;
    std::vector<std::vector<int>> point_and_edge_occupancy;

protected:
    int dim;
    double built_in_radius;

public: // constructor
    SpatialHash() { }

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
        voxel.clear();
        point_and_edge_occupancy.clear();
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

protected: // helper functions
    void query_point_for_points(int vi, unordered_set<int>& vert_inds) const;

    void query_point_for_edges(int vi, unordered_set<int>& edge_inds) const;

    void query_point_for_triangles(int vi, unordered_set<int>& tri_inds) const;

    // will only put edges with larger than ei index into edge_inds
    void query_edge_for_edges(int eai, unordered_set<int>& edge_inds) const;

    void query_edge_for_triangles(int ei, unordered_set<int>& tri_inds) const;

    int locate_voxel_index(const VectorMax3d& p) const;

    void locate_voxel_axis_index(
        const VectorMax3d& p, ArrayMax3i& voxel_axis_index) const;

    void locate_box_voxel_axis_index(
        ArrayMax3d min_corner,
        ArrayMax3d max_corner,
        ArrayMax3i& min_index,
        ArrayMax3i& max_index,
        const double inflation_radius = 0) const;

    int voxelAxisIndex2VoxelIndex(const ArrayMax3i& voxel_axis_index) const;

    int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const;
};

} // namespace ipc
