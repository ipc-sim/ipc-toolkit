// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#pragma once

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <memory>
#include <vector>

namespace ipc {

/// @brief Spatial hash broad phase collision detection.
class SpatialHash : public BroadPhase {
public: // constructor
    SpatialHash();
    ~SpatialHash();

    /// @brief Get the name of the broad phase method.
    /// @return The name of the broad phase method.
    std::string name() const override { return "SpatialHash"; }

    /// @brief Build the spatial hash for static collision detection.
    /// @param vertices Vertex positions
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0) override
    {
        build(vertices, edges, faces, inflation_radius, /*voxel_size=*/-1);
    }

    /// @brief Build the spatial hash for continuous collision detection.
    /// @param vertices_t0 Starting vertices of the vertices.
    /// @param vertices_t1 Ending vertices of the vertices.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0) override
    {
        build(
            vertices_t0, vertices_t1, edges, faces, inflation_radius,
            /*voxel_size=*/-1);
    }

    /// @brief Build the spatial hash for static collision detection.
    /// @param vertices Ending vertices of the vertices.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    /// @param voxel_size Size of the voxels used in the spatial hash.
    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        double inflation_radius,
        double voxel_size)
    {
        build(vertices, vertices, edges, faces, inflation_radius, voxel_size);
    }

    /// @brief Build the spatial hash for continuous collision detection.
    /// @param vertices_t0 Starting vertices of the vertices.
    /// @param vertices_t1 Ending vertices of the vertices.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    /// @param voxel_size Size of the voxels used in the spatial hash.
    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        double inflation_radius,
        double voxel_size);

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

private: // helper functions
    int locate_voxel_index(Eigen::ConstRef<VectorMax3d> p) const;

    ArrayMax3i locate_voxel_axis_index(Eigen::ConstRef<VectorMax3d> p) const;

    void locate_box_voxel_axis_index(
        ArrayMax3d min_corner,
        ArrayMax3d max_corner,
        Eigen::Ref<ArrayMax3i> min_index,
        Eigen::Ref<ArrayMax3i> max_index,
        const double inflation_radius = 0) const;

    int voxel_axis_index_to_voxel_index(
        Eigen::ConstRef<ArrayMax3i> voxel_axis_index) const;

    int voxel_axis_index_to_voxel_index(int ix, int iy, int iz) const;

    // --- Data members -------------------------------------------------------

    /// @brief The left bottom corner of the world bounding box.
    ArrayMax3d left_bottom_corner;

    /// @brief The right top corner of the world bounding box.
    ArrayMax3d right_top_corner;

    /// @brief The number of voxels in each dimension.
    ArrayMax3i voxel_count;

    /// @note Use the Pimpl idiom to hide unordered_map and unordered_set from the public API.
    struct Impl;
    /// @brief Pointer to the implementation details.
    std::unique_ptr<Impl> impl;

    /// @brief 1.0 / voxel_size
    double one_div_voxelSize;

    /// @brief The number of voxels in the first two dimensions.
    int voxel_count_0x1;

    /// @brief The dimension of the space (2D or 3D).
    int dim;
};

} // namespace ipc
