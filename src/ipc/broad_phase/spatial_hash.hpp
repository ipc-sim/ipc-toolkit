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
    ArrayMax3d leftBottomCorner, rightTopCorner;
    ArrayMax3i voxelCount;
    double one_div_voxelSize;
    int voxelCount0x1;

    int edgeStartInd, triStartInd;

    unordered_map<int, std::vector<int>> voxel;
    std::vector<std::vector<int>> pointAndEdgeOccupancy;

protected:
    int dim;
    double builtInRadius;

public: // constructor
    SpatialHash() { }

    SpatialHash(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0,
        double voxelSize = -1)
    {
        build(vertices, edges, faces, inflation_radius, voxelSize);
    }

    SpatialHash(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0,
        double voxelSize = -1)
    {
        build(
            vertices_t0, vertices_t1, edges, faces, inflation_radius,
            voxelSize);
    }

public: // API
    void build(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0) override
    {
        build(vertices, edges, faces, inflation_radius, /*voxelSize=*/-1);
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
            /*voxelSize=*/-1);
    }

    void build(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius,
        double voxelSize);

    void build(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius,
        double voxelSize);

    void clear() override
    {
        BroadPhase::clear();
        voxel.clear();
        pointAndEdgeOccupancy.clear();
    }

    void queryPointForTriangles(
        const VectorMax3d& p,
        unordered_set<int>& triInds,
        double radius = 0) const;

    void queryPointForTriangles(
        const VectorMax3d& p_t0,
        const VectorMax3d& p_t1,
        unordered_set<int>& triInds,
        double radius = 0) const;

    void queryPointForPrimitives(
        const VectorMax3d& p_t0,
        const VectorMax3d& p_t1,
        unordered_set<int>& vertInds,
        unordered_set<int>& edgeInds,
        unordered_set<int>& triInds,
        double radius = 0) const;

    void queryEdgeForPE(
        const VectorMax3d& e0,
        const VectorMax3d& e1,
        std::vector<int>& vertInds,
        std::vector<int>& edgeInds,
        double radius = 0) const;

    void queryEdgeForEdges(
        const VectorMax3d& ea0,
        const VectorMax3d& ea1,
        std::vector<int>& edgeInds,
        double radius = 0,
        int eai = -1) const;

    void queryEdgeForEdgesWithBBoxCheck(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const VectorMax3d& ea0,
        const VectorMax3d& ea1,
        std::vector<int>& edgeInds,
        double radius = 0,
        int eai = -1) const;

    void queryEdgeForEdges(
        const VectorMax3d& ea0_t0,
        const VectorMax3d& ea1_t0,
        const VectorMax3d& ea0_t1,
        const VectorMax3d& ea1_t1,
        std::vector<int>& edgeInds,
        double radius = 0,
        int eai = -1) const;

    void queryTriangleForPoints(
        const VectorMax3d& t0,
        const VectorMax3d& t1,
        const VectorMax3d& t2,
        unordered_set<int>& pointInds,
        double radius = 0) const;

    void queryTriangleForPoints(
        const VectorMax3d& t0_t0,
        const VectorMax3d& t1_t0,
        const VectorMax3d& t2_t0,
        const VectorMax3d& t0_t1,
        const VectorMax3d& t1_t1,
        const VectorMax3d& t2_t1,
        unordered_set<int>& pointInds,
        double radius = 0) const;

    void queryTriangleForEdges(
        const VectorMax3d& t0,
        const VectorMax3d& t1,
        const VectorMax3d& t2,
        unordered_set<int>& edgeInds,
        double radius = 0) const;

    void queryEdgeForTriangles(
        const VectorMax3d& e0,
        const VectorMax3d& e1,
        unordered_set<int>& triInds,
        double radius = 0) const;

    void queryPointForPrimitives(
        int vi,
        unordered_set<int>& vertInds,
        unordered_set<int>& edgeInds,
        unordered_set<int>& triInds) const;

    void queryPointForEdges(int vi, unordered_set<int>& edgeInds) const;

    void queryPointForTriangles(int vi, unordered_set<int>& triInds) const;

    // will only put edges with larger than ei index into edgeInds
    void queryEdgeForEdges(int eai, unordered_set<int>& edgeInds) const;

    void queryEdgeForEdgesWithBBoxCheck(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        int eai,
        unordered_set<int>& edgeInds) const;

    void queryEdgeForTriangles(int ei, unordered_set<int>& triInds) const;

    ////////////////////////////////////////////////////////////////////////////
    // BroadPhase API

    /// @brief Find the candidate edge-vertex collisisons.
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
    int locateVoxelIndex(const VectorMax3d& p) const;

    void locateVoxelAxisIndex(
        const VectorMax3d& p, ArrayMax3i& voxelAxisIndex) const;

    void locateBoxVoxelAxisIndex(
        ArrayMax3d minCorner,
        ArrayMax3d maxCorner,
        ArrayMax3i& minIndex,
        ArrayMax3i& maxIndex,
        const double inflation_radius = 0) const;

    int voxelAxisIndex2VoxelIndex(const ArrayMax3i& voxelAxisIndex) const;

    int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const;
};

} // namespace ipc
