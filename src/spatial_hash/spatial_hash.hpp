// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#pragma once

#include <vector>

#include <Eigen/Core>
#include <ipc/utils/eigen_ext.hpp>

#include <ipc/spatial_hash/collision_candidate.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

namespace ipc {

class SpatialHash {
public: // data
    VectorMax3d leftBottomCorner, rightTopCorner;
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
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SpatialHash() {}

    SpatialHash(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0,
        double voxelSize = -1)
    {
        build(V, E, F, inflation_radius, voxelSize);
    }

    SpatialHash(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0,
        double voxelSize = -1)
    {
        build(V0, V1, E, F, inflation_radius, voxelSize);
    }

public: // API
    void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0,
        double voxelSize = -1);

    void build(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double inflation_radius = 0,
        double voxelSize = -1);

    void clear()
    {
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
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
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
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        int eai,
        unordered_set<int>& edgeInds) const;

    void queryMeshForCandidates(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        Candidates& candidates,
        bool queryEV = false,
        bool queryEE = true,
        bool queryFV = true) const;

    void queryMeshForCandidates(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        Candidates& candidates,
        bool queryEV = false,
        bool queryEE = true,
        bool queryFV = true) const;

public: // helper functions
    int locateVoxelIndex(const VectorMax3d& p) const;

    void locateVoxelAxisIndex(
        const VectorMax3d& p, ArrayMax3i& voxelAxisIndex) const;

    int voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]) const;

    int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const;
};

} // namespace ipc
