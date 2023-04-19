// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#include "spatial_hash.hpp"

#include <ipc/ccd/aabb.hpp>
#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/utils/merge_thread_local.hpp>

#include <ipc/config.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

namespace ipc {

void SpatialHash::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double inflation_radius,
    double voxelSize)
{
    build(vertices, vertices, edges, faces, inflation_radius, voxelSize);
}

void SpatialHash::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double inflation_radius,
    double voxelSize)
{
    const size_t num_vertices = vertices_t0.rows();
    dim = vertices_t0.cols();

    assert(vertices_t1.rows() == num_vertices && vertices_t1.cols() == dim);

    // also calls clear()
    BroadPhase::build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    builtInRadius = inflation_radius;

    if (voxelSize <= 0) {
        voxelSize = suggest_good_voxel_size(
            vertices_t0, vertices_t1, edges, inflation_radius);
    }

    leftBottomCorner = vertices_t0.colwise().minCoeff().cwiseMin(
        vertices_t1.colwise().minCoeff());
    rightTopCorner = vertices_t0.colwise().maxCoeff().cwiseMax(
        vertices_t1.colwise().maxCoeff());

    AABB::conservative_inflation(
        leftBottomCorner, rightTopCorner, inflation_radius);

    one_div_voxelSize = 1.0 / voxelSize;
    ArrayMax3d range = rightTopCorner - leftBottomCorner;
    voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxelCount.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxelCount.setOnes();
    }
    voxelCount0x1 = voxelCount[0] * voxelCount[1];

    edgeStartInd = num_vertices;
    triStartInd = edgeStartInd + edges.rows();

    // precompute vVAI
    std::vector<Eigen::Array3i> vertexMinVAI(
        num_vertices, Eigen::Array3i::Zero());
    std::vector<Eigen::Array3i> vertexMaxVAI(
        num_vertices, Eigen::Array3i::Zero());
    tbb::parallel_for(size_t(0), num_vertices, [&](size_t vi) {
        ArrayMax3d v_min = vertices_t0.row(vi).cwiseMin(vertices_t1.row(vi));
        ArrayMax3d v_max = vertices_t0.row(vi).cwiseMax(vertices_t1.row(vi));
        AABB::conservative_inflation(v_min, v_max, inflation_radius);

        ArrayMax3i vVAIMin, vVAIMax;
        locateVoxelAxisIndex(v_min, vVAIMin);
        locateVoxelAxisIndex(v_max, vVAIMax);

        vertexMinVAI[vi].head(dim) = vVAIMin;
        vertexMaxVAI[vi].head(dim) = vVAIMax;
    });

    pointAndEdgeOccupancy.resize(triStartInd);

    tbb::parallel_for(size_t(0), num_vertices, [&](size_t vi) {
        const Eigen::Array3i &mins = vertexMinVAI[vi], &maxs = vertexMaxVAI[vi];
        assert((mins <= maxs).all());
        pointAndEdgeOccupancy[vi].reserve((maxs - mins + 1).prod());
        for (int iz = mins[2]; iz <= maxs[2]; iz++) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; iy++) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                    pointAndEdgeOccupancy[vi].emplace_back(ix + yzOffset);
                }
            }
        }
    });

    tbb::parallel_for(size_t(0), size_t(edges.rows()), [&](size_t ei) {
        int eiInd = ei + edgeStartInd;

        Eigen::Array3i mins =
            vertexMinVAI[edges(ei, 0)].min(vertexMinVAI[edges(ei, 1)]);
        Eigen::Array3i maxs =
            vertexMaxVAI[edges(ei, 0)].max(vertexMaxVAI[edges(ei, 1)]);

        pointAndEdgeOccupancy[eiInd].reserve((maxs - mins + 1).prod());
        for (int iz = mins[2]; iz <= maxs[2]; iz++) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; iy++) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                    pointAndEdgeOccupancy[eiInd].emplace_back(ix + yzOffset);
                }
            }
        }
    });

    std::vector<std::vector<int>> voxelLoc_f(faces.rows());
    tbb::parallel_for(size_t(0), size_t(faces.rows()), [&](size_t fi) {
        Eigen::Array3i mins = vertexMinVAI[faces(fi, 0)]
                                  .min(vertexMinVAI[faces(fi, 1)])
                                  .min(vertexMinVAI[faces(fi, 2)]);
        Eigen::Array3i maxs = vertexMaxVAI[faces(fi, 0)]
                                  .max(vertexMaxVAI[faces(fi, 1)])
                                  .max(vertexMaxVAI[faces(fi, 2)]);

        for (int iz = mins[2]; iz <= maxs[2]; iz++) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; iy++) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                    voxelLoc_f[fi].emplace_back(ix + yzOffset);
                }
            }
        }
    });

    for (int i = 0; i < pointAndEdgeOccupancy.size(); i++) {
        for (const auto& voxelI : pointAndEdgeOccupancy[i]) {
            voxel[voxelI].emplace_back(i);
        }
    }
    for (int fi = 0; fi < voxelLoc_f.size(); fi++) {
        for (const auto& voxelI : voxelLoc_f[fi]) {
            voxel[voxelI].emplace_back(fi + triStartInd);
        }
    }
}

void SpatialHash::queryPointForTriangles(
    const VectorMax3d& p, unordered_set<int>& triInds, double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(p, p, mins, maxs, radius);

    triInds.clear();

    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= triStartInd) {
                            triInds.insert(indI - triStartInd);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryPointForTriangles(
    const VectorMax3d& p_t0,
    const VectorMax3d& p_t1,
    unordered_set<int>& triInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        p_t0.cwiseMin(p_t1), p_t0.cwiseMax(p_t1), mins, maxs, radius);

    triInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= triStartInd) {
                            triInds.insert(indI - triStartInd);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryPointForPrimitives(
    const VectorMax3d& p_t0,
    const VectorMax3d& p_t1,
    unordered_set<int>& vertInds,
    unordered_set<int>& edgeInds,
    unordered_set<int>& triInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        p_t0.cwiseMin(p_t1), p_t0.cwiseMax(p_t1), mins, maxs, radius);

    vertInds.clear();
    edgeInds.clear();
    triInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < edgeStartInd) {
                            vertInds.insert(indI);
                        } else if (indI < triStartInd) {
                            edgeInds.insert(indI - edgeStartInd);
                        } else {
                            triInds.insert(indI - triStartInd);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryEdgeForPE(
    const VectorMax3d& e0,
    const VectorMax3d& e1,
    std::vector<int>& vertInds,
    std::vector<int>& edgeInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        e0.cwiseMin(e1), e0.cwiseMax(e1), mins, maxs, radius);

    vertInds.clear();
    edgeInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < edgeStartInd) {
                            vertInds.emplace_back(indI);
                        } else if (indI < triStartInd) {
                            edgeInds.emplace_back(indI - edgeStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(
        std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
    std::sort(vertInds.begin(), vertInds.end());
    vertInds.erase(
        std::unique(vertInds.begin(), vertInds.end()), vertInds.end());
}

void SpatialHash::queryEdgeForEdges(
    const VectorMax3d& e0,
    const VectorMax3d& e1,
    std::vector<int>& edgeInds,
    double radius,
    int eai) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        e0.cwiseMin(e1), e0.cwiseMax(e1), mins, maxs, radius);

    edgeInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= edgeStartInd && indI < triStartInd
                            && indI - edgeStartInd > eai) {
                            edgeInds.emplace_back(indI - edgeStartInd);
                        }
                    }
                }
            }
        }
    }

    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(
        std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
}

void SpatialHash::queryEdgeForEdgesWithBBoxCheck(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const VectorMax3d& ea0,
    const VectorMax3d& ea1,
    std::vector<int>& edgeInds,
    double radius,
    int eai) const
{
    ArrayMax3d leftBottom = ea0.cwiseMin(ea1);
    ArrayMax3d rightTop = ea0.cwiseMax(ea1);
    AABB::conservative_inflation(leftBottom, rightTop, radius);

    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

    edgeInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= edgeStartInd && indI < triStartInd
                            && indI - edgeStartInd > eai) {
                            int ebi = indI - edgeStartInd;
                            const VectorMax3d& eb0 =
                                vertices.row(edges(ebi, 0));
                            const VectorMax3d& eb1 =
                                vertices.row(edges(ebi, 1));
                            ArrayMax3d bboxEBBottomLeft = eb0.cwiseMin(eb1);
                            ArrayMax3d bboxEBTopRight = eb0.cwiseMax(eb1);
                            if (!((bboxEBBottomLeft > rightTop).any()
                                  || (leftBottom > bboxEBTopRight).any())) {
                                edgeInds.emplace_back(indI - edgeStartInd);
                            }
                        }
                    }
                }
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(
        std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
}

void SpatialHash::queryEdgeForEdges(
    const VectorMax3d& ea0_t0,
    const VectorMax3d& ea1_t0,
    const VectorMax3d& ea0_t1,
    const VectorMax3d& ea1_t1,
    std::vector<int>& edgeInds,
    double radius,
    int eai) const
{
    VectorMax3d leftBottom =
        ea0_t0.cwiseMin(ea1_t0).cwiseMin(ea0_t1).cwiseMin(ea1_t1);
    VectorMax3d rightTop =
        ea0_t0.cwiseMax(ea1_t0).cwiseMax(ea0_t1).cwiseMax(ea1_t1);
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

    edgeInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= edgeStartInd && indI < triStartInd
                            && indI - edgeStartInd > eai) {
                            edgeInds.emplace_back(indI - edgeStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(
        std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
}

void SpatialHash::queryTriangleForPoints(
    const VectorMax3d& t0,
    const VectorMax3d& t1,
    const VectorMax3d& t2,
    unordered_set<int>& pointInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        t0.cwiseMin(t1).cwiseMin(t2), t0.cwiseMax(t1).cwiseMax(t2), mins, maxs,
        radius);

    pointInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < edgeStartInd) {
                            pointInds.insert(indI);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryTriangleForPoints(
    const VectorMax3d& t0_t0,
    const VectorMax3d& t1_t0,
    const VectorMax3d& t2_t0,
    const VectorMax3d& t0_t1,
    const VectorMax3d& t1_t1,
    const VectorMax3d& t2_t1,
    unordered_set<int>& pointInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    // clang-format off
    locateBoxVoxelAxisIndex(
        t0_t0.cwiseMin(t1_t0).cwiseMin(t2_t0).cwiseMin(t0_t1).cwiseMin(t1_t1).cwiseMin(t2_t1),
        t0_t0.cwiseMax(t1_t0).cwiseMax(t2_t0).cwiseMax(t0_t1).cwiseMax(t1_t1).cwiseMax(t2_t1),
        mins, maxs, radius);
    // clang-format on

    pointInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < edgeStartInd) {
                            pointInds.insert(indI);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryTriangleForEdges(
    const VectorMax3d& t0,
    const VectorMax3d& t1,
    const VectorMax3d& t2,
    unordered_set<int>& edgeInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        t0.cwiseMin(t1).cwiseMin(t2), t0.cwiseMax(t1).cwiseMax(t2), mins, maxs,
        radius);

    edgeInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= edgeStartInd && indI < triStartInd) {
                            edgeInds.insert(indI - edgeStartInd);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryEdgeForTriangles(
    const VectorMax3d& e0,
    const VectorMax3d& e1,
    unordered_set<int>& triInds,
    double radius) const
{
    ArrayMax3i mins, maxs;
    locateBoxVoxelAxisIndex(
        e0.cwiseMin(e1), e0.cwiseMax(e1), mins, maxs, radius);

    triInds.clear();
    int min_z = mins.size() >= 3 ? mins[2] : 0;
    int max_z = maxs.size() >= 3 ? maxs[2] : 0;
    for (int iz = min_z; iz <= max_z; iz++) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; iy++) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= triStartInd) {
                            triInds.insert(indI - triStartInd);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::queryPointForPrimitives(
    int vi,
    unordered_set<int>& vertInds,
    unordered_set<int>& edgeInds,
    unordered_set<int>& triInds) const
{
    vertInds.clear();
    edgeInds.clear();
    triInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[vi]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI < edgeStartInd) {
                vertInds.insert(indI);
            } else if (indI < triStartInd) {
                edgeInds.insert(indI - edgeStartInd);
            } else {
                triInds.insert(indI - triStartInd);
            }
        }
    }
}

void SpatialHash::queryPointForEdges(int vi, unordered_set<int>& edgeInds) const
{
    edgeInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[vi]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI >= edgeStartInd && indI < triStartInd) {
                edgeInds.insert(indI - edgeStartInd);
            }
        }
    }
}

void SpatialHash::queryPointForTriangles(
    int vi, unordered_set<int>& triInds) const
{
    triInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[vi]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI >= triStartInd) {
                triInds.insert(indI - triStartInd);
            }
        }
    }
}

// will only put edges with larger than eai index into edgeInds
void SpatialHash::queryEdgeForEdges(int eai, unordered_set<int>& edgeInds) const
{
    edgeInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[eai + edgeStartInd]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI >= edgeStartInd && indI < triStartInd
                && indI - edgeStartInd > eai) {
                edgeInds.insert(indI - edgeStartInd);
            }
        }
    }
}

void SpatialHash::queryEdgeForEdgesWithBBoxCheck(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    int eai,
    unordered_set<int>& edgeInds) const
{
    const VectorMax3d& ea0_t0 = vertices_t0.row(edges(eai, 0));
    const VectorMax3d& ea1_t0 = vertices_t0.row(edges(eai, 1));
    const VectorMax3d& ea0_t1 = vertices_t1.row(edges(eai, 0));
    const VectorMax3d& ea1_t1 = vertices_t1.row(edges(eai, 1));

    const ArrayMax3d bboxEABottomLeft =
        ea0_t0.cwiseMin(ea1_t0).cwiseMin(ea0_t1).cwiseMin(ea1_t1);
    const ArrayMax3d bboxEATopRight =
        ea0_t0.cwiseMax(ea1_t0).cwiseMax(ea0_t1).cwiseMax(ea1_t1);

    edgeInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[eai + edgeStartInd]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI >= edgeStartInd && indI < triStartInd
                && indI - edgeStartInd > eai) {
                int ebi = indI - edgeStartInd;
                const VectorMax3d& eb0_t0 = vertices_t0.row(edges(ebi, 0));
                const VectorMax3d& eb1_t0 = vertices_t0.row(edges(ebi, 1));
                const VectorMax3d& eb0_t1 = vertices_t1.row(edges(ebi, 0));
                const VectorMax3d& eb1_t1 = vertices_t1.row(edges(ebi, 1));

                const ArrayMax3d bboxEBBottomLeft =
                    eb0_t0.cwiseMin(eb1_t0).cwiseMin(eb0_t1).cwiseMin(eb1_t1);
                const ArrayMax3d bboxEBTopRight =
                    eb0_t0.cwiseMax(eb1_t0).cwiseMax(eb0_t1).cwiseMax(eb1_t1);

                if (!((bboxEBBottomLeft > bboxEATopRight).any()
                      || (bboxEABottomLeft > bboxEBTopRight).any())) {
                    edgeInds.insert(indI - edgeStartInd);
                }
            }
        }
    }
}

void SpatialHash::queryEdgeForTriangles(
    int ei, unordered_set<int>& triInds) const
{
    triInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[ei + edgeStartInd]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI >= triStartInd) {
                triInds.insert(indI - triStartInd);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////
// BroadPhase API

void SpatialHash::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<EdgeVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), vertex_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long vi = range.begin(); vi != range.end(); vi++) {
                const AABB& vertex_box = vertex_boxes[vi];

                unordered_set<int> edgeInds;
                queryPointForEdges(vi, edgeInds);

                for (const auto& ei : edgeInds) {
                    if (!can_edge_vertex_collide(ei, vi)) {
                        continue;
                    }

                    const AABB& edge_box = edge_boxes[ei];
                    if (vertex_box.intersects(edge_box)) {
                        local_candidates.emplace_back(ei, vi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<EdgeEdgeCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), edge_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long eai = range.begin(); eai != range.end(); eai++) {
                const AABB& edge_a_box = edge_boxes[eai];

                unordered_set<int> edgeInds;
                queryEdgeForEdges(eai, edgeInds);

                for (const auto& ebi : edgeInds) {
                    if (eai >= ebi || !can_edges_collide(eai, ebi)) {
                        continue;
                    }

                    const AABB& edge_b_box = edge_boxes[ebi];
                    if (edge_a_box.intersects(edge_b_box)) {
                        local_candidates.emplace_back(eai, ebi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<FaceVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), vertex_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long vi = range.begin(); vi != range.end(); vi++) {
                const AABB& vertex_box = vertex_boxes[vi];

                unordered_set<int> triInds;
                queryPointForTriangles(vi, triInds);

                for (const auto& fi : triInds) {
                    if (!can_face_vertex_collide(fi, vi)) {
                        continue;
                    }

                    const AABB& face_box = face_boxes[fi];
                    if (vertex_box.intersects(face_box)) {
                        local_candidates.emplace_back(fi, vi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<EdgeFaceCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), edge_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long ei = range.begin(); ei != range.end(); ei++) {
                const AABB& edge_box = edge_boxes[ei];

                unordered_set<int> triInds;
                queryEdgeForTriangles(ei, triInds);

                for (const auto& fi : triInds) {
                    if (!can_edge_face_collide(ei, fi)) {
                        continue;
                    }

                    const AABB& face_box = face_boxes[fi];
                    if (edge_box.intersects(face_box)) {
                        local_candidates.emplace_back(ei, fi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

////////////////////////////////////////////////////////////////////////////////

int SpatialHash::locateVoxelIndex(const VectorMax3d& p) const
{
    ArrayMax3i voxelAxisIndex;
    locateVoxelAxisIndex(p, voxelAxisIndex);
    return voxelAxisIndex2VoxelIndex(voxelAxisIndex);
}

void SpatialHash::SpatialHash::locateVoxelAxisIndex(
    const VectorMax3d& p, ArrayMax3i& voxelAxisIndex) const
{
    voxelAxisIndex = ((p.array() - leftBottomCorner) * one_div_voxelSize)
                         .floor()
                         .template cast<int>();
}

void SpatialHash::locateBoxVoxelAxisIndex(
    ArrayMax3d minCorner,
    ArrayMax3d maxCorner,
    ArrayMax3i& minIndex,
    ArrayMax3i& maxIndex,
    const double inflation_radius) const
{
    AABB::conservative_inflation(minCorner, maxCorner, inflation_radius);
    locateVoxelAxisIndex(minCorner, minIndex);
    locateVoxelAxisIndex(maxCorner, maxIndex);
    minIndex = minIndex.max(ArrayMax3i::Zero(dim));
    maxIndex = maxIndex.min(voxelCount - 1);
}

int SpatialHash::voxelAxisIndex2VoxelIndex(
    const ArrayMax3i& voxelAxisIndex) const
{
    return voxelAxisIndex2VoxelIndex(
        voxelAxisIndex[0], voxelAxisIndex[1],
        voxelAxisIndex.size() >= 3 ? voxelAxisIndex[2] : 0);
}

int SpatialHash::voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const
{
    assert(ix >= 0 && ix < voxelCount[0]);
    assert(iy >= 0 && iy < voxelCount[1]);
    assert(iz >= 0 && iz < (voxelCount.size() >= 3 ? voxelCount[2] : 1));
    return ix + iy * voxelCount[0] + iz * voxelCount0x1;
}

} // namespace ipc
