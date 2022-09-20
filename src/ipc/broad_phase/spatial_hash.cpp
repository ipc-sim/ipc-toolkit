// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#include "spatial_hash.hpp"

#include <ipc/ccd/aabb.hpp>
#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/utils/merge_thread_local_vectors.hpp>

#include <ipc/config.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

// Uncomment this to construct spatial hash in parallel.
// #define IPC_TOOLKIT_PARALLEL_SH_CONSTRUCT

namespace ipc {

void SpatialHash::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double inflation_radius,
    double voxelSize)
{
    build(V, V, E, F, inflation_radius, voxelSize);
}

void SpatialHash::build(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double inflation_radius,
    double voxelSize)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    BroadPhase::build(V0, V1, E, F, inflation_radius); // also calls clear()

    dim = V0.cols();
    builtInRadius = inflation_radius;

    if (voxelSize <= 0) {
        voxelSize = suggest_good_voxel_size(V0, V1, E, inflation_radius);
    }

    leftBottomCorner =
        V0.colwise().minCoeff().cwiseMin(V1.colwise().minCoeff()).array()
        - inflation_radius;
    rightTopCorner =
        V0.colwise().maxCoeff().cwiseMax(V1.colwise().maxCoeff()).array()
        + inflation_radius;
    one_div_voxelSize = 1.0 / voxelSize;
    ArrayMax3d range = rightTopCorner - leftBottomCorner;
    voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxelCount.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxelCount.setOnes();
    }
    voxelCount0x1 = voxelCount[0] * voxelCount[1];

    edgeStartInd = V0.rows();
    triStartInd = edgeStartInd + E.rows();

    // precompute vVAI
    std::vector<Eigen::Array3i> vertexMinVAI(V0.rows(), Eigen::Array3i::Zero());
    std::vector<Eigen::Array3i> vertexMaxVAI(V0.rows(), Eigen::Array3i::Zero());
    tbb::parallel_for(size_t(0), size_t(V0.rows()), [&](size_t vi) {
        ArrayMax3i vVAIMin, vVAIMax;
        locateVoxelAxisIndex(
            V0.row(vi).cwiseMin(V1.row(vi)).array() - inflation_radius,
            vVAIMin);
        locateVoxelAxisIndex(
            V0.row(vi).cwiseMax(V1.row(vi)).array() + inflation_radius,
            vVAIMax);
        vertexMinVAI[vi].head(dim) = vVAIMin;
        vertexMaxVAI[vi].head(dim) = vVAIMax;
    });

    // #ifdef IPC_TOOLKIT_PARALLEL_SH_CONSTRUCT
    //     std::vector<std::pair<int, int>> voxel_tmp;
    //
    //     for (int vi = 0; vi < V.rows(); vi++) {
    //         voxel_tmp.emplace_back(locateVoxelIndex(V.row(vi)), vi);
    //     }
    // #else
    //     for (int vi = 0; vi < V.rows(); vi++) {
    //         voxel[locateVoxelIndex(V.row(vi))].emplace_back(vi);
    //     }
    // #endif

    pointAndEdgeOccupancy.resize(triStartInd);

    tbb::parallel_for(size_t(0), size_t(V0.rows()), [&](size_t vi) {
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

    tbb::parallel_for(size_t(0), size_t(E.rows()), [&](size_t ei) {
        int eiInd = ei + edgeStartInd;

        Eigen::Array3i mins =
            vertexMinVAI[E(ei, 0)].min(vertexMinVAI[E(ei, 1)]);
        Eigen::Array3i maxs =
            vertexMaxVAI[E(ei, 0)].max(vertexMaxVAI[E(ei, 1)]);

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

    std::vector<std::vector<int>> voxelLoc_f(F.rows());
    tbb::parallel_for(size_t(0), size_t(F.rows()), [&](size_t fi) {
        Eigen::Array3i mins = vertexMinVAI[F(fi, 0)]
                                  .min(vertexMinVAI[F(fi, 1)])
                                  .min(vertexMinVAI[F(fi, 2)]);
        Eigen::Array3i maxs = vertexMaxVAI[F(fi, 0)]
                                  .max(vertexMaxVAI[F(fi, 1)])
                                  .max(vertexMaxVAI[F(fi, 2)]);

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

    // #ifdef IPC_TOOLKIT_PARALLEL_SH_CONSTRUCT
    //     for (int ei = 0; ei < voxelLoc_e.size(); ei++) {
    //         for (const auto& voxelI : voxelLoc_e[ei]) {
    //             voxel_tmp.emplace_back(voxelI, ei + edgeStartInd);
    //         }
    //     }
    //
    //     for (int fi = 0; fi < voxelLoc_f.size(); fi++) {
    //         for (const auto& voxelI : voxelLoc_f[fi]) {
    //             voxel_tmp.emplace_back(voxelI, fi + triStartInd);
    //         }
    //     }
    //
    //     // Sort the voxels based on the voxel indices
    //     tbb::parallel_sort(
    //         voxel_tmp.begin(), voxel_tmp.end(),
    //         [](const std::pair<int, int>& f, const std::pair<int, int>& s) {
    //             return f.first < s.first;
    //         });
    //
    //     std::vector<std::pair<int, std::vector<int>>> voxel_tmp_merged;
    //     voxel_tmp_merged.reserve(voxel_tmp.size());
    //     int current_voxel = -1;
    //     for (const auto& v : voxel_tmp) {
    //         if (current_voxel != v.first) {
    //             assert(current_voxel < v.first);
    //             voxel_tmp_merged.emplace_back();
    //             voxel_tmp_merged.back().first = v.first;
    //             current_voxel = v.first;
    //         }
    //
    //         voxel_tmp_merged.back().second.push_back(v.second);
    //     }
    //     assert(voxel_tmp_merged.size() <= voxel_tmp.size());
    //
    //     voxel.insert(voxel_tmp_merged.begin(), voxel_tmp_merged.end());
    // #else
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
    // #endif
}

void SpatialHash::queryPointForTriangles(
    const VectorMax3d& p, unordered_set<int>& triInds, double radius) const
{
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(p.array() - radius, mins);
    locateVoxelAxisIndex(p.array() + radius, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    locateVoxelAxisIndex(p_t0.array().min((p_t1).array()) - radius, mins);
    locateVoxelAxisIndex(p_t0.array().max((p_t1).array()) + radius, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    locateVoxelAxisIndex(p_t0.array().min((p_t1).array()) - radius, mins);
    locateVoxelAxisIndex(p_t0.array().max((p_t1).array()) + radius, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    VectorMax3d leftBottom = e0.array().min(e1.array()) - radius;
    VectorMax3d rightTop = e0.array().max(e1.array()) + radius;
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    VectorMax3d leftBottom = e0.array().min(e1.array()) - radius;
    VectorMax3d rightTop = e0.array().max(e1.array()) + radius;
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

void SpatialHash::queryEdgeForEdgesWithBBoxCheck(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const VectorMax3d& ea0,
    const VectorMax3d& ea1,
    std::vector<int>& edgeInds,
    double radius,
    int eai) const
{
    VectorMax3d leftBottom = ea0.array().min(ea1.array()) - radius;
    VectorMax3d rightTop = ea0.array().max(ea1.array()) + radius;
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
                            const VectorMax3d& eb0 = V.row(E(ebi, 0));
                            const VectorMax3d& eb1 = V.row(E(ebi, 1));
                            ArrayMax3d bboxEBTopRight =
                                eb0.array().max(eb1.array());
                            ArrayMax3d bboxEBBottomLeft =
                                eb0.array().min(eb1.array());
                            if (!((bboxEBBottomLeft - rightTop.array() > 0.0)
                                      .any()
                                  || (leftBottom.array() - bboxEBTopRight > 0.0)
                                         .any())) {
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
    VectorMax3d leftBottom = ea0_t0.array()
                                 .min(ea1_t0.array())
                                 .min(ea0_t1.array())
                                 .min(ea1_t1.array());
    VectorMax3d rightTop = ea0_t0.array()
                               .max(ea1_t0.array())
                               .max(ea0_t1.array())
                               .max(ea1_t1.array());
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
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
    VectorMax3d leftBottom =
        t0.array().min(t1.array()).min(t2.array()) - radius;
    VectorMax3d rightTop = t0.array().max(t1.array()).max(t2.array()) + radius;
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    VectorMax3d leftBottom = t0_t0.array()
                                 .min(t1_t0.array())
                                 .min(t2_t0.array())
                                 .min(t0_t1.array())
                                 .min(t1_t1.array())
                                 .min(t2_t1.array())
        - radius;
    VectorMax3d rightTop = t0_t0.array()
                               .max(t1_t0.array())
                               .max(t2_t0.array())
                               .max(t0_t1.array())
                               .max(t1_t1.array())
                               .max(t2_t1.array())
        + radius;
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom.array(), mins);
    locateVoxelAxisIndex(rightTop.array(), maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    VectorMax3d leftBottom = t0.array().min(t1.array()).min(t2.array());
    VectorMax3d rightTop = t0.array().max(t1.array()).max(t2.array());
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
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
    VectorMax3d leftBottom = e0.array().min(e1.array());
    VectorMax3d rightTop = e0.array().max(e1.array());
    ArrayMax3i mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
    mins = mins.max(ArrayMax3i::Zero(dim));
    maxs = maxs.min(voxelCount - 1);

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
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    int eai,
    unordered_set<int>& edgeInds) const
{
    const VectorMax3d& ea0_t0 = V0.row(E(eai, 0));
    const VectorMax3d& ea1_t0 = V0.row(E(eai, 1));
    const VectorMax3d& ea0_t1 = V1.row(E(eai, 0));
    const VectorMax3d& ea1_t1 = V1.row(E(eai, 1));

    ArrayMax3d bboxEATopRight = ea0_t0.array()
                                    .max(ea1_t0.array())
                                    .max(ea0_t1.array())
                                    .max(ea1_t1.array());
    ArrayMax3d bboxEABottomLeft = ea0_t0.array()
                                      .min(ea1_t0.array())
                                      .min(ea0_t1.array())
                                      .min(ea1_t1.array());

    edgeInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[eai + edgeStartInd]) {
        const auto& voxelI = voxel.at(voxelInd);
        for (const auto& indI : voxelI) {
            if (indI >= edgeStartInd && indI < triStartInd
                && indI - edgeStartInd > eai) {
                int ebi = indI - edgeStartInd;
                const VectorMax3d& eb0_t0 = V0.row(E(ebi, 0));
                const VectorMax3d& eb1_t0 = V0.row(E(ebi, 1));
                const VectorMax3d& eb0_t1 = V1.row(E(ebi, 0));
                const VectorMax3d& eb1_t1 = V1.row(E(ebi, 1));

                ArrayMax3d bboxEBTopRight = eb0_t0.array()
                                                .max(eb1_t0.array())
                                                .max(eb0_t1.array())
                                                .max(eb1_t1.array());
                ArrayMax3d bboxEBBottomLeft = eb0_t0.array()
                                                  .min(eb1_t0.array())
                                                  .min(eb0_t1.array())
                                                  .min(eb1_t1.array());

                if (!((bboxEBBottomLeft - bboxEATopRight > 0.0).any()
                      || (bboxEABottomLeft - bboxEBTopRight > 0.0).any())) {
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
        tbb::blocked_range<size_t>(0ul, vertex_boxes.size()),
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
        tbb::blocked_range<size_t>(0ul, edge_boxes.size()),
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
        tbb::blocked_range<size_t>(0ul, vertex_boxes.size()),
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
        tbb::blocked_range<size_t>(0ul, edge_boxes.size()),
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
    voxelAxisIndex = ((p - leftBottomCorner) * one_div_voxelSize)
                         .array()
                         .floor()
                         .template cast<int>();
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
