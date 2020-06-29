#pragma once

#include <Utils/MESHIO.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace JGSL {

template <class T, int dim>
class SPATIAL_HASH {
public: // data
    Eigen::Matrix<T, dim, 1> leftBottomCorner, rightTopCorner;
    T one_div_voxelSize;
    Eigen::Array<int, dim, 1> voxelCount;
    int voxelCount0x1;

    int surfEdgeStartInd, surfTriStartInd;

    std::unordered_map<int, std::vector<int>> voxel;
    std::vector<std::vector<int>> pointAndEdgeOccupancy; // for CCD

public: // constructor
    SPATIAL_HASH(void) {}

public: // API
    void Build(MESH_NODE<T, dim>& X,
        const std::vector<int>& boundaryNode,
        const std::vector<VECTOR<int, 2>>& boundaryEdge,
        const std::vector<VECTOR<int, 3>>& boundaryTri, 
        T voxelSize)
    {
        BASE_STORAGE<int> threadsP(boundaryNode.size());
        for (int i = 0; i < boundaryNode.size(); ++i) {
            threadsP.Append(i);
        }
        BASE_STORAGE<int> threadsE(boundaryEdge.size());
        for (int i = 0; i < boundaryEdge.size(); ++i) {
            threadsE.Append(i);
        }
        BASE_STORAGE<int> threadsT;
        if constexpr (dim == 3) {
            threadsT.Reserve(boundaryTri.size());
            for (int i = 0; i < boundaryTri.size(); ++i) {
                threadsT.Append(i);
            }
        }

        Eigen::Matrix<T, Eigen::Dynamic, dim> V(X.size, dim);
        X.Par_Each([&](int id, auto data) {
            auto &[x] = data;
            V(id, 0) = x[0];
            V(id, 1) = x[1];
            if constexpr (dim == 3) {
                V(id, 2) = x[2];
            }
        });

        Eigen::Matrix<T, Eigen::Dynamic, 1> eLen(boundaryEdge.size(), 1);
        threadsE.Par_Each([&](int seCount, auto data) {
            const auto& seI = boundaryEdge[seCount];
            eLen[seCount] = (V.row(seI[0]) - V.row(seI[1])).norm();
        });
        voxelSize *= eLen.mean();

        leftBottomCorner = V.colwise().minCoeff();
        rightTopCorner = V.colwise().maxCoeff();
        one_div_voxelSize = 1.0 / voxelSize;
        Eigen::Array<T, dim, 1> range = rightTopCorner - leftBottomCorner;
        voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
        if (voxelCount.minCoeff() <= 0) {
            // cast overflow due to huge search direction or tiny voxelSize
            one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
            voxelCount.setOnes();
        }
        voxelCount0x1 = voxelCount[0] * voxelCount[1];

        surfEdgeStartInd = boundaryNode.size();
        surfTriStartInd = surfEdgeStartInd + boundaryEdge.size();

        // precompute svVAI
        std::vector<Eigen::Array<int, dim, 1>> svVoxelAxisIndex(boundaryNode.size());
        std::vector<int> vI2SVI(X.size);
        threadsP.Par_Each([&](int svI, auto data) {
            int vI = boundaryNode[svI];
            Locate_Voxel_Axis_Index(V.row(vI).transpose(), svVoxelAxisIndex[svI]);
            vI2SVI[vI] = svI;
        });

        voxel.clear();

#ifdef PARALLEL_SH_CONSTRUCT
        std::vector<std::pair<int, int>> voxel_tmp;

        for (int svI = 0; svI < boundaryNode.size(); ++svI) {
            voxel_tmp.emplace_back(Locate_Voxel_Index(V.row(boundaryNode[svI]).transpose()), svI);
        }
#else
        for (int svI = 0; svI < boundaryNode.size(); ++svI) {
            voxel[Locate_Voxel_Index(V.row(boundaryNode[svI]).transpose())].emplace_back(svI);
        }
#endif

        std::vector<std::vector<int>> voxelLoc_e(boundaryEdge.size());
        threadsE.Par_Each([&](int seCount, auto data) {
            const auto& seI = boundaryEdge[seCount];
            const Eigen::Array<int, dim, 1>& voxelAxisIndex_first = svVoxelAxisIndex[vI2SVI[seI[0]]];
            const Eigen::Array<int, dim, 1>& voxelAxisIndex_second = svVoxelAxisIndex[vI2SVI[seI[1]]];
            Eigen::Array<int, dim, 1> mins = voxelAxisIndex_first.min(voxelAxisIndex_second);
            Eigen::Array<int, dim, 1> maxs = voxelAxisIndex_first.max(voxelAxisIndex_second);
            if constexpr (dim == 3) {
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_e[seCount].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
            else {
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yOffset = iy * voxelCount[0];
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        voxelLoc_e[seCount].emplace_back(ix + yOffset);
                    }
                }
            }
        });

        std::vector<std::vector<int>> voxelLoc_sf;
        if constexpr (dim == 3) {
            voxelLoc_sf.resize(boundaryTri.size());
            threadsT.Par_Each([&](int sfI, auto data) {
                const Eigen::Array<int, dim, 1>& voxelAxisIndex0 = svVoxelAxisIndex[vI2SVI[boundaryTri[sfI][0]]];
                const Eigen::Array<int, dim, 1>& voxelAxisIndex1 = svVoxelAxisIndex[vI2SVI[boundaryTri[sfI][1]]];
                const Eigen::Array<int, dim, 1>& voxelAxisIndex2 = svVoxelAxisIndex[vI2SVI[boundaryTri[sfI][2]]];
                Eigen::Array<int, dim, 1> mins = voxelAxisIndex0.min(voxelAxisIndex1).min(voxelAxisIndex2);
                Eigen::Array<int, dim, 1> maxs = voxelAxisIndex0.max(voxelAxisIndex1).max(voxelAxisIndex2);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            });
        }

#ifdef PARALLEL_SH_CONSTRUCT
        for (int seCount = 0; seCount < voxelLoc_e.size(); ++seCount) {
            for (const auto& voxelI : voxelLoc_e[seCount]) {
                voxel_tmp.emplace_back(voxelI, seCount + surfEdgeStartInd);
            }
        }

        if constexpr (dim == 3) {
            for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
                for (const auto& voxelI : voxelLoc_sf[sfI]) {
                    voxel_tmp.emplace_back(voxelI, sfI + surfTriStartInd);
                }
            }
        }

#ifdef USE_TBB
        tbb::parallel_sort(voxel_tmp.begin(), voxel_tmp.end(), [](const std::pair<int, int>& f, const std::pair<int, int>& s) { return f.first < s.first; });
#else
        std::sort(voxel_tmp.begin(), voxel_tmp.end(), [](const std::pair<int, int>& f, const std::pair<int, int>& s) { return f.first < s.first; });
#endif
        std::vector<std::pair<int, std::vector<int>>> voxel_tmp_merged;
        voxel_tmp_merged.reserve(voxel_tmp.size());
        int current_voxel = -1;
        for (const auto& v : voxel_tmp) {
            if (current_voxel != v.first) {
                assert(current_voxel < v.first);
                voxel_tmp_merged.emplace_back();
                voxel_tmp_merged.back().first = v.first;
                current_voxel = v.first;
            }

            voxel_tmp_merged.back().second.push_back(v.second);
        }
        assert(voxel_tmp_merged.size() <= voxel_tmp.size());

        voxel.insert(voxel_tmp_merged.begin(), voxel_tmp_merged.end());
#else
        for (int seCount = 0; seCount < voxelLoc_e.size(); ++seCount) {
            for (const auto& voxelI : voxelLoc_e[seCount]) {
                voxel[voxelI].emplace_back(seCount + surfEdgeStartInd);
            }
        }
        if constexpr (dim == 3) {
            for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
                for (const auto& voxelI : voxelLoc_sf[sfI]) {
                    voxel[voxelI].emplace_back(sfI + surfTriStartInd);
                }
            }
        }
#endif
    }

    void Query_Point_For_Triangles(const Eigen::Matrix<T, 3, 1>& pos,
        T radius, std::unordered_set<int>& triInds) const
    {
        Eigen::Array<int, 3, 1> mins, maxs;
        Locate_Voxel_Axis_Index(pos.array() - radius, mins);
        Locate_Voxel_Axis_Index(pos.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 3, 1>::Zero());
        maxs = maxs.min(voxelCount - 1);

        triInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfTriStartInd) {
                                triInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void Query_Edge_For_Edges(const Eigen::Matrix<T, dim, 1>& vBegin,
        const Eigen::Matrix<T, dim, 1>& vEnd,
        T radius, std::vector<int>& edgeInds, int eIq = -1) const
    {
        Eigen::Matrix<T, dim, 1> leftBottom = vBegin.array().min(vEnd.array()) - radius;
        Eigen::Matrix<T, dim, 1> rightTop = vBegin.array().max(vEnd.array()) + radius;
        Eigen::Array<int, dim, 1> mins, maxs;
        Locate_Voxel_Axis_Index(leftBottom, mins);
        Locate_Voxel_Axis_Index(rightTop, maxs);
        mins = mins.max(Eigen::Array<int, dim, 1>::Zero());
        maxs = maxs.min(voxelCount - 1);

        edgeInds.resize(0);
        if constexpr (dim == 3) {
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        const auto voxelI = voxel.find(ix + yzOffset);
                        if (voxelI != voxel.end()) {
                            for (const auto& indI : voxelI->second) {
                                if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > eIq) {
                                    edgeInds.emplace_back(indI - surfEdgeStartInd);
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yOffset = iy * voxelCount[0];
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI - surfEdgeStartInd > eIq) {
                                edgeInds.emplace_back(indI - surfEdgeStartInd);
                            }
                        }
                    }
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
    }

    void Query_Point_For_Edges(const Eigen::Matrix<T, dim, 1>& pos,
        T radius, std::unordered_set<int>& eInds) const
    {
        Eigen::Array<int, dim, 1> mins, maxs;
        Locate_Voxel_Axis_Index(pos.array() - radius, mins);
        Locate_Voxel_Axis_Index(pos.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, dim, 1>::Zero());
        maxs = maxs.min(voxelCount - 1);

        eInds.clear();
        if constexpr (dim == 3) {
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        const auto voxelI = voxel.find(ix + yzOffset);
                        if (voxelI != voxel.end()) {
                            for (const auto& indI : voxelI->second) {
                                if (indI >= surfEdgeStartInd && indI < surfTriStartInd) {
                                    eInds.insert(indI - surfEdgeStartInd);
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yOffset = iy * voxelCount[0];
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd) {
                                eInds.insert(indI - surfEdgeStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    // for CCD:
    void Build(MESH_NODE<T, dim>& X,
        const std::vector<int>& boundaryNode,
        const std::vector<VECTOR<int, 2>>& boundaryEdge,
        const std::vector<VECTOR<int, 3>>& boundaryTri, 
        const std::vector<T>& searchDir,
        T& curMaxStepSize, T voxelSize)
    {
        BASE_STORAGE<int> threadsP(boundaryNode.size());
        for (int i = 0; i < boundaryNode.size(); ++i) {
            threadsP.Append(i);
        }
        BASE_STORAGE<int> threadsE(boundaryEdge.size());
        for (int i = 0; i < boundaryEdge.size(); ++i) {
            threadsE.Append(i);
        }
        BASE_STORAGE<int> threadsT;
        if constexpr (dim == 3) {
            threadsT.Reserve(boundaryTri.size());
            for (int i = 0; i < boundaryTri.size(); ++i) {
                threadsT.Append(i);
            }
        }

        Eigen::Matrix<T, Eigen::Dynamic, 1> eLen(boundaryEdge.size(), 1);
        threadsE.Par_Each([&](int seCount, auto data) {
            const auto& seI = boundaryEdge[seCount];
            const VECTOR<T, dim>& v0 = std::get<0>(X.Get_Unchecked(seI[0]));
            const VECTOR<T, dim>& v1 = std::get<0>(X.Get_Unchecked(seI[1]));
            eLen[seCount] = (v0 - v1).length();
        });
        voxelSize *= eLen.mean();

        T pSize = 0;
        for (int svI = 0; svI < boundaryNode.size(); ++svI) {
            int vI = boundaryNode[svI];
            pSize += std::abs(searchDir[vI * dim]);
            pSize += std::abs(searchDir[vI * dim + 1]);
            if constexpr (dim == 3) {
                pSize += std::abs(searchDir[vI * dim + 2]);
            }
        }
        pSize /= boundaryNode.size() * dim;

        const T spanSize = pSize / voxelSize;
        std::cout << "span size = " << spanSize << std::endl;
        if (spanSize > 1) {
            curMaxStepSize /= spanSize;
            std::cout << "curMaxStepSize reduced" << std::endl;
        }

        Eigen::Matrix<T, Eigen::Dynamic, dim> SV(boundaryNode.size(), dim);
        Eigen::Matrix<T, Eigen::Dynamic, dim> SVt(boundaryNode.size(), dim);
        std::unordered_map<int, int> vI2SVI;
        for (int svI = 0; svI < boundaryNode.size(); ++svI) {
            int vI = boundaryNode[svI];
            vI2SVI[vI] = svI;

            const VECTOR<T, dim>& v = std::get<0>(X.Get_Unchecked(vI));
            SV(svI, 0) = v[0];
            SVt(svI, 0) = v[0] + curMaxStepSize * searchDir[vI * dim];
            SV(svI, 1) = v[1];
            SVt(svI, 1) = v[1] + curMaxStepSize * searchDir[vI * dim + 1];
            if constexpr (dim == 3) {
                SV(svI, 2) = v[2];
                SVt(svI, 2) = v[2] + curMaxStepSize * searchDir[vI * dim + 2];
            }
        }

        leftBottomCorner = SV.colwise().minCoeff().array().min(SVt.colwise().minCoeff().array());
        rightTopCorner = SV.colwise().maxCoeff().array().max(SVt.colwise().maxCoeff().array());
        one_div_voxelSize = 1.0 / voxelSize;
        Eigen::Array<double, dim, 1> range = rightTopCorner - leftBottomCorner;
        voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
        if (voxelCount.minCoeff() <= 0) {
            // cast overflow due to huge search direction
            one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
            voxelCount.setOnes();
        }
        voxelCount0x1 = voxelCount[0] * voxelCount[1];

        surfEdgeStartInd = boundaryNode.size();
        surfTriStartInd = surfEdgeStartInd + boundaryEdge.size();

        // precompute svVAI
        std::vector<Eigen::Array<int, dim, 1>> svMinVAI(boundaryNode.size());
        std::vector<Eigen::Array<int, dim, 1>> svMaxVAI(boundaryNode.size());
        threadsP.Par_Each([&](int svI, auto data) {
            int vI = boundaryNode[svI];
            Eigen::Matrix<T, dim, 1> minCoord = SV.row(svI).array().min(SVt.row(svI).array());
            Eigen::Matrix<T, dim, 1> maxCoord = SV.row(svI).array().max(SVt.row(svI).array());
            Locate_Voxel_Axis_Index(minCoord, svMinVAI[svI]);
            Locate_Voxel_Axis_Index(maxCoord, svMaxVAI[svI]);
        });

        voxel.clear(); //TODO: parallel insert
        pointAndEdgeOccupancy.resize(0);
        pointAndEdgeOccupancy.resize(surfTriStartInd);

        threadsP.Par_Each([&](int svI, auto data) {
            const Eigen::Array<int, dim, 1>& mins = svMinVAI[svI];
            const Eigen::Array<int, dim, 1>& maxs = svMaxVAI[svI];
            pointAndEdgeOccupancy[svI].reserve((maxs - mins + 1).prod());
            if constexpr (dim == 3) {
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            pointAndEdgeOccupancy[svI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
            else {
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yOffset = iy * voxelCount[0];
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        pointAndEdgeOccupancy[svI].emplace_back(ix + yOffset);
                    }
                }
            }
        });

        threadsE.Par_Each([&](int seCount, auto data) {
            int seIInd = seCount + surfEdgeStartInd;
            const auto& seI = boundaryEdge[seCount];

            Eigen::Array<int, dim, 1> mins = svMinVAI[vI2SVI[seI[0]]].min(svMinVAI[vI2SVI[seI[1]]]);
            Eigen::Array<int, dim, 1> maxs = svMaxVAI[vI2SVI[seI[0]]].max(svMaxVAI[vI2SVI[seI[1]]]);
            pointAndEdgeOccupancy[seIInd].reserve((maxs - mins + 1).prod());
            if constexpr (dim == 3) {
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            pointAndEdgeOccupancy[seIInd].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
            else {
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yOffset = iy * voxelCount[0];
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        pointAndEdgeOccupancy[seIInd].emplace_back(ix + yOffset);
                    }
                }
            }
        });

        std::vector<std::vector<int>> voxelLoc_sf;
        if constexpr (dim == 3) {
            voxelLoc_sf.resize(boundaryTri.size());
            threadsT.Par_Each([&](int sfI, auto data) {
                Eigen::Array<int, dim, 1> mins = svMinVAI[vI2SVI[boundaryTri[sfI][0]]].min(svMinVAI[vI2SVI[boundaryTri[sfI][1]]]).min(svMinVAI[vI2SVI[boundaryTri[sfI][2]]]);
                Eigen::Array<int, dim, 1> maxs = svMaxVAI[vI2SVI[boundaryTri[sfI][0]]].max(svMaxVAI[vI2SVI[boundaryTri[sfI][1]]]).max(svMaxVAI[vI2SVI[boundaryTri[sfI][2]]]);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            });
        }

        for (int i = 0; i < pointAndEdgeOccupancy.size(); ++i) {
            for (const auto& voxelI : pointAndEdgeOccupancy[i]) {
                voxel[voxelI].emplace_back(i);
            }
        }
        if constexpr (dim == 3) {
            for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
                for (const auto& voxelI : voxelLoc_sf[sfI]) {
                    voxel[voxelI].emplace_back(sfI + surfTriStartInd);
                }
            }
        }
    }

    void Query_Point_For_Primitives(int svI, 
        std::unordered_set<int>& sEdgeInds,
        std::unordered_set<int>& sTriInds) const
    {
        sTriInds.clear();
        sEdgeInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfTriStartInd) {
                    sTriInds.insert(indI - surfTriStartInd);
                }
                else if (indI >= surfEdgeStartInd) {
                    sEdgeInds.insert(indI - surfEdgeStartInd);
                }
            }
        }
    }

    // will only put edges with larger than seI index into sEdgeInds
    void Query_Edge_For_Edges(int seI, std::unordered_set<int>& sEdgeInds) const
    {
        sEdgeInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[seI + surfEdgeStartInd]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > seI) {
                    sEdgeInds.insert(indI - surfEdgeStartInd);
                }
            }
        }
    }

    void Query_Point_For_Edges(int svI, std::unordered_set<int>& sEdgeInds) const
    {
        sEdgeInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfEdgeStartInd && indI < surfTriStartInd) {
                    sEdgeInds.insert(indI - surfEdgeStartInd);
                }
            }
        }
    }

public: // helper functions
    int Locate_Voxel_Index(const Eigen::Matrix<T, dim, 1>& pos) const
    {
        Eigen::Array<int, dim, 1> voxelAxisIndex;
        Locate_Voxel_Axis_Index(pos, voxelAxisIndex);
        if constexpr (dim == 3) {
            return voxelAxisIndex[0] + voxelAxisIndex[1] * voxelCount[0] + voxelAxisIndex[2]* voxelCount0x1;
        }
        else {
            return voxelAxisIndex[0] + voxelAxisIndex[1] * voxelCount[0];
        }
    }
    void Locate_Voxel_Axis_Index(const Eigen::Matrix<T, dim, 1>& pos,
        Eigen::Array<int, dim, 1>& voxelAxisIndex) const
    {
        voxelAxisIndex = ((pos - leftBottomCorner) * one_div_voxelSize).array().floor().template cast<int>();
    }
};

}