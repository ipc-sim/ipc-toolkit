#pragma once

// #include <ccdgpu/record.hpp>
#include <string_view>
#include <typeinfo>

#include "Broadphase/Algorithms/BF/BF_Base.h"
#include "Broadphase/Algorithms/BaseAlgorithm.h"
#include "Broadphase/Algorithms/DBVT/DBVT.h"
#include "Broadphase/Algorithms/iSAP/AxisSweep.h"
#include "Broadphase/ObjectPair.h"
#include <algorithm>
#include <ipc/ccd/ccd.hpp>
#include <set>

#include <ipc/broad_phase/aabb.hpp>

// #include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/candidates/candidates.hpp>

namespace ipc {

class InterfaceBase {
public:
    Candidates candidates;
    int m_Objects;
    long long m_narrowPhase = 0;
    long long m_broadPhase = 0;
    long long m_truePositive = 0;
    virtual void Clear();
    virtual ~InterfaceBase() = default; // placeholder to create vtable
};

template <typename T> class Interface : public InterfaceBase {

public:
    T* algo;

    // broadphase::Overlaps overlaps;

    Interface();
    // virtual ~Interface() = default;
    void FilterOverlaps(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces);
    void CalcOverlaps(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        std::vector<broadmark::Aabb>& broadmark_aabbs,
        bool init = false);
    void ConstructBoxes(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        std::vector<broadmark::Aabb>& broadmark_aabbs);
};

template <>
void Interface<DBVT_F>::FilterOverlaps(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces);

template <>
void Interface<DBVT_D>::FilterOverlaps(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces);

template <>
void Interface<AxisSweep>::FilterOverlaps(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces);

void to_aabbs(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<broadmark::Aabb>& broadmark_aabbs);

// void to_simplices(
//     const Eigen::MatrixXd& V0,
//     const Eigen::MatrixXi& E,
//     std::vector<simplexType>& simplices);

void mesh_to_aabbs(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<ipc::AABB>& edge_aabbs,
    std::vector<ipc::AABB>& face_aabbs,
    std::vector<ipc::AABB>& vertex_aabbs,
    const double inflation_radius = 0);

// void toSimplexAabbs(
//     const Eigen::MatrixXd& V0,
//     const Eigen::MatrixXd& V1,
//     const Eigen::MatrixXi& E,
//     const Eigen::MatrixXi& F,
//     SimplexAabbs & S
// );

void combine_aabbs(
    std::vector<ipc::AABB>& edge_aabbs,
    std::vector<ipc::AABB>& face_aabbs,
    std::vector<ipc::AABB>& vertex_aabbs,
    std::vector<ipc::AABB>& aabbs);

broadmark::Aabb buildWorldAabb(const std::vector<broadmark::Aabb>& aabbs);

void growAabbs(std::vector<broadmark::Aabb>& aabbs, const Vec3& amount);

} // namespace ipc

#include "interface.tpp"