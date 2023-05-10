#pragma once

#include <ipc/config.hpp>
#ifdef IPC_TOOLKIT_WITH_BROADMARK

#include <ipc/candidates/candidates.hpp>

namespace broadmark {
class Aabb;
} // namespace broadmark
class DBVT_F;
class DBVT_D;
class AxisSweep;
class Vec3;

namespace ipc {

class InterfaceBase {
public:
    virtual ~InterfaceBase() = default;
    virtual void Clear();

    Candidates candidates;
    int m_Objects;
    long long m_narrowPhase = 0;
    long long m_broadPhase = 0;
    long long m_truePositive = 0;
};

template <typename T> class Interface : public InterfaceBase {
public:
    T* algo;

    Interface();

    void FilterOverlaps(
        const long num_vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        Candidates& candidates) const;

    void CalcOverlaps(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        std::vector<broadmark::Aabb>& broadmark_aabbs,
        bool init = false);
};

template <>
void Interface<DBVT_F>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    Candidates& candidates) const;

template <>
void Interface<DBVT_D>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    Candidates& candidates) const;

template <>
void Interface<AxisSweep>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    Candidates& candidates) const;

void to_aabbs(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<broadmark::Aabb>& broadmark_aabbs,
    const double inflation_radius = 0);

void mesh_to_aabbs(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<ipc::AABB>& edge_aabbs,
    std::vector<ipc::AABB>& face_aabbs,
    std::vector<ipc::AABB>& vertex_aabbs,
    const double inflation_radius = 0);

void combine_aabbs(
    const std::vector<ipc::AABB>& edge_aabbs,
    const std::vector<ipc::AABB>& face_aabbs,
    const std::vector<ipc::AABB>& vertex_aabbs,
    std::vector<ipc::AABB>& aabbs);

void to_broadmark_aabbs(
    const std::vector<ipc::AABB>& ipc_aabbs,
    std::vector<broadmark::Aabb>& broadmark_aabbs);

broadmark::Aabb buildWorldAabb(const std::vector<broadmark::Aabb>& aabbs);

void growAabbs(std::vector<broadmark::Aabb>& aabbs, const Vec3& amount);

} // namespace ipc

#include "broadmark_util.tpp"
#endif