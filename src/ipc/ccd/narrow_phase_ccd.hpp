#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class NarrowPhaseCCD {
public:
    NarrowPhaseCCD() = default;

    virtual ~NarrowPhaseCCD() = default;

    virtual bool point_point_ccd(
        const VectorMax3d& p0_t0,
        const VectorMax3d& p1_t0,
        const VectorMax3d& p0_t1,
        const VectorMax3d& p1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;

    virtual bool point_edge_ccd(
        const VectorMax3d& p_t0,
        const VectorMax3d& e0_t0,
        const VectorMax3d& e1_t0,
        const VectorMax3d& p_t1,
        const VectorMax3d& e0_t1,
        const VectorMax3d& e1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;

    virtual bool point_triangle_ccd(
        const Eigen::Vector3d& p_t0,
        const Eigen::Vector3d& t0_t0,
        const Eigen::Vector3d& t1_t0,
        const Eigen::Vector3d& t2_t0,
        const Eigen::Vector3d& p_t1,
        const Eigen::Vector3d& t0_t1,
        const Eigen::Vector3d& t1_t1,
        const Eigen::Vector3d& t2_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;

    virtual bool edge_edge_ccd(
        const Eigen::Vector3d& ea0_t0,
        const Eigen::Vector3d& ea1_t0,
        const Eigen::Vector3d& eb0_t0,
        const Eigen::Vector3d& eb1_t0,
        const Eigen::Vector3d& ea0_t1,
        const Eigen::Vector3d& ea1_t1,
        const Eigen::Vector3d& eb0_t1,
        const Eigen::Vector3d& eb1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;
};

} // namespace ipc