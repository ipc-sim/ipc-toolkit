#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Narrow Phase Continuous Collision Detection (CCD) interface.
class NarrowPhaseCCD {
public:
    NarrowPhaseCCD() = default;

    virtual ~NarrowPhaseCCD() = default;

    /// @brief Perform narrow phase CCD between two points.
    /// @param p0_t0 The starting position of the first point.
    /// @param p1_t0 The starting position of the second point.
    /// @param p0_t1 The ending position of the first point.
    /// @param p1_t1 The ending position of the second point.
    /// @param toi The time of impact.
    /// @param min_distance The minimum distance between the two points.
    /// @param tmax The maximum time to check for collision.
    /// @return True if a collision was detected, false otherwise.
    virtual bool point_point_ccd(
        Eigen::ConstRef<VectorMax3d> p0_t0,
        Eigen::ConstRef<VectorMax3d> p1_t0,
        Eigen::ConstRef<VectorMax3d> p0_t1,
        Eigen::ConstRef<VectorMax3d> p1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;

    /// @brief Perform narrow phase CCD between a point and a linear edge.
    /// @param p_t0 The starting position of the point.
    /// @param e0_t0 The starting position of the first endpoint of the edge.
    /// @param e1_t0 The starting position of the second endpoint of the edge.
    /// @param p_t1 The ending position of the point.
    /// @param e0_t1 The ending position of the first endpoint of the edge.
    /// @param e1_t1 The ending position of the second endpoint of the edge.
    /// @param toi The time of impact.
    /// @param min_distance The minimum distance between the point and the edge.
    /// @param tmax The maximum time to check for collision.
    /// @return True if a collision was detected, false otherwise.
    virtual bool point_edge_ccd(
        Eigen::ConstRef<VectorMax3d> p_t0,
        Eigen::ConstRef<VectorMax3d> e0_t0,
        Eigen::ConstRef<VectorMax3d> e1_t0,
        Eigen::ConstRef<VectorMax3d> p_t1,
        Eigen::ConstRef<VectorMax3d> e0_t1,
        Eigen::ConstRef<VectorMax3d> e1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;

    /// @brief Perform narrow phase CCD between a point and a linear triangle.
    /// @param p_t0 The starting position of the point.
    /// @param t0_t0 The starting position of the first vertex of the triangle.
    /// @param t1_t0 The starting position of the second vertex of the triangle.
    /// @param t2_t0 The starting position of the third vertex of the triangle.
    /// @param p_t1 The ending position of the point.
    /// @param t0_t1 The ending position of the first vertex of the triangle.
    /// @param t1_t1 The ending position of the second vertex of the triangle.
    /// @param t2_t1 The ending position of the third vertex of the triangle.
    /// @param toi The time of impact.
    /// @param min_distance The minimum distance between the point and the triangle.
    /// @param tmax The maximum time to check for collision.
    /// @return True if a collision was detected, false otherwise.
    virtual bool point_triangle_ccd(
        Eigen::ConstRef<Eigen::Vector3d> p_t0,
        Eigen::ConstRef<Eigen::Vector3d> t0_t0,
        Eigen::ConstRef<Eigen::Vector3d> t1_t0,
        Eigen::ConstRef<Eigen::Vector3d> t2_t0,
        Eigen::ConstRef<Eigen::Vector3d> p_t1,
        Eigen::ConstRef<Eigen::Vector3d> t0_t1,
        Eigen::ConstRef<Eigen::Vector3d> t1_t1,
        Eigen::ConstRef<Eigen::Vector3d> t2_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;

    /// @brief Perform narrow phase CCD between two linear edges.
    /// @param ea0_t0 The starting position of the first edge's first endpoint.
    /// @param ea1_t0 The starting position of the first edge's second endpoint.
    /// @param eb0_t0 The starting position of the second edge's first endpoint.
    /// @param eb1_t0 The starting position of the second edge's second endpoint.
    /// @param ea0_t1 The ending position of the first edge's first endpoint.
    /// @param ea1_t1 The ending position of the first edge's second endpoint.
    /// @param eb0_t1 The ending position of the second edge's first endpoint.
    /// @param eb1_t1 The ending position of the second edge's second endpoint.
    /// @param toi The time of impact.
    /// @param min_distance The minimum distance between the two edges.
    /// @param tmax The maximum time to check for collision.
    /// @return True if a collision was detected, false otherwise.
    virtual bool edge_edge_ccd(
        Eigen::ConstRef<Eigen::Vector3d> ea0_t0,
        Eigen::ConstRef<Eigen::Vector3d> ea1_t0,
        Eigen::ConstRef<Eigen::Vector3d> eb0_t0,
        Eigen::ConstRef<Eigen::Vector3d> eb1_t0,
        Eigen::ConstRef<Eigen::Vector3d> ea0_t1,
        Eigen::ConstRef<Eigen::Vector3d> ea1_t1,
        Eigen::ConstRef<Eigen::Vector3d> eb0_t1,
        Eigen::ConstRef<Eigen::Vector3d> eb1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const = 0;
};

} // namespace ipc