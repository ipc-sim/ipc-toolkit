#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <cstdint>

namespace ipc {

/// @brief How a narrow-phase CCD query is used.
enum class NarrowPhaseCCDMode : uint8_t {
    /// A standalone query: the method applies its own conservative strategy
    /// (minimum-distance inflation, small-TOI reruns, and rescaling).
    STANDALONE,
    /// The per-piece linear CCD inside a piecewise driver (e.g.,
    /// NonlinearCCD), which supplies its own conservative margins (the
    /// linearization error as the minimum separation, an early exit on
    /// decreased distance, and subdivision on zero TOI). The method skips its
    /// standalone strategy and returns the raw time of impact.
    COMPONENT,
    /// Like COMPONENT, but the caller cannot accept a time of impact of
    /// exactly zero (it knows the start of the piece is separated).
    COMPONENT_NO_ZERO_TOI,
};

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
    /// @param mode How the query is used (see NarrowPhaseCCDMode).
    /// @return True if a collision was detected, false otherwise.
    virtual bool point_point_ccd(
        Eigen::ConstRef<VectorMax3d> p0_t0,
        Eigen::ConstRef<VectorMax3d> p1_t0,
        Eigen::ConstRef<VectorMax3d> p0_t1,
        Eigen::ConstRef<VectorMax3d> p1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const = 0;

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
    /// @param mode How the query is used (see NarrowPhaseCCDMode).
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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const = 0;

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
    /// @param mode How the query is used (see NarrowPhaseCCDMode).
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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const = 0;

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
    /// @param mode How the query is used (see NarrowPhaseCCDMode).
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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const = 0;

    /// @brief The conservative rescaling value used to avoid taking steps exactly to impact.
    ///
    /// This is used to creat a effective minimum distance between objects,
    /// which can help prevent numerical issues when objects are very close to
    /// each other.
    ///
    /// I.e., min_effective_distance = (1.0 - conservative_rescaling) *
    /// (initial_distance - min_distance)
    ///
    /// If the time of impact is small (e.g. 1e-6), then the CCD will be rerun
    /// without a effective minimum distance to try to get a more accurate time
    /// of impact. In this case, the conservative rescaling will be used to
    /// scale the time of impact to avoid taking steps exactly to impact.
    double conservative_rescaling;
};

} // namespace ipc