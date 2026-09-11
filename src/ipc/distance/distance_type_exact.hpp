#pragma once

#include <ipc/config.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

constexpr double PARALLEL_THRESHOLD {
    2.5e-16
}; // TODO set to zero eventually (requires geogram)

/// @brief Runtime switch between the standard analytic distance-type routines
/// and the exact predicate-based implementations. Controls
/// point_edge_distance_type_exact, point_triangle_distance_type_exact, and
/// edge_edge_distance_type_exact. Defaults to predicate (exact).
class DistanceTypeConfig {
public:
    static DistanceTypeConfig& instance()
    {
        static DistanceTypeConfig cfg;
        return cfg;
    }

    bool use_standard() const { return m_use_standard; }
    void set_use_standard(bool v) { m_use_standard = v; }

    DistanceTypeConfig(const DistanceTypeConfig&) = delete;
    DistanceTypeConfig& operator=(const DistanceTypeConfig&) = delete;

private:
    DistanceTypeConfig() = default;
    bool m_use_standard = false;
};

/// @brief Determine the closest pair between a point and edge, using exact
/// (geogram) predicates when available and falling back to the standard
/// analytic implementation otherwise (or when DistanceTypeConfig prefers it).
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The distance type of the point-edge pair.
PointEdgeDistanceType point_edge_distance_type_exact(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1);

/// @brief Determine the closest pair between a point and triangle, using
/// exact (geogram) predicates when available.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance type of the point-triangle pair.
PointTriangleDistanceType point_triangle_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

/// @brief Determine the closest pair between two edges, using exact
/// (geogram) predicates when available.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
EdgeEdgeDistanceType edge_edge_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

/// @brief Determine whether two edges are (nearly) parallel.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
bool is_parallel_edge_edge(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

} // namespace ipc
