#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <array>

namespace ipc {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    AABB() { }

    AABB(const ArrayMax3d& min, const ArrayMax3d& max);

    AABB(const AABB& aabb1, const AABB& aabb2)
        : AABB(aabb1.min.min(aabb2.min), aabb1.max.max(aabb2.max))
    {
    }

    AABB(const AABB& aabb1, const AABB& aabb2, const AABB& aabb3)
        : AABB(
            aabb1.min.min(aabb2.min).min(aabb3.min),
            aabb1.max.max(aabb2.max).max(aabb3.max))
    {
    }

    /// @brief Compute a AABB for a static point.
    static AABB
    from_point(const VectorMax3d& p, const double inflation_radius = 0);

    /// @brief Compute a AABB for a moving point (i.e. temporal edge).
    static AABB from_point(
        const VectorMax3d& p_t0,
        const VectorMax3d& p_t1,
        const double inflation_radius = 0)
    {
        return AABB(
            from_point(p_t0, inflation_radius),
            from_point(p_t1, inflation_radius));
    }

    bool intersects(const AABB& other) const;

    /// @brief Compute a conservative inflation of the AABB.
    static void conservative_inflation(
        ArrayMax3d& min, ArrayMax3d& max, const double inflation_radius);

public:
    ArrayMax3d min;
    ArrayMax3d max;
    std::array<long, 3> vertex_ids;
};

void build_vertex_boxes(
    const Eigen::MatrixXd& vertices,
    std::vector<AABB>& vertex_boxes,
    const double inflation_radius = 0);

void build_vertex_boxes(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    std::vector<AABB>& vertex_boxes,
    const double inflation_radius = 0);

void build_edge_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& edges,
    std::vector<AABB>& edge_boxes);

void build_face_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& faces,
    std::vector<AABB>& face_boxes);

} // namespace ipc
