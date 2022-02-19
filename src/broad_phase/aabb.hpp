#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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

    /// @brief Compute a AABB for a moving point (i.e. temporal edge).
    static AABB from_point(const VectorMax3d& p, double inflation_radius = 0)
    {
        return AABB(p.array() - inflation_radius, p.array() + inflation_radius);
    }

    /// @brief Compute a AABB for a moving point (i.e. temporal edge).
    static AABB from_point(
        const VectorMax3d& p_t0,
        const VectorMax3d& p_t1,
        double inflation_radius = 0)
    {
        return AABB(
            from_point(p_t0, inflation_radius),
            from_point(p_t1, inflation_radius));
    }

    bool intersects(const AABB& other) const;

    ArrayMax3d min;
    ArrayMax3d max;
    std::array<int, 3> vertex_ids;
    // ArrayMax3d half_extent;
    // ArrayMax3d center;
};

void build_vertex_boxes(
    const Eigen::MatrixXd& V,
    std::vector<AABB>& vertex_boxes,
    double inflation_radius = 0);

void build_vertex_boxes(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<AABB>& vertex_boxes,
    double inflation_radius = 0);

void build_edge_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& E,
    std::vector<AABB>& edge_boxes);

void build_face_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& F,
    std::vector<AABB>& face_boxes);

} // namespace ipc
