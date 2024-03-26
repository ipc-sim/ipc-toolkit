#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <array>

namespace ipc {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    AABB() = default;

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

    /// @brief Construct an AABB for a static point.
    /// @param p The point's position.
    /// @param inflation_radius Radius of a sphere around the point which the AABB encloses.
    /// @return The constructed AABB.
    static AABB
    from_point(const VectorMax3d& p, const double inflation_radius = 0);

    /// @brief Construct an AABB for a moving point (i.e. temporal edge).
    /// @param p_t0 The point's position at time t=0.
    /// @param p_t1 The point's position at time t=1.
    /// @param inflation_radius Radius of a capsule around the temporal edge which the AABB encloses.
    /// @return The constructed AABB.
    static AABB from_point(
        const VectorMax3d& p_t0,
        const VectorMax3d& p_t1,
        const double inflation_radius = 0)
    {
        return AABB(
            from_point(p_t0, inflation_radius),
            from_point(p_t1, inflation_radius));
    }

    /// @brief Check if another AABB intersects with this one.
    /// @param other The other AABB.
    /// @return If the two AABBs intersect.
    bool intersects(const AABB& other) const;

    /// @brief Compute a conservative inflation of the AABB.
    static void conservative_inflation(
        ArrayMax3d& min, ArrayMax3d& max, const double inflation_radius);

public:
    /// @brief Minimum corner of the AABB.
    ArrayMax3d min;
    /// @brief Maximum corner of the AABB.
    ArrayMax3d max;
    /// @brief Vertex IDs attached to the AABB.
    std::array<long, 3> vertex_ids;
};

/// @brief Build one AABB per vertex position (row of V).
/// @param[in] vertices Vertex positions (rowwise).
/// @param[out] vertex_boxes Vertex AABBs.
/// @param[in] inflation_radius Radius of a sphere around the points which the AABBs enclose.
void build_vertex_boxes(
    const Eigen::MatrixXd& vertices,
    std::vector<AABB>& vertex_boxes,
    const double inflation_radius = 0);

/// @brief Build one AABB per vertex position moving linearly from t=0 to t=1.
/// @param vertices_t0 Vertex positions at t=0 (rowwise).
/// @param vertices_t1 Vertex positions at t=1 (rowwise).
/// @param vertex_boxes Vertex AABBs.
/// @param inflation_radius Radius of a capsule around the temporal edges which the AABBs enclose.
void build_vertex_boxes(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    std::vector<AABB>& vertex_boxes,
    const double inflation_radius = 0);

/// @brief Build one AABB per edge.
/// @param vertex_boxes Vertex AABBs.
/// @param edges Edges (rowwise).
/// @param edge_boxes Edge AABBs.
void build_edge_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& edges,
    std::vector<AABB>& edge_boxes);

/// @brief Build one AABB per face.
/// @param vertex_boxes Vertex AABBs.
/// @param faces Faces (rowwise).
/// @param face_boxes Face AABBs.
void build_face_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& faces,
    std::vector<AABB>& face_boxes);

} // namespace ipc
