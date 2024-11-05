#include "aabb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cfenv>

namespace ipc {

AABB::AABB(const ArrayMax3d& _min, const ArrayMax3d& _max)
    : min(_min)
    , max(_max)
{
    assert(min.size() == max.size());
    assert((min <= max).all());
    // half_extent = (max() - min()) / 2;
    // center = min() + half_extent();
}

AABB AABB::from_point(const VectorMax3d& p, const double inflation_radius)
{
    ArrayMax3d min = p.array(), max = p.array();
    conservative_inflation(min, max, inflation_radius);
    return AABB(min, max);
}

bool AABB::intersects(const AABB& other) const
{
    assert(this->min.size() == other.max.size());
    assert(this->max.size() == other.min.size());

    // NOTE: This is a faster check (https://gamedev.stackexchange.com/a/587),
    // but it is inexact and can result in false negatives.
    // return (std::abs(a.center.x() - b.center.x())
    //         <= (a.half_extent.x() + b.half_extent.x()))
    //     && (std::abs(a.center.y() - b.center.y())
    //         <= (a.half_extent.y() + b.half_extent.y()))
    //     && (a.min.size() == 2
    //         || std::abs(a.center.z() - b.center.z())
    //             <= (a.half_extent.z() + b.half_extent.z()));

    // This on the otherhand, is exact because there is no rounding.
    return (this->min <= other.max).all() && (other.min <= this->max).all();
}

void AABB::conservative_inflation(
    ArrayMax3d& min, ArrayMax3d& max, const double inflation_radius)
{
#pragma STDC FENV_ACCESS ON
    const int current_round = std::fegetround();

    std::fesetround(FE_DOWNWARD);
    min -= inflation_radius;

    std::fesetround(FE_UPWARD);
    max += inflation_radius;

    std::fesetround(current_round);
}

void build_vertex_boxes(
    const Eigen::MatrixXd& vertices,
    std::vector<AABB>& vertex_boxes,
    const double inflation_radius)
{
    vertex_boxes.resize(vertices.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, vertices.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                vertex_boxes[i] =
                    AABB::from_point(vertices.row(i), inflation_radius);
                vertex_boxes[i].vertex_ids = { { long(i), -1, -1 } };
            }
        });
}

void build_vertex_boxes(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    std::vector<AABB>& vertex_boxes,
    const double inflation_radius)
{
    vertex_boxes.resize(vertices_t0.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, vertices_t0.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                vertex_boxes[i] = AABB::from_point(
                    vertices_t0.row(i), vertices_t1.row(i), inflation_radius);
                vertex_boxes[i].vertex_ids = { { long(i), -1, -1 } };
            }
        });
}

void build_edge_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& edges,
    std::vector<AABB>& edge_boxes)
{
    edge_boxes.resize(edges.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, edges.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                edge_boxes[i] =
                    AABB(vertex_boxes[edges(i, 0)], vertex_boxes[edges(i, 1)]);
                edge_boxes[i].vertex_ids = { { edges(i, 0), edges(i, 1), -1 } };
            }
        });
}

void build_face_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& faces,
    std::vector<AABB>& face_boxes)
{
    face_boxes.resize(faces.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, faces.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                face_boxes[i] = AABB(
                    vertex_boxes[faces(i, 0)], vertex_boxes[faces(i, 1)],
                    vertex_boxes[faces(i, 2)]);
                face_boxes[i].vertex_ids = { { faces(i, 0), faces(i, 1),
                                               faces(i, 2) } };
            }
        });
}

} // namespace ipc