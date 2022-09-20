#include "aabb.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {

AABB::AABB(const ArrayMax3d& min, const ArrayMax3d& max)
    : min(min)
    , max(max)
{
    assert(min.size() == max.size());
    assert((min <= max).all());
    // half_extent = (max() - min()) / 2;
    // center = min() + half_extent();
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
};

void build_vertex_boxes(
    const Eigen::MatrixXd& V,
    std::vector<AABB>& vertex_boxes,
    double inflation_radius)
{
    vertex_boxes.resize(V.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, V.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                vertex_boxes[i] = AABB::from_point(V.row(i), inflation_radius);
                vertex_boxes[i].vertex_ids = { { long(i), -1, -1 } };
            }
        });
}

void build_vertex_boxes(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<AABB>& vertex_boxes,
    double inflation_radius)
{
    vertex_boxes.resize(V0.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, V0.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                vertex_boxes[i] =
                    AABB::from_point(V0.row(i), V1.row(i), inflation_radius);
                vertex_boxes[i].vertex_ids = { { long(i), -1, -1 } };
            }
        });
}

void build_edge_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& E,
    std::vector<AABB>& edge_boxes)
{
    edge_boxes.resize(E.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, E.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                edge_boxes[i] =
                    AABB(vertex_boxes[E(i, 0)], vertex_boxes[E(i, 1)]);
                edge_boxes[i].vertex_ids = { { E(i, 0), E(i, 1), -1 } };
            }
        });
}

void build_face_boxes(
    const std::vector<AABB>& vertex_boxes,
    const Eigen::MatrixXi& F,
    std::vector<AABB>& face_boxes)
{
    face_boxes.resize(F.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, F.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                face_boxes[i] = AABB(
                    vertex_boxes[F(i, 0)], vertex_boxes[F(i, 1)],
                    vertex_boxes[F(i, 2)]);
                face_boxes[i].vertex_ids = { { F(i, 0), F(i, 1), F(i, 2) } };
            }
        });
}

} // namespace ipc