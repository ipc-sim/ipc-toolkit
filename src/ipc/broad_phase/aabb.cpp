#include "aabb.hpp"

#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cfenv>

namespace ipc {

AABB::AABB(Eigen::ConstRef<ArrayMax3d> _min, Eigen::ConstRef<ArrayMax3d> _max)
    : min(Eigen::Array3d::Zero())
    , max(Eigen::Array3d::Zero())
{
    assert(_min.size() == _max.size());
    assert((_min <= _max).all());
    min.head(_min.size()) = _min;
    max.head(_max.size()) = _max;
}

AABB AABB::from_point(
    Eigen::ConstRef<VectorMax3d> p, const double inflation_radius)
{
    ArrayMax3d min = p.array(), max = p.array();
    conservative_inflation(min, max, inflation_radius);
    return AABB(min, max);
}

bool AABB::intersects(const AABB& other) const
{
    // This is exact because there is no rounding.
    return (this->min <= other.max).all() && (other.min <= this->max).all();
}

void AABB::conservative_inflation(
    Eigen::Ref<ArrayMax3d> min,
    Eigen::Ref<ArrayMax3d> max,
    const double inflation_radius)
{
    assert(min.size() == max.size());
    assert((min <= max).all());
    assert(inflation_radius >= 0);

    // Nudge the bounds outward to ensure conservativity.

    min = min.unaryExpr([inflation_radius](double v) {
        return std::nextafter(
            v - inflation_radius, -std::numeric_limits<double>::infinity());
    });

    max = max.unaryExpr([inflation_radius](double v) {
        return std::nextafter(
            v + inflation_radius, std::numeric_limits<double>::infinity());
    });
}

void build_vertex_boxes(
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    AABBs& vertex_boxes,
    const double inflation_radius)
{
    IPC_TOOLKIT_PROFILE_BLOCK("build_vertex_boxes");

    vertex_boxes.resize(vertices.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, vertices.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                vertex_boxes[i] =
                    AABB::from_point(vertices.row(i), inflation_radius);
                vertex_boxes[i].vertex_ids = { { index_t(i), -1, -1 } };
            }
        });
}

void build_vertex_boxes(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    AABBs& vertex_boxes,
    const double inflation_radius)
{
    IPC_TOOLKIT_PROFILE_BLOCK("build_vertex_boxes");

    vertex_boxes.resize(vertices_t0.rows());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, vertices_t0.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                vertex_boxes[i] = AABB::from_point(
                    vertices_t0.row(i), vertices_t1.row(i), inflation_radius);
                vertex_boxes[i].vertex_ids = { { index_t(i), -1, -1 } };
            }
        });
}

void build_edge_boxes(
    const AABBs& vertex_boxes,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    AABBs& edge_boxes)
{
    IPC_TOOLKIT_PROFILE_BLOCK("build_edge_boxes");

    if (edge_boxes.size() != edges.rows()) {
        IPC_TOOLKIT_PROFILE_BLOCK("resize_edge_boxes");
        edge_boxes.resize(edges.rows());
    }

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, edges.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const int e0 = edges(i, 0), e1 = edges(i, 1);
                edge_boxes[i] = AABB(vertex_boxes[e0], vertex_boxes[e1]);
                edge_boxes[i].vertex_ids = { { e0, e1, -1 } };
            }
        });
}

void build_face_boxes(
    const AABBs& vertex_boxes,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    AABBs& face_boxes)
{
    IPC_TOOLKIT_PROFILE_BLOCK("build_face_boxes");

    if (face_boxes.size() != faces.rows()) {
        IPC_TOOLKIT_PROFILE_BLOCK("resize_face_boxes");
        face_boxes.resize(faces.rows());
    }

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, faces.rows()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const int f0 = faces(i, 0), f1 = faces(i, 1), f2 = faces(i, 2);
                face_boxes[i] =
                    AABB(vertex_boxes[f0], vertex_boxes[f1], vertex_boxes[f2]);
                face_boxes[i].vertex_ids = { { f0, f1, f2 } };
            }
        });
}

} // namespace ipc