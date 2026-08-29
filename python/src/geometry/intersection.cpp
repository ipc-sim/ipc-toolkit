#include <common.hpp>

#include <ipc/geometry/intersection.hpp>

#include <igl/predicates/segment_segment_intersect.h>

using namespace ipc;

void define_intersection(py::module_& m)
{
    m.def(
        "is_edge_intersecting_triangle", &is_edge_intersecting_triangle, "e0"_a,
        "e1"_a, "t0"_a, "t1"_a, "t2"_a);

    m.def(
        "edge_triangle_intersection",
        [](Eigen::ConstRef<Eigen::Vector3d> e0,
           Eigen::ConstRef<Eigen::Vector3d> e1,
           Eigen::ConstRef<Eigen::Vector3d> t0,
           Eigen::ConstRef<Eigen::Vector3d> t1,
           Eigen::ConstRef<Eigen::Vector3d> t2) {
            double u, v, t;
            const bool hit =
                edge_triangle_intersection(e0, e1, t0, t1, t2, u, v, t);
            return std::make_tuple(hit, u, v, t);
        },
        R"ipc_Qu8mg5v7(
        Edge-triangle intersection test that also returns the hit location.

        Uses the same robust orient3d plane-side gate as
        is_edge_intersecting_triangle, then solves for the barycentric (u, v) on
        the triangle and the parameter t along the edge. Boundary hits count as
        intersections.

        Note:
            The coordinates are only meaningful when intersects is True, and are
            NaN otherwise. The one exception is a degenerate configuration (edge
            coplanar with the triangle, or a degenerate triangle), where this
            conservatively reports True but the coordinates are not uniquely
            defined and are left as NaN.

        Parameters:
            e0: Edge start point.
            e1: Edge end point.
            t0: Triangle vertex 0.
            t1: Triangle vertex 1.
            t2: Triangle vertex 2.

        Returns:
            Tuple of (intersects, u, v, t) where (u, v) are the triangle
            barycentric coordinates (along t1-t0, t2-t0) and t is the edge
            parameter in [0, 1].
        )ipc_Qu8mg5v7",
        "e0"_a, "e1"_a, "t0"_a, "t1"_a, "t2"_a);

    m.def(
        "segment_segment_intersect",
        [](Eigen::ConstRef<Eigen::Vector2d> A,
           Eigen::ConstRef<Eigen::Vector2d> B,
           Eigen::ConstRef<Eigen::Vector2d> C,
           Eigen::ConstRef<Eigen::Vector2d> D) -> bool {
            igl::predicates::exactinit();
            return igl::predicates::segment_segment_intersect(A, B, C, D);
        },
        R"ipc_Qu8mg5v7(
        Given two segments in 2d test whether they intersect each other using predicates orient2d

        Parameters:
            A: 1st endpoint of segment 1
            B: 2st endpoint of segment 1
            C: 1st endpoint of segment 2
            D: 2st endpoint of segment 2

        Returns:
            true if they intersect
        )ipc_Qu8mg5v7",
        "A"_a, "B"_a, "C"_a, "D"_a);
}
