#include <common.hpp>

#include <ipc/utils/intersection.hpp>

#include <igl/predicates/segment_segment_intersect.h>

using namespace ipc;

void define_intersection(py::module_& m)
{
    m.def(
        "is_edge_intersecting_triangle", &is_edge_intersecting_triangle, "e0"_a,
        "e1"_a, "t0"_a, "t1"_a, "t2"_a);

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
