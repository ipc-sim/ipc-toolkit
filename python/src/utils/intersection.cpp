#include <common.hpp>

#include <ipc/utils/intersection.hpp>

#include <igl/predicates/segment_segment_intersect.h>

namespace py = pybind11;
using namespace ipc;

void define_intersection(py::module_& m)
{
    m.def(
        "is_edge_intersecting_triangle", &is_edge_intersecting_triangle,
        py::arg("e0"), py::arg("e1"), py::arg("t0"), py::arg("t1"),
        py::arg("t2"));

    m.def(
        "segment_segment_intersect",
        [](const Eigen::Vector2d& A, const Eigen::Vector2d& B,
           const Eigen::Vector2d& C, const Eigen::Vector2d& D) -> bool {
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
        py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"));
}
