#include <common.hpp>

#include <ipc/distance/point_line.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_line_distance(py::module_& m)
{
    m.def(
        "point_line_distance", &point_line_distance,
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and line in 2D or 3D.

        Parameters:
            p: point
            e0: first vertex of the edge defining the line
            e1: second vertex of the edge defining the line

        Returns:
            The distance between the point and line.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_line_distance_gradient", &point_line_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and line.

        Parameters:
            p: point
            e0: first vertex of the edge defining the line.
            e1: second vertex of the edge defining the line.

        Returns:
            The gradient of the distance wrt p, e0, and e1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_line_distance_hessian", &point_line_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and line.

        Parameters:
            p: point
            e0: first vertex of the edge defining the line
            e1: second vertex of the edge defining the line

        Returns:
            The hessian of the distance wrt p, e0, and e1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));
}
