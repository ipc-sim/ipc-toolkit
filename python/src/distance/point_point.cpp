#include <common.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_point_distance(py::module_& m)
{
    m.def(
        "point_point_distance", &point_point_distance,
        R"ipc_Qu8mg5v7(
        Compute the distance between two points.

        Parameters:
            p0: first point
            p1: second point

        Returns:
            The distance between p0 and p1

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_point_distance_gradient", &point_point_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between two points.

        Parameters:
            p0: first point
            p1: second point

        Returns:
            The gradient of the distance wrt p0 and p1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_point_distance_hessian", &point_point_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and point.

        Parameters:
            p0: first point
            p1: second point

        Returns:
            The hessian of the distance wrt p0 and p1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p0"), py::arg("p1"));
}
