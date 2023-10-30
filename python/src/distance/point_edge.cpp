#include <common.hpp>

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_edge_distance(py::module_& m)
{
    m.def(
        "point_edge_distance", &point_edge_distance,
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and edge in 2D or 3D.

        Note:
            The distance is actually squared distance.

        Parameters:
            p: The point.
            e0: The first vertex of the edge.
            e1: The second vertex of the edge.
            dtype: The point edge distance type to compute.

        Returns:
            The distance between the point and edge.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = PointEdgeDistanceType::AUTO);

    m.def(
        "point_edge_distance_gradient", &point_edge_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and edge.

        Note:
            The distance is actually squared distance.

        Parameters:
            p: The point.
            e0: The first vertex of the edge.
            e1: The second vertex of the edge.
            dtype: The point edge distance type to compute.

        Returns:
            grad The gradient of the distance wrt p, e0, and e1.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = PointEdgeDistanceType::AUTO);

    m.def(
        "point_edge_distance_hessian", &point_edge_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and edge.

        Note:
            The distance is actually squared distance.

        Parameters:
            p: The point.
            e0: The first vertex of the edge.
            e1: The second vertex of the edge.
            dtype: The point edge distance type to compute.

        Returns:
            hess The hessian of the distance wrt p, e0, and e1.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = PointEdgeDistanceType::AUTO);
}
