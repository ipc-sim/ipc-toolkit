#include <common.hpp>

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_distance(py::module_& m)
{
    m.def(
        "edge_edge_distance", &edge_edge_distance,
        R"ipc_Qu8mg5v7(
        Compute the distance between a two lines segments in 3D.

        Parameters:
            ea0: first vertex of the first edge
            ea1: second vertex of the first edge
            eb0: first vertex of the second edge
            eb1: second vertex of the second edge
            dtype: (optional) edge-edge distance type to compute

        Returns:
            The distance between the two edges.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("dtype") = EdgeEdgeDistanceType::AUTO);

    m.def(
        "edge_edge_distance_gradient", &edge_edge_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines segments.

        Parameters:
            ea0: first vertex of the first edge
            ea1: second vertex of the first edge
            eb0: first vertex of the second edge
            eb1: second vertex of the second edge
            dtype: (optional) point edge distance type to compute

        Returns:
            The gradient of the distance wrt ea0, ea1, eb0, and eb1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("dtype") = EdgeEdgeDistanceType::AUTO);

    m.def(
        "edge_edge_distance_hessian", &edge_edge_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines segments.

        Parameters:
            ea0: first vertex of the first edge
            ea1: second vertex of the first edge
            eb0: first vertex of the second edge
            eb1: second vertex of the second edge
            dtype: (optional) point edge distance type to compute

        Returns:
            The hessian of the distance wrt ea0, ea1, eb0, and eb1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("dtype") = EdgeEdgeDistanceType::AUTO);
}
