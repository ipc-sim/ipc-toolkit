#include <common.hpp>

#include <ipc/friction/relative_velocity.hpp>

namespace py = pybind11;
using namespace ipc;

void define_relative_velocity(py::module_& m)
{
    m.def(
        "point_point_relative_velocity", &point_point_relative_velocity,
        R"ipc_Qu8mg5v7(
        Compute the relative velocity of two points

        Parameters:
            dp0: Velocity of the first point
            dp1: Velocity of the second point

        Returns:
            The relative velocity of the two points
        )ipc_Qu8mg5v7",
        py::arg("dp0"), py::arg("dp1"));

    m.def(
        "point_point_relative_velocity_matrix",
        &point_point_relative_velocity_matrix,
        R"ipc_Qu8mg5v7(
        Compute the relative velocity premultiplier matrix

        Parameters:
            dim: Dimension (2 or 3)

        Returns:
            The relative velocity premultiplier matrix
        )ipc_Qu8mg5v7",
        py::arg("dim"));

    m.def(
        "point_point_relative_velocity_matrix_jacobian",
        &point_point_relative_velocity_matrix_jacobian,
        R"ipc_Qu8mg5v7(
        Compute the jacobian of the relative velocity premultiplier matrix

        Parameters:
            dim: Dimension (2 or 3)

        Returns:
            The jacobian of the relative velocity premultiplier matrix
        )ipc_Qu8mg5v7",
        py::arg("dim"));

    m.def(
        "point_edge_relative_velocity", &point_edge_relative_velocity,
        R"ipc_Qu8mg5v7(
        Compute the relative velocity of a point and an edge

        Parameters:
            dp: Velocity of the point
            de0: Velocity of the first endpoint of the edge
            de1: Velocity of the second endpoint of the edge
            alpha: Parametric coordinate of the closest point on the edge

        Returns:
            The relative velocity of the point and the edge
        )ipc_Qu8mg5v7",
        py::arg("dp"), py::arg("de0"), py::arg("de1"), py::arg("alpha"));

    m.def(
        "point_edge_relative_velocity_matrix",
        &point_edge_relative_velocity_matrix, py::arg("dim"), py::arg("alpha"));

    m.def(
        "point_edge_relative_velocity_matrix_jacobian",
        &point_edge_relative_velocity_matrix_jacobian, py::arg("dim"),
        py::arg("alpha"));

    m.def(
        "edge_edge_relative_velocity", &edge_edge_relative_velocity,
        R"ipc_Qu8mg5v7(
        Compute the relative velocity of the edges.

        Parameters:
            dea0: Velocity of the first endpoint of the first edge
            dea1: Velocity of the second endpoint of the first edge
            deb0: Velocity of the first endpoint of the second edge
            deb1: Velocity of the second endpoint of the second edge
            coords: Two parametric coordinates of the closest points on the edges

        Returns:
            The relative velocity of the edges
        )ipc_Qu8mg5v7",
        py::arg("dea0"), py::arg("dea1"), py::arg("deb0"), py::arg("deb1"),
        py::arg("coords"));

    m.def(
        "edge_edge_relative_velocity_matrix",
        &edge_edge_relative_velocity_matrix, py::arg("dim"), py::arg("coords"));

    m.def(
        "edge_edge_relative_velocity_matrix_jacobian",
        &edge_edge_relative_velocity_matrix_jacobian, py::arg("dim"),
        py::arg("coords"));

    m.def(
        "point_triangle_relative_velocity", &point_triangle_relative_velocity,
        R"ipc_Qu8mg5v7(
        Compute the relative velocity of the point to the triangle.

        Parameters:
            dp: Velocity of the point
            dt0: Velocity of the first vertex of the triangle
            dt1: Velocity of the second vertex of the triangle
            dt2: Velocity of the third vertex of the triangle
            coords: Baricentric coordinates of the closest point on the triangle

        Returns:
            The relative velocity of the point to the triangle
        )ipc_Qu8mg5v7",
        py::arg("dp"), py::arg("dt0"), py::arg("dt1"), py::arg("dt2"),
        py::arg("coords"));

    m.def(
        "point_triangle_relative_velocity_matrix",
        &point_triangle_relative_velocity_matrix, py::arg("dim"),
        py::arg("coords"));

    m.def(
        "point_triangle_relative_velocity_matrix_jacobian",
        &point_triangle_relative_velocity_matrix_jacobian, py::arg("dim"),
        py::arg("coords"));
}
