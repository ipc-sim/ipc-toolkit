#include <common.hpp>

#include <ipc/tangent/relative_velocity.hpp>

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
        "dp0"_a, "dp1"_a);

    m.def(
        "point_point_relative_velocity_matrix",
        &point_point_relative_velocity_matrix,
        R"ipc_Qu8mg5v7(
        Compute the point-point relative velocity premultiplier matrix

        Parameters:
            dim: Dimension (2 or 3)

        Returns:
            The relative velocity premultiplier matrix
        )ipc_Qu8mg5v7",
        "dim"_a);

    m.def(
        "point_point_relative_velocity_matrix_jacobian",
        &point_point_relative_velocity_matrix_jacobian,
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the relative velocity premultiplier matrix

        Parameters:
            dim: Dimension (2 or 3)

        Returns:
            The Jacobian of the relative velocity premultiplier matrix
        )ipc_Qu8mg5v7",
        "dim"_a);

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
        "dp"_a, "de0"_a, "de1"_a, "alpha"_a);

    m.def(
        "point_edge_relative_velocity_matrix",
        &point_edge_relative_velocity_matrix,
        R"ipc_Qu8mg5v7(
        Compute the point-edge relative velocity premultiplier matrix

        Parameters:
            dim: Dimension (2 or 3)
            alpha: Parametric coordinate of the closest point on the edge

        Returns:
            The relative velocity premultiplier matrix
        )ipc_Qu8mg5v7",
        "dim"_a, "alpha"_a);

    m.def(
        "point_edge_relative_velocity_matrix_jacobian",
        &point_edge_relative_velocity_matrix_jacobian,
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the relative velocity premultiplier matrix

        Parameters:
            dim: Dimension (2 or 3)
            alpha: Parametric coordinate of the closest point on the edge

        Returns:
            The Jacobian of the relative velocity premultiplier matrix
        )ipc_Qu8mg5v7",
        "dim"_a, "alpha"_a);

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
        "dea0"_a, "dea1"_a, "deb0"_a, "deb1"_a, "coords"_a);

    m.def(
        "edge_edge_relative_velocity_matrix",
        &edge_edge_relative_velocity_matrix,
        R"ipc_Qu8mg5v7(
        Compute the edge-edge relative velocity matrix.

        Parameters:
            dim: Dimension (2 or 3)
            coords: Two parametric coordinates of the closest points on the edges

        Returns:
            The relative velocity matrix
        )ipc_Qu8mg5v7",
        "dim"_a, "coords"_a);

    m.def(
        "edge_edge_relative_velocity_matrix_jacobian",
        &edge_edge_relative_velocity_matrix_jacobian,
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the edge-edge relative velocity matrix.

        Parameters:
            dim: Dimension (2 or 3)
            coords: Two parametric coordinates of the closest points on the edges

        Returns:
            The Jacobian of the relative velocity matrix
        )ipc_Qu8mg5v7",
        "dim"_a, "coords"_a);

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
        "dp"_a, "dt0"_a, "dt1"_a, "dt2"_a, "coords"_a);

    m.def(
        "point_triangle_relative_velocity_matrix",
        &point_triangle_relative_velocity_matrix,
        R"ipc_Qu8mg5v7(
        Compute the point-triangle relative velocity matrix.

        Parameters:
            dim: Dimension (2 or 3)
            coords: Baricentric coordinates of the closest point on the triangle

        Returns:
            The relative velocity matrix
        )ipc_Qu8mg5v7",
        "dim"_a, "coords"_a);

    m.def(
        "point_triangle_relative_velocity_matrix_jacobian",
        &point_triangle_relative_velocity_matrix_jacobian,
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the point-triangle relative velocity matrix.

        Parameters:
            dim: Dimension (2 or 3)
            coords: Baricentric coordinates of the closest point on the triangle

        Returns:
            The Jacobian of the relative velocity matrix
        )ipc_Qu8mg5v7",
        "dim"_a, "coords"_a);
}
