#include <common.hpp>

#include <ipc/distance/signed/point_line.hpp>
#include <ipc/distance/signed/point_plane.hpp>
#include <ipc/distance/signed/line_line.hpp>

using namespace ipc;

void define_signed_distance(py::module_& m)
{
    // Point-line (2D) signed distance
    m.def(
        "point_line_signed_distance", point_line_signed_distance,
        R"ipc_Qu8mg5v7(
        Compute the signed distance from a point to a directed line segment (2D).

        Parameters:
            p: The query point (2D).
            e0: The first endpoint of the directed edge (2D).
            e1: The second endpoint of the directed edge (2D).

        Returns:
            The signed scalar distance from `p` to the (infinite) line through `e0` and `e1`.
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "point_line_signed_distance_gradient",
        point_line_signed_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the signed point-to-line distance (2D).

        Returns a 6-vector ordered as [dp, de0, de1].
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "point_line_signed_distance_hessian",
        point_line_signed_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the Hessian of the signed point-to-line distance (2D).

        Returns a 6x6 Hessian matrix with variables ordered as [p, e0, e1].
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    // Point-plane (3D) signed distance
    m.def(
        "point_plane_signed_distance", point_plane_signed_distance,
        R"ipc_Qu8mg5v7(
        Compute the signed distance from a point to the plane of a triangle (3D).

        Parameters:
            p: The query point (3D).
            t0: First vertex of the triangle (3D).
            t1: Second vertex of the triangle (3D).
            t2: Third vertex of the triangle (3D).

        Returns:
            The signed distance from `p` to the triangle plane.
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);

    m.def(
        "point_plane_signed_distance_gradient",
        point_plane_signed_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the signed point-to-plane distance (3D).

        Returns a 12-vector ordered as [dp, dt0, dt1, dt2].
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);

    m.def(
        "point_plane_signed_distance_hessian",
        point_plane_signed_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the Hessian of the signed point-to-plane distance (3D).

        Returns a 12x12 Hessian matrix with variables ordered as [p, t0, t1, t2].
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);

    // Line-line (3D) signed distance
    m.def(
        "line_line_signed_distance", line_line_signed_distance,
        R"ipc_Qu8mg5v7(
        Compute the signed distance between two lines in 3D.

        Parameters:
            ea0, ea1: Two points on the first line.
            eb0, eb1: Two points on the second line.

        Returns:
            The signed distance along the common normal between the two lines.
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "line_line_signed_distance_gradient",
        line_line_signed_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the signed line-line distance (3D).

        Returns a 12-vector ordered as [d/d(ea0); d/d(ea1); d/d(eb0); d/d(eb1)].
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "line_line_signed_distance_hessian", line_line_signed_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the Hessian of the signed line-line distance (3D).

        Returns a 12x12 Hessian matrix with variables ordered as [ea0, ea1, eb0, eb1].
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);
}
