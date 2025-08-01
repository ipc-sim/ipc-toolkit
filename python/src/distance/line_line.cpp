#include <common.hpp>

#include <ipc/distance/line_line.hpp>

using namespace ipc;

void define_line_line_distance(py::module_& m)
{
    m.def(
        "line_line_distance", &line_line_distance,
        R"ipc_Qu8mg5v7(
        Compute the distance between a two infinite lines in 3D.

        Note:
            The distance is actually squared distance.

        Warning:
            If the lines are parallel this function returns a distance of zero.

        Parameters:
            ea0: The first vertex of the edge defining the first line.
            ea1: The second vertex of the edge defining the first line.
            eb0: The first vertex of the edge defining the second line.
            eb1: The second vertex of the edge defining the second line.

        Returns:
            The distance between the two lines.
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "line_line_distance_gradient", &line_line_distance_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines in 3D.

        Note:
            The distance is actually squared distance.

        Warning:
            If the lines are parallel this function returns a distance of zero.

        Parameters:
            ea0: The first vertex of the edge defining the first line.
            ea1: The second vertex of the edge defining the first line.
            eb0: The first vertex of the edge defining the second line.
            eb1: The second vertex of the edge defining the second line.

        Returns:
            The gradient of the distance wrt ea0, ea1, eb0, and eb1.
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "line_line_distance_hessian", &line_line_distance_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines in 3D.

        Note:
            The distance is actually squared distance.

        Warning:
            If the lines are parallel this function returns a distance of zero.

        Parameters:
            ea0: The first vertex of the edge defining the first line.
            ea1: The second vertex of the edge defining the first line.
            eb0: The first vertex of the edge defining the second line.
            eb1: The second vertex of the edge defining the second line.

        Returns:
            The hessian of the distance wrt ea0, ea1, eb0, and eb1.
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);
}
