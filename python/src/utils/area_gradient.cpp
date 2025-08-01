#include <common.hpp>

#include <ipc/utils/area_gradient.hpp>
#include <ipc/utils/eigen_ext.hpp>

using namespace ipc;

void define_area_gradient(py::module_& m)
{
    m.def(
        "edge_length_gradient", &edge_length_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of an edge's length.

        Parameters:
            e0: The first vertex of the edge.
            e1: The second vertex of the edge.

        Returns:
            The gradient of the edge's length wrt e0, and e1.
        )ipc_Qu8mg5v7",
        "e0"_a, "e1"_a);

    m.def(
        "triangle_area_gradient", &triangle_area_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the area of a triangle.

        Parameters:
            t0: The first vertex of the triangle.
            t1: The second vertex of the triangle.
            t2: The third vertex of the triangle.

        Returns:
            The gradient of the triangle's area t0, t1, and t2.
        )ipc_Qu8mg5v7",
        "t0"_a, "t1"_a, "t2"_a);
}
