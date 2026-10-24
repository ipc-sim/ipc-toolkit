#include <common.hpp>

#include <ipc/distance/point_line.hpp>

using namespace ipc;

void define_point_line_distance(py::module_& m)
{
    m.def(
        "point_line_distance",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            return point_line_distance(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and line in 2D or 3D.

        Note:
            The distance is actually squared distance.

        Parameters:
            p: The point.
            e0: The first vertex of the edge defining the line.
            e1: The second vertex of the edge defining the line.

        Returns:
            The distance between the point and line.
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "point_line_distance_gradient",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            return point_line_distance_gradient(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and line.

        Note:
            The distance is actually squared distance.

        Parameters:
            p: The point.
            e0: The first vertex of the edge defining the line.
            e1: The second vertex of the edge defining the line.

        Returns:
            The gradient of the distance wrt p, e0, and e1.
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "point_line_distance_hessian",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            return point_line_distance_hessian(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and line.

        Note:
            The distance is actually squared distance.

        Parameters:
            p: The point.
            e0: The first vertex of the edge defining the line.
            e1: The second vertex of the edge defining the line.

        Returns:
            The hessian of the distance wrt p, e0, and e1.
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);
}
