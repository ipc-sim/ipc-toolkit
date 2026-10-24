#include <common.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

using namespace ipc;

void define_point_point_distance(py::module_& m)
{
    m.def(
        "point_point_distance",
        [](Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1) {
            return point_point_distance(p0, p1);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between two points.

        Note:
            The distance is actually squared distance.

        Parameters:
            p0: The first point.
            p1: The second point.

        Returns:
            The distance between p0 and p1.
        )ipc_Qu8mg5v7",
        "p0"_a, "p1"_a);

    m.def(
        "point_point_distance_gradient",
        [](Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1) {
            return point_point_distance_gradient(p0, p1);
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between two points.

        Note:
            The distance is actually squared distance.

        Parameters:
            p0: The first point.
            p1: The second point.

        Returns:
            The computed gradient.
        )ipc_Qu8mg5v7",
        "p0"_a, "p1"_a);

    m.def(
        "point_point_distance_hessian",
        [](Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1) {
            return point_point_distance_hessian(p0, p1);
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between two points.

        Note:
            The distance is actually squared distance.

        Parameters:
            p0: The first point.
            p1: The second point.

        Returns:
            The computed hessian.
        )ipc_Qu8mg5v7",
        "p0"_a, "p1"_a);
}
