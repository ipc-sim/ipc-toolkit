#include <common.hpp>

#include <ipc/tangent/closest_point.hpp>

using namespace ipc;

void define_closest_point(py::module_& m)
{
    m.def(
        "point_edge_closest_point",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_closest_point(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the barycentric coordinate of the closest point on the edge.

        Parameters:
            p: Point
            e0: First edge point
            e1: Second edge point

        Returns:
            barycentric coordinates of the closest point
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "point_edge_closest_point_jacobian",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_closest_point_jacobian(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the closest point on the edge.

        Parameters:
            p: Point
            e0: First edge point
            e1: Second edge point

        Returns:
            Jacobian of the closest point
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "edge_edge_closest_point",
        [](Eigen::ConstRef<Eigen::Vector3d> ea0,
           Eigen::ConstRef<Eigen::Vector3d> ea1,
           Eigen::ConstRef<Eigen::Vector3d> eb0,
           Eigen::ConstRef<Eigen::Vector3d> eb1) {
            return edge_edge_closest_point(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute the barycentric coordinates of the closest points between two edges.

        Parameters:
            ea0: First point of the first edge
            ea1: Second point of the first edge
            eb0: First point of the second edge
            eb1: Second point of the second edge

        Returns:
            Barycentric coordinates of the closest points
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "edge_edge_closest_point_jacobian",
        [](Eigen::ConstRef<Eigen::Vector3d> ea0,
           Eigen::ConstRef<Eigen::Vector3d> ea1,
           Eigen::ConstRef<Eigen::Vector3d> eb0,
           Eigen::ConstRef<Eigen::Vector3d> eb1) {
            return edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the closest points between two edges.

        Parameters:
            ea0: First point of the first edge
            ea1: Second point of the first edge
            eb0: First point of the second edge
            eb1: Second point of the second edge

        Returns:
            Jacobian of the closest points
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "point_triangle_closest_point",
        [](Eigen::ConstRef<Eigen::Vector3d> p,
           Eigen::ConstRef<Eigen::Vector3d> t0,
           Eigen::ConstRef<Eigen::Vector3d> t1,
           Eigen::ConstRef<Eigen::Vector3d> t2) {
            return point_triangle_closest_point(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Compute the barycentric coordinates of the closest point on the triangle.

        Parameters:
            p: Point
            t0: Triangle's first vertex
            t1: Triangle's second vertex
            t2: Triangle's third vertex

        Returns:
            Barycentric coordinates of the closest point
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);

    m.def(
        "point_triangle_closest_point_jacobian",
        [](Eigen::ConstRef<Eigen::Vector3d> p,
           Eigen::ConstRef<Eigen::Vector3d> t0,
           Eigen::ConstRef<Eigen::Vector3d> t1,
           Eigen::ConstRef<Eigen::Vector3d> t2) {
            return point_triangle_closest_point_jacobian(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the closest point on the triangle.

        Parameters:
            p: Point
            t0: Triangle's first vertex
            t1: Triangle's second vertex
            t2: Triangle's third vertex

        Returns:
            Jacobian of the closest point
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);
}
