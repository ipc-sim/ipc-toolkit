#include <common.hpp>

#include <ipc/tangent/tangent_basis.hpp>

using namespace ipc;

void define_tangent_basis(py::module_& m)
{
    m.def(
        "point_point_tangent_basis",
        [](Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1) {
            assert_2D_or_3D_vector(p0, "p0");
            assert_2D_or_3D_vector(p1, "p1");
            return point_point_tangent_basis(p0, p1);
        },
        R"ipc_Qu8mg5v7(
        Compute a basis for the space tangent to the point-point pair.

        Parameters:
            p0: First point
            p1: Second point

        Returns:
            A 3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "p0"_a, "p1"_a);

    m.def(
        "point_point_tangent_basis_jacobian",
        [](Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1) {
            assert_2D_or_3D_vector(p0, "p0");
            assert_2D_or_3D_vector(p1, "p1");
            return point_point_tangent_basis_jacobian(p0, p1);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the tangent basis for the point-point pair.

        Parameters:
            p0: First point
            p1: Second point

        Returns:
            A 6*3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "p0"_a, "p1"_a);

    m.def(
        "point_edge_tangent_basis",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_tangent_basis(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute a basis for the space tangent to the point-edge pair.

        Parameters:
            p: Point
            e0: First edge point
            e1: Second edge point

        Returns:
            A 3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "point_edge_tangent_basis_jacobian",
        [](Eigen::ConstRef<VectorMax3d> p, Eigen::ConstRef<VectorMax3d> e0,
           Eigen::ConstRef<VectorMax3d> e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_tangent_basis_jacobian(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the tangent basis for the point-edge pair.

        Parameters:
            p: Point
            e0: First edge point
            e1: Second edge point

        Returns:
            A 9*3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "p"_a, "e0"_a, "e1"_a);

    m.def(
        "edge_edge_tangent_basis",
        [](Eigen::ConstRef<Eigen::Vector3d> ea0,
           Eigen::ConstRef<Eigen::Vector3d> ea1,
           Eigen::ConstRef<Eigen::Vector3d> eb0,
           Eigen::ConstRef<Eigen::Vector3d> eb1) {
            return edge_edge_tangent_basis(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute a basis for the space tangent to the edge-edge pair.

        Parameters:
            ea0: First point of the first edge
            ea1: Second point of the first edge
            eb0: First point of the second edge
            eb1: Second point of the second edge

        Returns:
            A 3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "edge_edge_tangent_basis_jacobian",
        [](Eigen::ConstRef<Eigen::Vector3d> ea0,
           Eigen::ConstRef<Eigen::Vector3d> ea1,
           Eigen::ConstRef<Eigen::Vector3d> eb0,
           Eigen::ConstRef<Eigen::Vector3d> eb1) {
            return edge_edge_tangent_basis_jacobian(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the tangent basis for the edge-edge pair.

        Parameters:
            ea0: First point of the first edge
            ea1: Second point of the first edge
            eb0: First point of the second edge
            eb1: Second point of the second edge

        Returns:
            A 12*3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a);

    m.def(
        "point_triangle_tangent_basis",
        [](Eigen::ConstRef<Eigen::Vector3d> p,
           Eigen::ConstRef<Eigen::Vector3d> t0,
           Eigen::ConstRef<Eigen::Vector3d> t1,
           Eigen::ConstRef<Eigen::Vector3d> t2) {
            return point_triangle_tangent_basis(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Compute a basis for the space tangent to the point-triangle pair.

        .. math::

            \begin{bmatrix}
            \frac{t_1 - t_0}{\|t_1 - t_0\|} & \frac{((t_1 - t_0)\times(t_2 - t_0))
            \times(t_1 - t_0)}{\|((t_1 - t_0)\times(t_2 - t_0))\times(t_1 - t_0)\|}
            \end{bmatrix}

        Parameters:
            p: Point
            t0: Triangle's first vertex
            t1: Triangle's second vertex
            t2: Triangle's third vertex

        Returns:
            A 3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);

    m.def(
        "point_triangle_tangent_basis_jacobian",
        [](Eigen::ConstRef<Eigen::Vector3d> p,
           Eigen::ConstRef<Eigen::Vector3d> t0,
           Eigen::ConstRef<Eigen::Vector3d> t1,
           Eigen::ConstRef<Eigen::Vector3d> t2) {
            return point_triangle_tangent_basis_jacobian(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Compute the Jacobian of the tangent basis for the point-triangle pair.

        Parameters:
            p: Point
            t0: Triangle's first vertex
            t1: Triangle's second vertex
            t2: Triangle's third vertex

        Returns:
            A 12*3x2 matrix whose columns are the basis vectors.
        )ipc_Qu8mg5v7",
        "p"_a, "t0"_a, "t1"_a, "t2"_a);
}
