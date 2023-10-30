#include <common.hpp>

#include <ipc/friction/tangent_basis.hpp>

namespace py = pybind11;
using namespace ipc;

void define_tangent_basis(py::module_& m)
{
    m.def(
        "point_point_tangent_basis",
        [](const VectorMax3d& p0, const VectorMax3d& p1) {
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
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_point_tangent_basis_jacobian",
        [](const VectorMax3d& p0, const VectorMax3d& p1) {
            assert_2D_or_3D_vector(p0, "p0");
            assert_2D_or_3D_vector(p1, "p1");
            return point_point_tangent_basis_jacobian(p0, p1);
        },
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_edge_tangent_basis",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_tangent_basis(p, e0, e1);
        },
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_edge_tangent_basis_jacobian",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_tangent_basis_jacobian(p, e0, e1);
        },
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_edge_tangent_basis",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_tangent_basis(ea0, ea1, eb0, eb1);
        },
        "Compute a basis for the space tangent to the edge-edge pair.",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_tangent_basis_jacobian",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_tangent_basis_jacobian(ea0, ea1, eb0, eb1);
        },
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "point_triangle_tangent_basis",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
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
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "point_triangle_tangent_basis_jacobian",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_triangle_tangent_basis_jacobian(p, t0, t1, t2);
        },
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));
}
