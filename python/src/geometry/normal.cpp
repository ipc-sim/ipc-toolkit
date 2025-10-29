#include <common.hpp>

#include <ipc/geometry/normal.hpp>

using namespace ipc;

void define_normal(py::module_& m)
{
    m.def(
        "normalization_and_jacobian", &normalization_and_jacobian,
        R"ipc_qu8mg5v7(
        Computes the normalization and Jacobian of a vector.

        Parameters
        ----------
        x: The input vector.

        Returns
        -------
        A tuple containing the normalized vector and its Jacobian.
        )ipc_qu8mg5v7",
        py::arg("x"));

    m.def(
        "normalization_jacobian", &normalization_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the normalization operation.

        Parameters
        ----------
        x: The input vector.

        Returns
        -------
        The Jacobian of the normalization operation.
        )ipc_qu8mg5v7",
        py::arg("x"));

    m.def(
        "normalization_and_jacobian_and_hessian",
        &normalization_and_jacobian_and_hessian,
        R"ipc_qu8mg5v7(
        Computes the normalization, Jacobian, and Hessian of a vector.

        Parameters
        ----------
        x: The input vector.

        Returns
        -------
        A tuple containing the normalized vector, its Jacobian, and its Hessian.
        )ipc_qu8mg5v7",
        py::arg("x"));

    m.def(
        "normalization_hessian", &normalization_hessian,
        R"ipc_qu8mg5v7(
        Computes the Hessian of the normalization operation.

        Parameters
        ----------
        x: The input vector.

        Returns
        -------
        The Hessian of the normalization operation.
        )ipc_qu8mg5v7",
        py::arg("x"));

    m.def(
        "cross_product_matrix", &cross_product_matrix,
        R"ipc_qu8mg5v7(
        Cross product matrix for 3D vectors.

        Parameters
        ----------
        v: Vector to create the cross product matrix for.

        Returns
        -------
        The cross product matrix of the vector.
        )ipc_qu8mg5v7",
        py::arg("v"));

    m.def(
        "cross_product_matrix_jacobian", &cross_product_matrix_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the cross product matrix.
        Returns
        -------
        The Jacobian of the cross product matrix.
        )ipc_qu8mg5v7");

    m.def(
        "edge_vertex_unnormalized_normal", &edge_vertex_unnormalized_normal,
        R"ipc_qu8mg5v7(
        Computes the unnormalized normal vector of an edge-vertex pair.

        Parameters
        ----------
        v: The vertex position.
        e0: The start position of the edge.
        e1: The end position of the edge.

        Returns
        -------
        The unnormalized normal vector.
        )ipc_qu8mg5v7",
        py::arg("v"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_vertex_normal", &edge_vertex_normal,
        R"ipc_qu8mg5v7(
        Computes the normal vector of an edge-vertex pair.

        Parameters
        ----------
        v: The vertex position.
        e0: The start position of the edge.
        e1: The end position of the edge.

        Returns
        -------
        The normal vector.
        )ipc_qu8mg5v7",
        py::arg("v"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_vertex_unnormalized_normal_jacobian",
        &edge_vertex_unnormalized_normal_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the unnormalized normal vector of an edge-vertex pair.

        Parameters
        ----------
        v: The vertex position.
        e0: The start position of the edge.
        e1: The end position of the edge.

        Returns
        -------
        The Jacobian of the unnormalized normal vector.
        )ipc_qu8mg5v7",
        py::arg("v"), py::arg("e0"), py::arg("e1"));

    m.def(
        "triangle_unnormalized_normal", &triangle_unnormalized_normal,
        R"ipc_qu8mg5v7(
        Computes the unnormalized normal vector of a triangle.

        Parameters
        ----------
        a: The first vertex of the triangle.
        b: The second vertex of the triangle.
        c: The third vertex of the triangle.

        Returns
        -------
        The unnormalized normal vector of the triangle.
        )ipc_qu8mg5v7",
        py::arg("a"), py::arg("b"), py::arg("c"));

    m.def(
        "triangle_normal", &triangle_normal,
        R"ipc_qu8mg5v7(
        Computes the normal vector of a triangle.

        Parameters
        ----------
        a: The first vertex of the triangle.
        b: The second vertex of the triangle.
        c: The third vertex of the triangle.

        Returns
        -------
        The normal vector of the triangle.
        )ipc_qu8mg5v7",
        py::arg("a"), py::arg("b"), py::arg("c"));

    m.def(
        "triangle_unnormalized_normal_jacobian",
        &triangle_unnormalized_normal_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the unnormalized normal vector of a triangle.

        Parameters
        ----------
        a: The first vertex of the triangle.
        b: The second vertex of the triangle.
        c: The third vertex of the triangle.

        Returns
        -------
        The Jacobian of the unnormalized normal vector of the triangle.
        )ipc_qu8mg5v7",
        py::arg("a"), py::arg("b"), py::arg("c"));

    m.def(
        "triangle_unnormalized_normal_hessian",
        &triangle_unnormalized_normal_hessian,
        R"ipc_qu8mg5v7(
        Computes the Hessian of the unnormalized normal vector of a triangle.

        Parameters
        ----------
        a: The first vertex of the triangle.
        b: The second vertex of the triangle.
        c: The third vertex of the triangle.

        Returns
        -------
        The Hessian of the unnormalized normal vector of the triangle.
        )ipc_qu8mg5v7",
        py::arg("a"), py::arg("b"), py::arg("c"));

    m.def(
        "triangle_normal_jacobian", &triangle_normal_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the normal vector of a triangle.

        Parameters
        ----------
        a: The first vertex of the triangle.
        b: The second vertex of the triangle.
        c: The third vertex of the triangle.

        Returns
        -------
        The Jacobian of the normal vector of the triangle.
        )ipc_qu8mg5v7",
        py::arg("a"), py::arg("b"), py::arg("c"));

    m.def(
        "edge_edge_unnormalized_normal", &edge_edge_unnormalized_normal,
        R"ipc_qu8mg5v7(
        Computes the unnormalized normal vector of two edges.

        Parameters
        ----------
        ea0: The first vertex of the first edge.
        ea1: The second vertex of the first edge.
        eb0: The first vertex of the second edge.
        eb1: The second vertex of the second edge.

        Returns
        -------
        The unnormalized normal vector of the two edges.
        )ipc_qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_normal", &edge_edge_normal,
        R"ipc_qu8mg5v7(
        Computes the normal vector of two edges.

        Parameters
        ----------
        ea0: The first vertex of the first edge.
        ea1: The second vertex of the first edge.
        eb0: The first vertex of the second edge.
        eb1: The second vertex of the second edge.

        Returns
        -------
        The normal vector of the two edges.
        )ipc_qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_unnormalized_normal_jacobian",
        &edge_edge_unnormalized_normal_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the unnormalized normal vector of two edges.

        Parameters
        ----------
        ea0: The first vertex of the first edge.
        ea1: The second vertex of the first edge.
        eb0: The first vertex of the second edge.
        eb1: The second vertex of the second edge.

        Returns
        -------
        The Jacobian of the unnormalized normal vector of the two edges.
        )ipc_qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_normal_jacobian", &edge_edge_normal_jacobian,
        R"ipc_qu8mg5v7(
        Computes the Jacobian of the normal vector of two edges.

        Parameters
        ----------
        ea0: The first vertex of the first edge.
        ea1: The second vertex of the first edge.
        eb0: The first vertex of the second edge.
        eb1: The second vertex of the second edge.

        Returns
        -------
        The Jacobian of the normal vector of the two edges.
        )ipc_qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));
}
