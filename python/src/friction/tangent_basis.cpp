#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/friction/tangent_basis.hpp>

namespace py = pybind11;
using namespace ipc;

void define_tangent_basis_members(py::module_& m)
{
    m.def(
        "point_point_tangent_basis",
        &point_point_tangent_basis<VectorMax3d, VectorMax3d>, "", py::arg("p0"),
        py::arg("p1"));

    m.def(
        "point_edge_tangent_basis",
        &point_edge_tangent_basis<VectorMax3d, VectorMax3d, VectorMax3d>, "",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_edge_tangent_basis",
        &edge_edge_tangent_basis<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "Compute a basis for the space tangent to the edge-edge pair.",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "point_triangle_tangent_basis",
        &point_triangle_tangent_basis<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "Compute a basis for the space tangent to the point-triangle pair.",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "point_point_tangent_basis_jacobian",
        &point_point_tangent_basis_jacobian<VectorMax3d, VectorMax3d>, "",
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_edge_tangent_basis_jacobian",
        &point_edge_tangent_basis_jacobian<
            VectorMax3d, VectorMax3d, VectorMax3d>,
        "", py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_edge_tangent_basis_jacobian",
        &edge_edge_tangent_basis_jacobian<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "", py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "point_triangle_tangent_basis_jacobian",
        &point_triangle_tangent_basis_jacobian<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "", py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));
}
