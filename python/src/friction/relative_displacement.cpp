#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/friction/relative_displacement.hpp>

namespace py = pybind11;
using namespace ipc;

void define_relative_displacement_members(py::module_& m)
{
    m.def(
        "point_point_relative_displacement",
        &point_point_relative_displacement<VectorMax3d, VectorMax3d>, "",
        py::arg("dp0"), py::arg("dp1"));

    m.def(
        "point_point_relative_displacement_matrix",
        &point_point_relative_displacement_matrix<double>, "", py::arg("dim"));

    m.def(
        "point_point_relative_displacement_matrix_jacobian",
        &point_point_relative_displacement_matrix_jacobian<double>, "",
        py::arg("dim"));

    m.def(
        "point_edge_relative_displacement",
        &point_edge_relative_displacement<
            VectorMax3d, VectorMax3d, VectorMax3d, double>,
        "", py::arg("dp"), py::arg("de0"), py::arg("de1"), py::arg("alpha"));

    m.def(
        "point_edge_relative_displacement_matrix",
        &point_edge_relative_displacement_matrix<double>, "", py::arg("dim"),
        py::arg("alpha"));

    m.def(
        "point_edge_relative_displacement_matrix_jacobian",
        &point_edge_relative_displacement_matrix_jacobian<double>, "",
        py::arg("dim"), py::arg("alpha"));

    m.def(
        "edge_edge_relative_displacement",
        &edge_edge_relative_displacement<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d,
            Eigen::Vector2d>,
        "Compute the relative displacement of the edges.", py::arg("dea0"),
        py::arg("dea1"), py::arg("deb0"), py::arg("deb1"), py::arg("coords"));

    m.def(
        "edge_edge_relative_displacement_matrix",
        &edge_edge_relative_displacement_matrix<Eigen::Vector2d>, "",
        py::arg("dim"), py::arg("coords"));

    m.def(
        "edge_edge_relative_displacement_matrix_jacobian",
        &edge_edge_relative_displacement_matrix_jacobian<Eigen::Vector2d>, "",
        py::arg("dim"), py::arg("coords"));

    m.def(
        "point_triangle_relative_displacement",
        &point_triangle_relative_displacement<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d,
            Eigen::Vector2d>,
        "Compute the relative displacement of the point to the triangle.",
        py::arg("dp"), py::arg("dt0"), py::arg("dt1"), py::arg("dt2"),
        py::arg("coords"));

    m.def(
        "point_triangle_relative_displacement_matrix",
        &point_triangle_relative_displacement_matrix<Eigen::Vector2d>, "",
        py::arg("dim"), py::arg("coords"));

    m.def(
        "point_triangle_relative_displacement_matrix_jacobian",
        &point_triangle_relative_displacement_matrix_jacobian<Eigen::Vector2d>,
        "", py::arg("dim"), py::arg("coords"));
}
