#include "../common.hpp"

#include <ipc/friction/relative_displacement.hpp>

namespace py = pybind11;
using namespace ipc;

void define_relative_displacement(py::module_& m)
{
    m.def(
        "point_point_relative_displacement",
        [](const VectorMax3d& dp0, const VectorMax3d& dp1) {
            assert_2D_or_3D_vector(dp0, "dp0");
            assert_2D_or_3D_vector(dp1, "dp1");
            return point_point_relative_displacement(dp0, dp1);
        },
        "", py::arg("dp0"), py::arg("dp1"));

    m.def(
        "point_point_relative_displacement_matrix",
        &point_point_relative_displacement_matrix<double>, "", py::arg("dim"));

    m.def(
        "point_point_relative_displacement_matrix_jacobian",
        &point_point_relative_displacement_matrix_jacobian<double>, "",
        py::arg("dim"));

    m.def(
        "point_edge_relative_displacement",
        [](const VectorMax3d& dp, const VectorMax3d& de0,
           const VectorMax3d& de1, const double alpha) {
            assert_2D_or_3D_vector(dp, "dp");
            assert_2D_or_3D_vector(de0, "de0");
            assert_2D_or_3D_vector(de1, "de1");
            return point_edge_relative_displacement(dp, de0, de1, alpha);
        },
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
        [](const Eigen::Vector3d& dea0, const Eigen::Vector3d& dea1,
           const Eigen::Vector3d& deb0, const Eigen::Vector3d& deb1,
           const Eigen::Vector2d& coords) {
            return edge_edge_relative_displacement(
                dea0, dea1, deb0, deb1, coords);
        },
        "Compute the relative displacement of the edges.", py::arg("dea0"),
        py::arg("dea1"), py::arg("deb0"), py::arg("deb1"), py::arg("coords"));

    m.def(
        "edge_edge_relative_displacement_matrix",
        [](const int dim, const Eigen::Vector2d& coords) {
            return edge_edge_relative_displacement_matrix(dim, coords);
        },
        "", py::arg("dim"), py::arg("coords"));

    m.def(
        "edge_edge_relative_displacement_matrix_jacobian",
        [](const int dim, const Eigen::Vector2d& coords) {
            return edge_edge_relative_displacement_matrix_jacobian(dim, coords);
        },
        "", py::arg("dim"), py::arg("coords"));

    m.def(
        "point_triangle_relative_displacement",
        [](const Eigen::Vector3d& dp, const Eigen::Vector3d& dt0,
           const Eigen::Vector3d& dt1, const Eigen::Vector3d& dt2,
           const Eigen::Vector2d& coords) {
            return point_triangle_relative_displacement(
                dp, dt0, dt1, dt2, coords);
        },
        "Compute the relative displacement of the point to the triangle.",
        py::arg("dp"), py::arg("dt0"), py::arg("dt1"), py::arg("dt2"),
        py::arg("coords"));

    m.def(
        "point_triangle_relative_displacement_matrix",
        [](const int dim, const Eigen::Vector2d& coords) {
            return point_triangle_relative_displacement_matrix(dim, coords);
        },
        "", py::arg("dim"), py::arg("coords"));

    m.def(
        "point_triangle_relative_displacement_matrix_jacobian",
        [](const int dim, const Eigen::Vector2d& coords) {
            return point_triangle_relative_displacement_matrix_jacobian(
                dim, coords);
        },
        "", py::arg("dim"), py::arg("coords"));
}
