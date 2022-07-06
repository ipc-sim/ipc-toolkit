#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/friction/closest_point.hpp>

namespace py = pybind11;
using namespace ipc;

void define_closest_point_members(py::module_& m)
{
    m.def(
        "point_edge_closest_point",
        &point_edge_closest_point<VectorMax3d, VectorMax3d, VectorMax3d>, "",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_edge_closest_point",
        &edge_edge_closest_point<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "Compute the barycentric coordinates of the closest points",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "point_triangle_closest_point",
        &point_triangle_closest_point<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "Compute the barycentric coordinates of the closest point on the "
        "triangle.",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "point_edge_closest_point_jacobian",
        &point_edge_closest_point_jacobian<
            VectorMax3d, VectorMax3d, VectorMax3d>,
        "", py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_edge_closest_point_jacobian",
        &edge_edge_closest_point_jacobian<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "", py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "point_triangle_closest_point_jacobian",
        &point_triangle_closest_point_jacobian<
            Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>,
        "", py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));
}
