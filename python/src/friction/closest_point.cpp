#include <common.hpp>

#include <ipc/friction/closest_point.hpp>

namespace py = pybind11;
using namespace ipc;

void define_closest_point(py::module_& m)
{
    m.def(
        "point_edge_closest_point",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_closest_point(p, e0, e1);
        },
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_edge_closest_point_jacobian",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_closest_point_jacobian(p, e0, e1);
        },
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "edge_edge_closest_point",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_closest_point(ea0, ea1, eb0, eb1);
        },
        "Compute the barycentric coordinates of the closest points",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_closest_point_jacobian",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1);
        },
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "point_triangle_closest_point",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_triangle_closest_point(p, t0, t1, t2);
        },
        "Compute the barycentric coordinates of the closest point on the triangle.",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "point_triangle_closest_point_jacobian",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_triangle_closest_point_jacobian(p, t0, t1, t2);
        },
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));
}
