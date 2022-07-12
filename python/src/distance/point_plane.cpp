#include "../common.hpp"

#include <ipc/distance/point_plane.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_plane_distance(py::module_& m)
{
    m.def(
        "point_plane_distance",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_plane_distance(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and a plane.

        Parameters:
            p: point
            t0: first vertex of the triangle
            t1: second vertex of the triangle
            t2: third vertex of the triangle

        Returns:
            The distance between the point and plane.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t1"));

    m.def(
        "point_plane_distance_gradient",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            Vector<double, 12> grad;
            point_plane_distance_gradient(p, t0, t1, t2, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and a plane.

        Parameters:
            p: point
            t0: first vertex of the triangle
            t1: second vertex of the triangle
            t2: third vertex of the triangle

        Returns:
            The gradient of the distance wrt p, t0, t1, and t2.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t1"));

    m.def(
        "point_plane_distance_hessian",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            Eigen::Matrix<double, 12, 12> hess;
            point_plane_distance_hessian(p, t0, t1, t2, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and a plane.

        Parameters:
            p: point
            t0: first vertex of the triangle
            t1: second vertex of the triangle
            t2: third vertex of the triangle

        Returns:
            The hessian of the distance wrt p, t0, t1, and t2.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t1"));
}
