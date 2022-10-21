#include "../common.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_triangle_distance(py::module_& m)
{
    m.def(
        "point_triangle_distance",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
           const PointTriangleDistanceType dtype) {
            return point_triangle_distance(p, t0, t1, t2, dtype);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a two lines segments in 3D.

        Parameters:
            p : point.
            t0: first vertex of the triangle.
            t1: second vertex of the triangle.
            t2: third vertex of the triangle.
            dtype: (Optional) point-triangle distance type to compute

        Returns:
            The distance between the point and triangle.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dtype") = PointTriangleDistanceType::AUTO);

    m.def(
        "point_triangle_distance_gradient",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
           const PointTriangleDistanceType dtype) {
            Vector<double, 12> grad;
            point_triangle_distance_gradient(p, t0, t1, t2, grad, dtype);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines segments.

        Parameters:
            p : point
            t0: first vertex of the triangle
            t1: second vertex of the triangle
            t2: third vertex of the triangle
            dtype: (optional) point-triangle distance type to compute

        Returns:
            The gradient of the distance wrt ea0, ea1, eb0, and eb1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dtype") = PointTriangleDistanceType::AUTO);

    m.def(
        "point_triangle_distance_hessian",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
           const PointTriangleDistanceType dtype) {
            Eigen::Matrix<double, 12, 12> hess;
            point_triangle_distance_hessian(p, t0, t1, t2, hess, dtype);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines segments.

        Parameters:
            p : point
            t0: first vertex of the triangle
            t1: second vertex of the triangle
            t2: third vertex of the triangle
            dtype: (optional) point-triangle distance type to compute

        Returns:
            The hessian of the distance wrt ea0, ea1, eb0, and eb1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dtype") = PointTriangleDistanceType::AUTO);
}
