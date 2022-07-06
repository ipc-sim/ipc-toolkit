#include "../common.hpp"

#include <ipc/distance/line_line.hpp>

namespace py = pybind11;
using namespace ipc;

void define_line_line_distance(py::module_& m)
{
    m.def(
        "line_line_distance",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return line_line_distance(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a two infinite lines in 3D.

        Parameters:
            ea0: first vertex of the edge defining the first line
            ea1: second vertex of the edge defining the first line
            eb0: first vertex of the edge defining the second line
            eb1: second vertex of the edge defining the second line

        Returns:
            The distance between the two lines.

        Note:
            The distance is actually squared distance.

        Warning:
            If the lines are parallel this function returns a distance of zero.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "line_line_distance_gradient",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            Vector<double, 12> grad;
            line_line_distance_gradient(ea0, ea1, eb0, eb1, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines in 3D.

        Parameters:
            ea0: first vertex of the edge defining the first line
            ea1: second vertex of the edge defining the first line
            eb0: first vertex of the edge defining the second line
            eb1: second vertex of the edge defining the second line

        Returns:
            The gradient of the distance wrt ea0, ea1, eb0, and eb1.

        Note:
            The distance is actually squared distance.

        Warning:
            If the lines are parallel this function returns a distance of zero.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "line_line_distance_hessian",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            Eigen::Matrix<double, 12, 12> hess;
            line_line_distance_hessian(ea0, ea1, eb0, eb1, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines in 3D.

        Parameters:
            ea0: first vertex of the edge defining the first line
            ea1: second vertex of the edge defining the first line
            eb0: first vertex of the edge defining the second line
            eb1: second vertex of the edge defining the second line

        Returns:
            The hessian of the distance wrt ea0, ea1, eb0, and eb1.

        Note:
            The distance is actually squared distance.

        Warning:
            If the lines are parallel this function returns a distance of zero.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));
}
