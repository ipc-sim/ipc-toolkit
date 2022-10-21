#include "../common.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_edge_distance(py::module_& m)
{
    m.def(
        "point_edge_distance",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1,
           const PointEdgeDistanceType dtype) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_distance(p, e0, e1, dtype);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and edge in 2D or 3D.

        Parameters:
            p: point
            e0: first vertex of the edge
            e1: second vertex of the edge
            dtype: (optional) point edge distance type to compute

        Returns:
            The distance between the point and edge

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = PointEdgeDistanceType::AUTO);

    m.def(
        "point_edge_distance_gradient",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1,
           const PointEdgeDistanceType dtype) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            VectorMax9d grad;
            point_edge_distance_gradient(p, e0, e1, grad, dtype);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and edge.

        Parameters:
            p: point
            e0: first vertex of the edge
            e1: second vertex of the edge
            dtype: (optional) point edge distance type to compute

        Returns:
            The gradient of the distance wrt p, e0, and e1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = PointEdgeDistanceType::AUTO);

    m.def(
        "point_edge_distance_hessian",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1,
           const PointEdgeDistanceType dtype) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            MatrixMax9d hess;
            point_edge_distance_hessian(p, e0, e1, hess, dtype);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and edge.

        Parameters:
            p: point
            e0: first vertex of the edge
            e1: second vertex of the edge
            dtype: (optional) point edge distance type to compute

        Returns:
            The hessian of the distance wrt p, e0, and e1.

        Note:
            The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = PointEdgeDistanceType::AUTO);
}
