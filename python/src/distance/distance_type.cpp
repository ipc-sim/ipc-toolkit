#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <ipc/distance/distance_type.hpp>

#include "../utils.hpp"

namespace py = pybind11;
using namespace ipc;

void define_distance_type_functions(py::module_& m)
{
    py::enum_<PointEdgeDistanceType>(m, "PointEdgeDistanceType")
        .value("P_E0", PointEdgeDistanceType::P_E0)
        .value("P_E1", PointEdgeDistanceType::P_E1)
        .value("P_E", PointEdgeDistanceType::P_E)
        .export_values();

    py::enum_<PointTriangleDistanceType>(m, "PointTriangleDistanceType")
        .value("P_T0", PointTriangleDistanceType::P_T0)
        .value("P_T1", PointTriangleDistanceType::P_T1)
        .value("P_T2", PointTriangleDistanceType::P_T2)
        .value("P_E0", PointTriangleDistanceType::P_E0)
        .value("P_E1", PointTriangleDistanceType::P_E1)
        .value("P_E2", PointTriangleDistanceType::P_E2)
        .value("P_T", PointTriangleDistanceType::P_T)
        .export_values();

    py::enum_<EdgeEdgeDistanceType>(m, "EdgeEdgeDistanceType")
        .value("EA0_EB0", EdgeEdgeDistanceType::EA0_EB0)
        .value("EA0_EB1", EdgeEdgeDistanceType::EA0_EB1)
        .value("EA1_EB0", EdgeEdgeDistanceType::EA1_EB0)
        .value("EA1_EB1", EdgeEdgeDistanceType::EA1_EB1)
        .value("EA_EB0", EdgeEdgeDistanceType::EA_EB0)
        .value("EA_EB1", EdgeEdgeDistanceType::EA_EB1)
        .value("EA0_EB", EdgeEdgeDistanceType::EA0_EB)
        .value("EA1_EB", EdgeEdgeDistanceType::EA1_EB)
        .value("EA_EB", EdgeEdgeDistanceType::EA_EB)
        .export_values();

    m.def(
        "point_edge_distance_type",
        [](const VectorMax3d& p, const VectorMax3d& e0, const VectorMax3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_distance_type(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Determine the closest pair between a point and edge.

        Parameters
        ----------
        p  : The point
        e0 : The first vertex of the edge
        e1 : The second vertex of the edge

        Returns
        -------
        The distance type of the point-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_triangle_distance_type",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_triangle_distance_type(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Determine the closest pair between a point and triangle.

        Parameters
        ----------
        p  : The point
        t0 : The first vertex of the triangle
        t1 : The second vertex of the triangle
        t2 : The third vertex of the triangle

        Returns
        -------
        The distance type of the point-triangle pair.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "edge_edge_distance_type",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_distance_type(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Determine the closest pair between two edges.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge

        Returns
        -------
        The distance type of the edge-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));
}