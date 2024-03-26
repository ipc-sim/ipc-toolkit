#include <common.hpp>

#include <ipc/distance/distance_type.hpp>

namespace py = pybind11;
using namespace ipc;

void define_distance_type(py::module_& m)
{
    py::enum_<PointEdgeDistanceType>(m, "PointEdgeDistanceType")
        .value(
            "P_E0", PointEdgeDistanceType::P_E0,
            "point is closest to edge vertex zero")
        .value(
            "P_E1", PointEdgeDistanceType::P_E1,
            "point is closest to edge vertex one")
        .value(
            "P_E", PointEdgeDistanceType::P_E,
            "point is closest to the interior of the edge")
        .value(
            "AUTO", PointEdgeDistanceType::AUTO,
            "automatically determine the closest point")
        .export_values();

    py::enum_<PointTriangleDistanceType>(m, "PointTriangleDistanceType")
        .value(
            "P_T0", PointTriangleDistanceType::P_T0,
            "point is closest to triangle vertex zero")
        .value(
            "P_T1", PointTriangleDistanceType::P_T1,
            "point is closest to triangle vertex one")
        .value(
            "P_T2", PointTriangleDistanceType::P_T2,
            "point is closest to triangle vertex two")
        .value(
            "P_E0", PointTriangleDistanceType::P_E0,
            "point is closest to triangle edge zero (vertex zero to one)")
        .value(
            "P_E1", PointTriangleDistanceType::P_E1,
            "point is closest to triangle edge one (vertex one to two)")
        .value(
            "P_E2", PointTriangleDistanceType::P_E2,
            "point is closest to triangle edge two (vertex two to zero)")
        .value(
            "P_T", PointTriangleDistanceType::P_T,
            "point is closest to the interior of the triangle")
        .value(
            "AUTO", PointTriangleDistanceType::AUTO,
            "automatically determine the closest point")
        .export_values();

    py::enum_<EdgeEdgeDistanceType>(m, "EdgeEdgeDistanceType")
        .value(
            "EA0_EB0", EdgeEdgeDistanceType::EA0_EB0,
            "edges are closest at vertex 0 of edge A and 0 of edge B")
        .value(
            "EA0_EB1", EdgeEdgeDistanceType::EA0_EB1,
            "edges are closest at vertex 0 of edge A and 1 of edge B")
        .value(
            "EA1_EB0", EdgeEdgeDistanceType::EA1_EB0,
            "edges are closest at vertex 1 of edge A and 0 of edge B")
        .value(
            "EA1_EB1", EdgeEdgeDistanceType::EA1_EB1,
            "edges are closest at vertex 1 of edge A and 1 of edge B")
        .value(
            "EA_EB0", EdgeEdgeDistanceType::EA_EB0,
            "edges are closest at the interior of edge A and vertex 0 "
            "of edge B")
        .value(
            "EA_EB1", EdgeEdgeDistanceType::EA_EB1,
            "edges are closest at the interior of edge A and vertex 1 "
            "of edge B")
        .value(
            "EA0_EB", EdgeEdgeDistanceType::EA0_EB,
            "edges are closest at vertex 0 of edge A and the interior "
            "of edge B")
        .value(
            "EA1_EB", EdgeEdgeDistanceType::EA1_EB,
            "edges are closest at vertex 1 of edge A and the interior "
            "of edge B")
        .value(
            "EA_EB", EdgeEdgeDistanceType::EA_EB,
            "edges are closest at an interior point of edge A and B")
        .value(
            "AUTO", EdgeEdgeDistanceType::AUTO,
            "automatically determine the closest point")
        .export_values();

    m.def(
        "point_edge_distance_type", &point_edge_distance_type,
        R"ipc_Qu8mg5v7(
        Determine the closest pair between a point and edge.

        Parameters:
            p: The point.
            e0: The first vertex of the edge.
            e1: The second vertex of the edge.

        Returns:
            The distance type of the point-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_triangle_distance_type", &point_triangle_distance_type,
        R"ipc_Qu8mg5v7(
        Determine the closest pair between a point and triangle.

        Parameters:
            p: The point.
            t0: The first vertex of the triangle.
            t1: The second vertex of the triangle.
            t2: The third vertex of the triangle.

        Returns:
            The distance type of the point-triangle pair.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "edge_edge_distance_type", &edge_edge_distance_type,
        R"ipc_Qu8mg5v7(
        Determine the closest pair between two edges.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The distance type of the edge-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_parallel_distance_type", &edge_edge_parallel_distance_type,
        R"ipc_Qu8mg5v7(
        Determine the closest pair between two parallel edges.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The distance type of the edge-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));
}
