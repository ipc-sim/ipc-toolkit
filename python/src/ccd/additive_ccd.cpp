#include <common.hpp>

#include <ipc/ccd/additive_ccd.hpp>

namespace py = pybind11;
using namespace ipc;

void define_additive_ccd(py::module_& m)
{
    using namespace ipc::additive_ccd;

    auto m_accd = m.def_submodule(
        "additive_ccd", "Additive CCD method of [Li et al. 2021].");

    m_accd.def(
        "point_point_ccd",
        [](const VectorMax3d& p0_t0, const VectorMax3d& p1_t0,
           const VectorMax3d& p0_t1, const VectorMax3d& p1_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_point_ccd(
                p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance, tmax,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between two points using continuous collision detection.

        Parameters:
            p0_t0: The initial position of the first point.
            p1_t0: The initial position of the second point.
            p0_t1: The final position of the first point.
            p1_t1: The final position of the second point.
            min_distance: The minimum distance between two objects.
            tmax: The maximum time to check for collisions.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the two points.
        )ipc_Qu8mg5v7",
        py::arg("p0_t0"), py::arg("p1_t0"), py::arg("p0_t1"), py::arg("p1_t1"),
        py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m_accd.def(
        "point_edge_ccd",
        [](const VectorMax3d& p_t0, const VectorMax3d& e0_t0,
           const VectorMax3d& e1_t0, const VectorMax3d& p_t1,
           const VectorMax3d& e0_t1, const VectorMax3d& e1_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_edge_ccd(
                p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi, min_distance, tmax,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between a point and an edge using continuous collision detection.

        Parameters:
            p_t0: The initial position of the point.
            e0_t0: The initial position of the first endpoint of the edge.
            e1_t0: The initial position of the second endpoint of the edge.
            p_t1: The final position of the point.
            e0_t1: The final position of the first endpoint of the edge.
            e1_t1: The final position of the second endpoint of the edge.
            min_distance: The minimum distance between two objects.
            tmax: The maximum time to check for collisions.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the point and the edge.
        )ipc_Qu8mg5v7",
        py::arg("p_t0"), py::arg("e0_t0"), py::arg("e1_t0"), py::arg("p_t1"),
        py::arg("e0_t1"), py::arg("e1_t1"), py::arg("min_distance") = 0.0,
        py::arg("tmax") = 1.0,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m_accd.def(
        "point_triangle_ccd",
        [](const Eigen::Vector3d& p_t0, const Eigen::Vector3d& t0_t0,
           const Eigen::Vector3d& t1_t0, const Eigen::Vector3d& t2_t0,
           const Eigen::Vector3d& p_t1, const Eigen::Vector3d& t0_t1,
           const Eigen::Vector3d& t1_t1, const Eigen::Vector3d& t2_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_triangle_ccd(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi,
                min_distance, tmax, conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between a point and a triangle using continuous collision detection.

        Parameters:
            p_t0: The initial position of the point.
            t0_t0: The initial position of the first vertex of the triangle.
            t1_t0: The initial position of the second vertex of the triangle.
            t2_t0: The initial position of the third vertex of the triangle.
            p_t1: The final position of the point.
            t0_t1: The final position of the first vertex of the triangle.
            t1_t1: The final position of the second vertex of the triangle.
            t2_t1: The final position of the third vertex of the triangle.
            min_distance: The minimum distance between two objects.
            tmax: The maximum time to check for collisions.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the point and the triangle.
        )ipc_Qu8mg5v7",
        py::arg("p_t0"), py::arg("t0_t0"), py::arg("t1_t0"), py::arg("t2_t0"),
        py::arg("p_t1"), py::arg("t0_t1"), py::arg("t1_t1"), py::arg("t2_t1"),
        py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m_accd.def(
        "edge_edge_ccd",
        [](const Eigen::Vector3d& ea0_t0, const Eigen::Vector3d& ea1_t0,
           const Eigen::Vector3d& eb0_t0, const Eigen::Vector3d& eb1_t0,
           const Eigen::Vector3d& ea0_t1, const Eigen::Vector3d& ea1_t1,
           const Eigen::Vector3d& eb0_t1, const Eigen::Vector3d& eb1_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = edge_edge_ccd(
                ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
                toi, min_distance, tmax, conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between two edges using continuous collision detection.

        Parameters:
            ea0_t0: The initial position of the first endpoint of the first edge.
            ea1_t0: The initial position of the second endpoint of the first edge.
            eb0_t0: The initial position of the first endpoint of the second edge.
            eb1_t0: The initial position of the second endpoint of the second edge.
            ea0_t1: The final position of the first endpoint of the first edge.
            ea1_t1: The final position of the second endpoint of the first edge.
            eb0_t1: The final position of the first endpoint of the second edge.
            eb1_t1: The final position of the second endpoint of the second edge.
            min_distance: The minimum distance between two objects.
            tmax: The maximum time to check for collisions.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the two edges.
        )ipc_Qu8mg5v7",
        py::arg("ea0_t0"), py::arg("ea1_t0"), py::arg("eb0_t0"),
        py::arg("eb1_t0"), py::arg("ea0_t1"), py::arg("ea1_t1"),
        py::arg("eb0_t1"), py::arg("eb1_t1"), py::arg("min_distance") = 0.0,
        py::arg("tmax") = 1.0,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m_accd.def(
        "additive_ccd",
        [](VectorMax12d x, const VectorMax12d& dx,
           const std::function<double(const VectorMax12d&)>& distance_squared,
           const double max_disp_mag, const double min_distance = 0.0,
           const double tmax = 1.0,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = ipc::additive_ccd::additive_ccd(
                x, dx, distance_squared, max_disp_mag, toi, min_distance, tmax,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between two objects using additive continuous collision detection.

        Parameters:
            distance_squared: A function that computes the squared distance between the two objects at a given time.
            min_distance: The minimum distance between the objects.
            tmax: The maximum time to check for collisions.
            conservative_rescaling: The amount to rescale the objects by to ensure conservative advancement.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the two objects.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("dx"), py::arg("distance_squared"),
        py::arg("max_disp_mag"), py::arg("min_distance") = 0.0,
        py::arg("tmax") = 1.0,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);
}
