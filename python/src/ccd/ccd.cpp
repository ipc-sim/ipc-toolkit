#include <common.hpp>

#include <ipc/ccd/ccd.hpp>

namespace py = pybind11;
using namespace ipc;

void define_ccd(py::module_& m)
{
    m.def(
        "point_point_ccd",
        [](const VectorMax3d& p0_t0, const VectorMax3d& p1_t0,
           const VectorMax3d& p0_t1, const VectorMax3d& p1_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double tolerance = DEFAULT_CCD_TOLERANCE,
           const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_point_ccd(
                p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance, tmax, tolerance,
                max_iterations, conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between two points using continuous collision detection.

        Parameters:
            p0_t0: The initial position of the first point.
            p1_t0: The initial position of the second point.
            p0_t1: The final position of the first point.
            p1_t1: The final position of the second point.
            min_distance: The minimum distance between the points.
            tmax: The maximum time to check for collisions.
            tolerance: The error tolerance for the time of impact.
            max_iterations: The maximum number of iterations to perform.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the two points.
        )ipc_Qu8mg5v7",
        py::arg("p0_t0"), py::arg("p1_t0"), py::arg("p0_t1"), py::arg("p1_t1"),
        py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
        py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m.def(
        "point_edge_ccd",
        [](const VectorMax3d& p_t0, const VectorMax3d& e0_t0,
           const VectorMax3d& e1_t0, const VectorMax3d& p_t1,
           const VectorMax3d& e0_t1, const VectorMax3d& e1_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double tolerance = DEFAULT_CCD_TOLERANCE,
           const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_edge_ccd(
                p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi, min_distance, tmax,
                tolerance, max_iterations, conservative_rescaling);
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
            min_distance: The minimum distance between the objects.
            tmax: The maximum time to check for collisions.
            tolerance: The error tolerance for the time of impact.
            max_iterations: The maximum number of iterations to perform.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the point and the edge.
        )ipc_Qu8mg5v7",
        py::arg("p_t0"), py::arg("e0_t0"), py::arg("e1_t0"), py::arg("p_t1"),
        py::arg("e0_t1"), py::arg("e1_t1"), py::arg("min_distance") = 0.0,
        py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m.def(
        "edge_edge_ccd",
        [](const Eigen::Vector3d& ea0_t0, const Eigen::Vector3d& ea1_t0,
           const Eigen::Vector3d& eb0_t0, const Eigen::Vector3d& eb1_t0,
           const Eigen::Vector3d& ea0_t1, const Eigen::Vector3d& ea1_t1,
           const Eigen::Vector3d& eb0_t1, const Eigen::Vector3d& eb1_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double tolerance = DEFAULT_CCD_TOLERANCE,
           const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = edge_edge_ccd(
                ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
                toi, min_distance, tmax, tolerance, max_iterations,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between two edges in 3D using continuous collision detection.

        Parameters:
            ea0_t0: The initial position of the first endpoint of the first edge.
            ea1_t0: The initial position of the second endpoint of the first edge.
            eb0_t0: The initial position of the first endpoint of the second edge.
            eb1_t0: The initial position of the second endpoint of the second edge.
            ea0_t1: The final position of the first endpoint of the first edge.
            ea1_t1: The final position of the second endpoint of the first edge.
            eb0_t1: The final position of the first endpoint of the second edge.
            eb1_t1: The final position of the second endpoint of the second edge.
            min_distance: The minimum distance between the objects.
            tmax: The maximum time to check for collisions.
            tolerance: The error tolerance for the time of impact.
            max_iterations: The maximum number of iterations to perform.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the two edges.
        )ipc_Qu8mg5v7",
        py::arg("ea0_t0"), py::arg("ea1_t0"), py::arg("eb0_t0"),
        py::arg("eb1_t0"), py::arg("ea0_t1"), py::arg("ea1_t1"),
        py::arg("eb0_t1"), py::arg("eb1_t1"), py::arg("min_distance") = 0.0,
        py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m.def(
        "point_triangle_ccd",
        [](const Eigen::Vector3d& p_t0, const Eigen::Vector3d& t0_t0,
           const Eigen::Vector3d& t1_t0, const Eigen::Vector3d& t2_t0,
           const Eigen::Vector3d& p_t1, const Eigen::Vector3d& t0_t1,
           const Eigen::Vector3d& t1_t1, const Eigen::Vector3d& t2_t1,
           const double min_distance = 0.0, const double tmax = 1.0,
           const double tolerance = DEFAULT_CCD_TOLERANCE,
           const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_triangle_ccd(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi,
                min_distance, tmax, tolerance, max_iterations,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Computes the time of impact between a point and a triangle in 3D using continuous collision detection.

        Parameters:
            p_t0: The initial position of the point.
            t0_t0: The initial position of the first vertex of the triangle.
            t1_t0: The initial position of the second vertex of the triangle.
            t2_t0: The initial position of the third vertex of the triangle.
            p_t1: The final position of the point.
            t0_t1: The final position of the first vertex of the triangle.
            t1_t1: The final position of the second vertex of the triangle.
            t2_t1: The final position of the third vertex of the triangle.
            min_distance: The minimum distance between the objects.
            tmax: The maximum time to check for collisions.
            tolerance: The error tolerance for the time of impact.
            max_iterations: The maximum number of iterations to perform.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            The time of impact between the point and the triangle.
        )ipc_Qu8mg5v7",
        py::arg("p_t0"), py::arg("t0_t0"), py::arg("t1_t0"), py::arg("t2_t0"),
        py::arg("p_t1"), py::arg("t0_t1"), py::arg("t1_t1"), py::arg("t2_t1"),
        py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
        py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);

    m.def(
        "ccd_strategy",
        [](const std::function<bool(long, double, bool, double&)>& ccd,
           const long max_iterations, const double min_distance,
           const double initial_distance, const double conservative_rescaling) {
            double toi;
            bool r = ccd_strategy(
                ccd, max_iterations, min_distance, initial_distance,
                conservative_rescaling, toi);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Perform the CCD strategy outlined by Li et al. [2020].

        Parameters:
            ccd: The continuous collision detection function.
            max_iterations: The maximum number of iterations to perform.
            min_distance: The minimum distance between the objects.
            initial_distance: The initial distance between the objects.
            conservative_rescaling: The conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            Output time of impact.
        )ipc_Qu8mg5v7",
        py::arg("ccd"), py::arg("max_iterations"), py::arg("min_distance"),
        py::arg("initial_distance"), py::arg("conservative_rescaling"));

    m.def(
        "check_initial_distance",
        [](const double initial_distance, const double min_distance) {
            double toi;
            bool r =
                check_initial_distance(initial_distance, min_distance, toi);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Helper function to check if the initial distance is less than the minimum distance.

        Parameters:
            initial_distance: The initial distance between the objects.
            min_distance: The minimum distance between the objects.

        Returns:
            Tuple of:
            True if the initial distance is less than the minimum distance, false otherwise.
            Set to 0 if the initial distance is less than the minimum distance.
        )ipc_Qu8mg5v7",
        py::arg("initial_distance"), py::arg("min_distance"));
}
