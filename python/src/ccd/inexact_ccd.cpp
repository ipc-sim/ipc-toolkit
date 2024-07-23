#include <common.hpp>

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
#include <ipc/ccd/inexact_ccd.hpp>
#else
namespace ipc {
}
#endif

namespace py = pybind11;
using namespace ipc;

void define_inexact_ccd(py::module_& m)
{
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
    py::class_<InexactCCD, NarrowPhaseCCD>(m, "InexactCCD")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a new AdditiveCCD object.

            Parameters:
                conservative_rescaling: The conservative rescaling of the time of impact.
            )ipc_Qu8mg5v7",
            py::arg("conservative_rescaling") =
                InexactCCD::DEFAULT_CONSERVATIVE_RESCALING)
        .def(
            "point_point_ccd_3D",
            [](const InexactCCD& self, const Eigen::Vector3d& p0_t0,
               const Eigen::Vector3d& p1_t0, const Eigen::Vector3d& p0_t1,
               const Eigen::Vector3d& p1_t1, const double min_distance = 0.0,
               const double tmax = 1.0) {
                double toi;
                bool r = self.point_point_ccd_3D(
                    p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Computes the time of impact between two points in 3D using continuous collision detection.

            Parameters:
                p0_t0: The initial position of the first point.
                p1_t0: The initial position of the second point.
                p0_t1: The final position of the first point.
                p1_t1: The final position of the second point.
                min_distance: The minimum distance between the objects.
                tmax: The maximum time to check for collisions.

            Returns:
                Tuple of:
                True if a collision was detected, false otherwise.
                The time of impact between the two points.
            )ipc_Qu8mg5v7",
            py::arg("p0_t0"), py::arg("p1_t0"), py::arg("p0_t1"),
            py::arg("p1_t1"), py::arg("min_distance") = 0.0,
            py::arg("tmax") = 1.0)
        .def(
            "point_edge_ccd_3D",
            [](const InexactCCD& self, const Eigen::Vector3d& p_t0,
               const Eigen::Vector3d& e0_t0, const Eigen::Vector3d& e1_t0,
               const Eigen::Vector3d& p_t1, const Eigen::Vector3d& e0_t1,
               const Eigen::Vector3d& e1_t1, const double min_distance = 0.0,
               const double tmax = 1.0) {
                double toi;
                bool r = self.point_edge_ccd_3D(
                    p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi, min_distance,
                    tmax);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Computes the time of impact between a point and an edge in 3D using continuous collision detection.

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
            py::arg("p_t0"), py::arg("e0_t0"), py::arg("e1_t0"),
            py::arg("p_t1"), py::arg("e0_t1"), py::arg("e1_t1"),
            py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0)
        .def_static(
            "ccd_strategy",
            [](const std::function<bool(double, double&)>& ccd,
               const double min_distance, const double initial_distance,
               const double conservative_rescaling) {
                double toi;
                static bool r = InexactCCD::ccd_strategy(
                    ccd, min_distance, initial_distance, conservative_rescaling,
                    toi);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform the CCD strategy outlined by Li et al. [2020].

            Parameters:
                ccd: The continuous collision detection function.
                min_distance: The minimum distance between the objects.
                initial_distance: The initial distance between the objects.
                conservative_rescaling: The conservative rescaling of the time of impact.

            Returns:
                Tuple of:
                True if a collision was detected, false otherwise.
                Output time of impact.
            )ipc_Qu8mg5v7",
            py::arg("ccd"), py::arg("min_distance"),
            py::arg("initial_distance"), py::arg("conservative_rescaling"))
        .def_readonly_static(
            "DEFAULT_CONSERVATIVE_RESCALING",
            &InexactCCD::DEFAULT_CONSERVATIVE_RESCALING,
            "The default conservative rescaling value used to avoid taking steps exactly to impact.")
        .def_readonly_static(
            "SMALL_TOI", &InexactCCD::SMALL_TOI,
            "Tolerance for small time of impact which triggers rerunning CCD without a minimum separation.")
        .def_readwrite(
            "conservative_rescaling", &InexactCCD::conservative_rescaling,
            "Conservative rescaling of the time of impact.");
#endif
}
