#include <common.hpp>

#include <ipc/ccd/nonlinear_ccd.hpp>

using namespace ipc;

class PyNonlinearTrajectory : public NonlinearTrajectory {
public:
    // clang-format off
    using NonlinearTrajectory::NonlinearTrajectory; // Inherit constructors
    VectorMax3d operator()(const double t) const override { PYBIND11_OVERRIDE_PURE_NAME(VectorMax3d, NonlinearTrajectory, "__call__", operator(), t); }
    double max_distance_from_linear(const double t0, const double t1) const override { PYBIND11_OVERRIDE_PURE(double, NonlinearTrajectory, max_distance_from_linear, t0, t1); }
    // clang-format on
};

#ifdef IPC_TOOLKIT_WITH_FILIB
class PyIntervalNonlinearTrajectory : public IntervalNonlinearTrajectory {
public:
    // clang-format off
    using IntervalNonlinearTrajectory::IntervalNonlinearTrajectory; // Inherit constructors
    VectorMax3d operator()(const double t) const override { PYBIND11_OVERRIDE_PURE_NAME(VectorMax3d, IntervalNonlinearTrajectory, "__call__", operator(), t); }
    VectorMax3I operator()(const filib::Interval& t) const override { PYBIND11_OVERRIDE_PURE_NAME(VectorMax3I, IntervalNonlinearTrajectory, "__call__", operator(), t); }
    double max_distance_from_linear(const double t0, const double t1) const override { PYBIND11_OVERRIDE(double, IntervalNonlinearTrajectory, max_distance_from_linear, t0, t1); }
    // clang-format on
};
#endif

void define_nonlinear_ccd(py::module_& m)
{
    py::class_<NonlinearTrajectory, PyNonlinearTrajectory>(
        m, "NonlinearTrajectory")
        .def(py::init<>())
        .def(
            "__call__", &NonlinearTrajectory::operator(),
            "Compute the point's position at time t", "t"_a)
        .def(
            "max_distance_from_linear",
            &NonlinearTrajectory::max_distance_from_linear,
            R"ipc_Qu8mg5v7(
            Compute the maximum distance from the nonlinear trajectory to a linearized trajectory

            Parameters:
                t0: Start time of the trajectory
                t1: End time of the trajectory
            )ipc_Qu8mg5v7",
            "t0"_a, "t1"_a);

#ifdef IPC_TOOLKIT_WITH_FILIB
    py::class_<
        IntervalNonlinearTrajectory, NonlinearTrajectory,
        PyIntervalNonlinearTrajectory>(m, "IntervalNonlinearTrajectory")
        .def(py::init<>())
        .def(
            "__call__",
            [](const IntervalNonlinearTrajectory& self, const double t) {
                return self(t);
            },
            "Compute the point's position at time t", "t"_a)
        .def(
            "__call__",
            py::overload_cast<const filib::Interval&>(
                &IntervalNonlinearTrajectory::operator(), py::const_),
            "Compute the point's position over a time interval t", "t"_a)
        .def(
            "max_distance_from_linear",
            &IntervalNonlinearTrajectory::max_distance_from_linear,
            R"ipc_Qu8mg5v7(
            Compute the maximum distance from the nonlinear trajectory to a linearized trajectory

            Note:
                This uses interval arithmetic to compute the maximum distance. If you know a tighter bound on the maximum distance, it is recommended to override this function.

            Parameters:
                t0: Start time of the trajectory
                t1: End time of the trajectory
            )ipc_Qu8mg5v7",
            "t0"_a, "t1"_a);
#endif

    py::class_<NonlinearCCD>(m, "NonlinearCCD")
        .def(
            py::init<const double, const long, const double>(),
            "tolerance"_a = NonlinearCCD::DEFAULT_TOLERANCE,
            "max_iterations"_a = NonlinearCCD::DEFAULT_MAX_ITERATIONS,
            "conservative_rescaling"_a =
                NonlinearCCD::DEFAULT_CONSERVATIVE_RESCALING)
        .def_readwrite(
            "tolerance", &NonlinearCCD::tolerance, "Solver tolerance.")
        .def_readwrite(
            "max_iterations", &NonlinearCCD::max_iterations,
            "Maximum number of iterations.")
        .def_readwrite(
            "conservative_rescaling", &NonlinearCCD::conservative_rescaling,
            "Conservative rescaling of the time of impact.")
        .def(
            "point_point_ccd",
            [](const NonlinearCCD& self, const NonlinearTrajectory& p0,
               const NonlinearTrajectory& p1, const double min_distance,
               const double tmax) {
                double toi;
                bool r = self.point_point_ccd(p0, p1, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform nonlinear CCD between two points moving along nonlinear trajectories.

            Parameters:
                p0: First point's trajectory
                p1: Second point's trajectory
                tmax: Maximum time to check for collision
                min_distance: Minimum separation distance between the two points

            Returns:
                Tuple of:
                True if the two points collide, false otherwise.
                Output time of impact
            )ipc_Qu8mg5v7",
            "p0"_a, "p1"_a, "min_distance"_a = 0, "tmax"_a = 1.0)
        .def(
            "point_edge_ccd",
            [](const NonlinearCCD& self, const NonlinearTrajectory& p,
               const NonlinearTrajectory& e0, const NonlinearTrajectory& e1,
               const double min_distance, const double tmax) {
                double toi;
                bool r =
                    self.point_edge_ccd(p, e0, e1, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform nonlinear CCD between a point and a linear edge moving along nonlinear trajectories.

            Parameters:
                p: Point's trajectory
                e0: Edge's first endpoint's trajectory
                e1: Edge's second endpoint's trajectory
                min_distance: Minimum separation distance between the point and the edge
                tmax: Maximum time to check for collision

            Returns:
                Tuple of:
                True if the point and edge collide, false otherwise.
                Output time of impact
            )ipc_Qu8mg5v7",
            "p"_a, "e0"_a, "e1"_a, "min_distance"_a = 0, "tmax"_a = 1.0)
        .def(
            "edge_edge_ccd",
            [](const NonlinearCCD& self, const NonlinearTrajectory& ea0,
               const NonlinearTrajectory& ea1, const NonlinearTrajectory& eb0,
               const NonlinearTrajectory& eb1, const double min_distance,
               const double tmax) {
                double toi;
                bool r = self.edge_edge_ccd(
                    ea0, ea1, eb0, eb1, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform nonlinear CCD between two linear edges moving along nonlinear trajectories.

            Parameters:
                ea0: First edge's first endpoint's trajectory
                ea1: First edge's second endpoint's trajectory
                eb0: Second edge's first endpoint's trajectory
                eb1: Second edge's second endpoint's trajectory
                min_distance: Minimum separation distance between the two edges
                tmax: Maximum time to check for collision

            Returns:
                Tuple of:
                True if the two edges collide, false otherwise.
                Output time of impact
            )ipc_Qu8mg5v7",
            "ea0"_a, "ea1"_a, "eb0"_a, "eb1"_a, "tmax"_a = 1.0,
            "min_distance"_a = 0)
        .def(
            "point_triangle_ccd",
            [](const NonlinearCCD& self, const NonlinearTrajectory& p,
               const NonlinearTrajectory& t0, const NonlinearTrajectory& t1,
               const NonlinearTrajectory& t2, const double min_distance,
               const double tmax) {
                double toi;
                bool r = self.point_triangle_ccd(
                    p, t0, t1, t2, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform nonlinear CCD between a point and a linear triangle moving along nonlinear trajectories.

            Parameters:
                p: Point's trajectory
                t0: Triangle's first vertex's trajectory
                t1: Triangle's second vertex's trajectory
                t2: Triangle's third vertex's trajectory
                min_distance: Minimum separation distance between the two edges
                tmax: Maximum time to check for collision

            Returns:
                Tuple of:
                True if the point and triangle collide, false otherwise.
                Output time of impact
            )ipc_Qu8mg5v7",
            "p"_a, "t0"_a, "t1"_a, "t2"_a, "min_distance"_a = 0, "tmax"_a = 1.0)
        .def_static(
            "conservative_piecewise_linear_ccd",
            [](const std::function<double(const double)>& distance,
               const std::function<double(const double, const double)>&
                   max_distance_from_linear,
               const std::function<bool(
                   const double, const double, const double, const bool,
                   double&)>& linear_ccd,
               const double min_distance, const double tmax,
               const double conservative_rescaling) {
                double toi;
                bool r = NonlinearCCD::conservative_piecewise_linear_ccd(
                    distance, max_distance_from_linear, linear_ccd, toi,
                    min_distance, tmax, conservative_rescaling);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform conservative piecewise linear CCD of a nonlinear trajectories.

            Parameters:
                distance: Return the distance for a given time in [0, 1].
                max_distance_from_linear: Return the maximum distance from the linearized trajectory for a given time interval.
                linear_ccd: Perform linear CCD on a given time interval.
                tmax: Maximum time to check for collision.
                min_distance: Minimum separation distance between the objects.
                conservative_rescaling: Conservative rescaling of the time of impact.

            Returns:
                Tuple of:

                Output time of impact.
            )ipc_Qu8mg5v7",
            "distance"_a, "max_distance_from_linear"_a, "linear_ccd"_a,
            "min_distance"_a = 0, "tmax"_a = 1.0,
            "conservative_rescaling"_a =
                TightInclusionCCD::DEFAULT_CONSERVATIVE_RESCALING);
}
