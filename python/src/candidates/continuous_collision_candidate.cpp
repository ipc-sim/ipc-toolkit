#include <common.hpp>

#include <ipc/candidates/continuous_collision_candidate.hpp>

namespace py = pybind11;
using namespace ipc;

void define_continuous_collision_candidate(py::module_& m)
{
    py::class_<ContinuousCollisionCandidate>(m, "ContinuousCollisionCandidate")
        .def(
            "ccd",
            [](const ContinuousCollisionCandidate& self,
               const VectorMax12d& vertices_t0, const VectorMax12d& vertices_t1,
               const double min_distance = 0.0, const double tmax = 1.0,
               const double tolerance = DEFAULT_CCD_TOLERANCE,
               const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
               const double conservative_rescaling =
                   DEFAULT_CCD_CONSERVATIVE_RESCALING) {
                double toi;
                bool r = self.ccd(
                    vertices_t0, vertices_t1, toi, min_distance, tmax,
                    tolerance, max_iterations, conservative_rescaling);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform narrow-phase CCD on the candidate.

            Parameters:
                vertices_t0: Stencil vertices at the start of the time step.
                vertices_t1: Stencil vertices at the end of the time step.
                toi: Computed time of impact (normalized).
                min_distance: Minimum separation distance between primitives.
                tmax: Maximum time (normalized) to look for collisions. Should be in [0, 1].
                tolerance: CCD tolerance used by Tight-Inclusion CCD.
                max_iterations: Maximum iterations used by Tight-Inclusion CCD.
                conservative_rescaling: Conservative rescaling value used to avoid taking steps exactly to impact.

            Returns:
                Tuple of:
                If the candidate had a collision over the time interval.
                Computed time of impact (normalized).
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query",
            [](const ContinuousCollisionCandidate& self,
               const VectorMax12d& vertices_t0,
               const VectorMax12d& vertices_t1) -> void {
                self.write_ccd_query(std::cout, vertices_t0, vertices_t1);
            },
            R"ipc_Qu8mg5v7(
            Print the CCD query to cout.

            Parameters:
                vertices_t0: Stencil vertices at the start of the time step.
                vertices_t1: Stencil vertices at the end of the time step.
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"));
}
