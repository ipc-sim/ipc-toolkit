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
               const double min_distance, const double tmax,
               const NarrowPhaseCCD& narrow_phase_ccd) {
                double toi;
                bool r = self.ccd(
                    vertices_t0, vertices_t1, toi, min_distance, tmax,
                    narrow_phase_ccd);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform narrow-phase CCD on the candidate.

            Parameters:
                vertices_t0: Stencil vertices at the start of the time step.
                vertices_t1: Stencil vertices at the end of the time step.
                min_distance: Minimum separation distance between primitives.
                tmax: Maximum time (normalized) to look for collisions. Should be in [0, 1].
                narrow_phase_ccd: The narrow phase CCD algorithm to use.

            Returns:
                Tuple of:
                If the candidate had a collision over the time interval.
                Computed time of impact (normalized).
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
            py::arg("narrow_phase_ccd") = DEFAULT_NARROW_PHASE_CCD)
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
