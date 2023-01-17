#include "../common.hpp"

#include <ipc/candidates/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_candidate(py::module_& m)
{
    py::class_<EdgeEdgeCandidate, ContinuousCollisionCandidate>(
        m, "EdgeEdgeCandidate")
        .def(
            py::init<long, long>(), "", py::arg("edge0_index"),
            py::arg("edge1_index"))
        .def(
            "__str__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge0_index, ev.edge1_index);
            })
        .def(
            "__repr__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format(
                    "EdgeEdgeCandidate({:d}, {:d})", ev.edge0_index,
                    ev.edge1_index);
            })
        .def("num_vertices", &EdgeEdgeCandidate::num_vertices, "")
        .def(
            "vertex_indices", &EdgeEdgeCandidate::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &EdgeEdgeCandidate::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        .def(
            "compute_distance_gradient",
            &EdgeEdgeCandidate::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"),
            py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        .def(
            "compute_distance_hessian",
            &EdgeEdgeCandidate::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"),
            py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        .def(
            "ccd",
            [](EdgeEdgeCandidate& self, const Eigen::MatrixXd& V0,
               const Eigen::MatrixXd& V1, const Eigen::MatrixXi& E,
               const Eigen::MatrixXi& F, const double min_distance = 0.0,
               const double tmax = 1.0,
               const double tolerance = DEFAULT_CCD_TOLERANCE,
               const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
               const double conservative_rescaling =
                   DEFAULT_CCD_CONSERVATIVE_RESCALING) {
                double toi;
                bool r = self.ccd(
                    V0, V1, E, F, toi, tmax, tolerance, max_iterations,
                    conservative_rescaling);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform narrow-phase CCD on the candidate.

            Parameters:
                V0: Mesh vertex positions at the start of the time step.
                V1: Mesh vertex positions at the end of the time step.
                E: Mesh edges as rows of indicies into V.
                F: Mesh triangular faces as rows of indicies into V.
                tmax: Maximum time (normalized) to look for collisions. Should be in [0, 1].
                tolerance: CCD tolerance used by Tight-Inclusion CCD.
                max_iterations: Maximum iterations used by Tight-Inclusion CCD.
                conservative_rescaling: Conservative rescaling value used to avoid taking steps exactly to impact.

            Returns:
                Tuple of:
                If the candidate had a collision over the time interval.
                Computed time of impact (normalized).
            )ipc_Qu8mg5v7",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"),
            py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &EdgeEdgeCandidate::print_ccd_query, "",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"))
        .def("__eq__", &EdgeEdgeCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &EdgeEdgeCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &EdgeEdgeCandidate::operator<,
            "Compare EdgeEdgeCandidates for sorting.", py::arg("other"))
        .def_readwrite("edge0_index", &EdgeEdgeCandidate::edge0_index, "")
        .def_readwrite("edge1_index", &EdgeEdgeCandidate::edge1_index, "");
}
