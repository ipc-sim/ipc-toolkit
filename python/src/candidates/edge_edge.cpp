#include <common.hpp>

#include <ipc/candidates/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_candidate(py::module_& m)
{
    py::class_<
        EdgeEdgeCandidate, CollisionStencil, ContinuousCollisionCandidate>(
        m, "EdgeEdgeCandidate")
        .def(
            py::init<long, long>(), "", py::arg("edge0_id"),
            py::arg("edge1_id"))
        .def(
            "__str__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format("[{:d}, {:d}]", ev.edge0_id, ev.edge1_id);
            })
        .def(
            "__repr__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format(
                    "EdgeEdgeCandidate({:d}, {:d})", ev.edge0_id, ev.edge1_id);
            })
        .def("num_vertices", &EdgeEdgeCandidate::num_vertices, "")
        .def(
            "vertex_ids", &EdgeEdgeCandidate::vertex_ids, "", py::arg("edges"),
            py::arg("faces"))
        .def(
            "ccd",
            [](EdgeEdgeCandidate& self, const Eigen::MatrixXd& vertices_t0,
               const Eigen::MatrixXd& vertices_t1, const Eigen::MatrixXi& edges,
               const Eigen::MatrixXi& faces, const double min_distance = 0.0,
               const double tmax = 1.0,
               const double tolerance = DEFAULT_CCD_TOLERANCE,
               const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
               const double conservative_rescaling =
                   DEFAULT_CCD_CONSERVATIVE_RESCALING) {
                double toi;
                bool r = self.ccd(
                    vertices_t0, vertices_t1, edges, faces, toi, min_distance,
                    tmax, tolerance, max_iterations, conservative_rescaling);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Perform narrow-phase CCD on the candidate.

            Parameters:
                vertices_t0: Mesh vertices at the start of the time step.
                vertices_t1: Mesh vertices at the end of the time step.
                edges: Collision mesh edges as rows of indicies into vertices.
                faces: Collision mesh triangular faces as rows of indicies into vertices.
                tmax: Maximum time (normalized) to look for collisions. Should be in [0, 1].
                tolerance: CCD tolerance used by Tight-Inclusion CCD.
                max_iterations: Maximum iterations used by Tight-Inclusion CCD.
                conservative_rescaling: Conservative rescaling value used to avoid taking steps exactly to impact.

            Returns:
                Tuple of:
                If the candidate had a collision over the time interval.
                Computed time of impact (normalized).
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"), py::arg("min_distance") = 0.0,
            py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &EdgeEdgeCandidate::print_ccd_query, "",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"))
        .def("__eq__", &EdgeEdgeCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &EdgeEdgeCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &EdgeEdgeCandidate::operator<,
            "Compare EdgeEdgeCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "edge0_id", &EdgeEdgeCandidate::edge0_id, "ID of the first edge")
        .def_readwrite(
            "edge1_id", &EdgeEdgeCandidate::edge1_id, "ID of the second edge");
}
