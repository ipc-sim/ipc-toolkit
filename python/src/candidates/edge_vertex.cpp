#include <common.hpp>

#include <ipc/candidates/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_candidate(py::module_& m)
{
    py::class_<
        EdgeVertexCandidate, CollisionStencil, ContinuousCollisionCandidate>(
        m, "EdgeVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("edge_id"),
            py::arg("vertex_id"))
        .def(
            "__str__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format("[{:d}, {:d}]", ev.edge_id, ev.vertex_id);
            })
        .def(
            "__repr__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "EdgeVertexCandidate({:d}, {:d})", ev.edge_id,
                    ev.vertex_id);
            })
        .def("num_vertices", &EdgeVertexCandidate::num_vertices, "")
        .def(
            "vertex_ids", &EdgeVertexCandidate::vertex_ids, "",
            py::arg("edges"), py::arg("faces"))
        .def(
            "ccd",
            [](EdgeVertexCandidate& self, const Eigen::MatrixXd& vertices_t0,
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
            "print_ccd_query", &EdgeVertexCandidate::print_ccd_query, "",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"))
        .def("__eq__", &EdgeVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &EdgeVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &EdgeVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "edge_id", &EdgeVertexCandidate::edge_id, "ID of the edge")
        .def_readwrite(
            "vertex_id", &EdgeVertexCandidate::vertex_id, "ID of the vertex");
}
