#include "../common.hpp"

#include <ipc/broad_phase/collision_candidate.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_candidate(py::module_& m)
{
    py::class_<ContinuousCollisionCandidate>(m, "ContinuousCollisionCandidate")
        .def(
            "ccd",
            [](ContinuousCollisionCandidate& self, const Eigen::MatrixXd& V0,
               const Eigen::MatrixXd& V1, const Eigen::MatrixXi& E,
               const Eigen::MatrixXi& F, const double tmax = 1.0,
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
            py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &ContinuousCollisionCandidate::print_ccd_query,
            "", py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"));

    py::class_<VertexVertexCandidate>(m, "VertexVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_index"),
            py::arg("vertex1_index"))
        .def(
            "__str__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.vertex0_index, ev.vertex1_index);
            })
        .def(
            "__repr__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "VertexVertexCandidate({:d}, {:d})", ev.vertex0_index,
                    ev.vertex1_index);
            })
        .def("__eq__", &VertexVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &VertexVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &VertexVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "vertex0_index", &VertexVertexCandidate::vertex0_index, "")
        .def_readwrite(
            "vertex1_index", &VertexVertexCandidate::vertex1_index, "");

    py::class_<EdgeVertexCandidate, ContinuousCollisionCandidate>(
        m, "EdgeVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("edge_index"),
            py::arg("vertex_index"))
        .def(
            "__str__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge_index, ev.vertex_index);
            })
        .def(
            "__repr__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "EdgeVertexCandidate({:d}, {:d})", ev.edge_index,
                    ev.vertex_index);
            })
        .def("__eq__", &EdgeVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &EdgeVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &EdgeVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def(
            "ccd",
            [](EdgeVertexCandidate& self, const Eigen::MatrixXd& V0,
               const Eigen::MatrixXd& V1, const Eigen::MatrixXi& E,
               const Eigen::MatrixXi& F, const double tmax = 1.0,
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
            py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &EdgeVertexCandidate::print_ccd_query, "",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"))
        .def_readwrite("edge_index", &EdgeVertexCandidate::edge_index, "")
        .def_readwrite("vertex_index", &EdgeVertexCandidate::vertex_index, "");

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
        .def("__eq__", &EdgeEdgeCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &EdgeEdgeCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &EdgeEdgeCandidate::operator<,
            "Compare EdgeEdgeCandidates for sorting.", py::arg("other"))
        .def(
            "ccd",
            [](EdgeEdgeCandidate& self, const Eigen::MatrixXd& V0,
               const Eigen::MatrixXd& V1, const Eigen::MatrixXi& E,
               const Eigen::MatrixXi& F, const double tmax = 1.0,
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
            py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &EdgeEdgeCandidate::print_ccd_query, "",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"))
        .def_readwrite("edge0_index", &EdgeEdgeCandidate::edge0_index, "")
        .def_readwrite("edge1_index", &EdgeEdgeCandidate::edge1_index, "");

    py::class_<EdgeFaceCandidate>(m, "EdgeFaceCandidate")
        .def(
            py::init<long, long>(), "", py::arg("edge_index"),
            py::arg("face_index"))
        .def(
            "__str__",
            [](const EdgeFaceCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge_index, ev.face_index);
            })
        .def(
            "__repr__",
            [](const EdgeFaceCandidate& ev) {
                return fmt::format(
                    "EdgeFaceCandidate({:d}, {:d})", ev.edge_index,
                    ev.face_index);
            })
        .def("__eq__", &EdgeFaceCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &EdgeFaceCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &EdgeFaceCandidate::operator<,
            "Compare EdgeFaceCandidate for sorting.", py::arg("other"))
        .def_readwrite("edge_index", &EdgeFaceCandidate::edge_index, "")
        .def_readwrite("face_index", &EdgeFaceCandidate::face_index, "");

    py::class_<FaceVertexCandidate, ContinuousCollisionCandidate>(
        m, "FaceVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("face_index"),
            py::arg("vertex_index"))
        .def(
            "__str__",
            [](const FaceVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.face_index, ev.vertex_index);
            })
        .def(
            "__repr__",
            [](const FaceVertexCandidate& ev) {
                return fmt::format(
                    "FaceVertexCandidate({:d}, {:d})", ev.face_index,
                    ev.vertex_index);
            })
        .def("__eq__", &FaceVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &FaceVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &FaceVertexCandidate::operator<,
            "Compare FaceVertexCandidate for sorting.", py::arg("other"))
        .def(
            "ccd",
            [](FaceVertexCandidate& self, const Eigen::MatrixXd& V0,
               const Eigen::MatrixXd& V1, const Eigen::MatrixXi& E,
               const Eigen::MatrixXi& F, const double tmax = 1.0,
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
            "", py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"),
            py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &FaceVertexCandidate::print_ccd_query, "",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"))
        .def_readwrite("face_index", &FaceVertexCandidate::face_index, "")
        .def_readwrite("vertex_index", &FaceVertexCandidate::vertex_index, "");

    py::class_<Candidates>(m, "Candidates")
        .def(py::init(), "")
        .def("__len__", &Candidates::size, "")
        .def("empty", &Candidates::empty, "")
        .def("clear", &Candidates::clear, "")
        .def(
            "__getitem__",
            [](Candidates& self, size_t idx) -> ContinuousCollisionCandidate* {
                return &self[idx];
            })
        .def(
            "save_obj", &Candidates::save_obj, "", py::arg("filename"),
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def_readwrite("ev_candidates", &Candidates::ev_candidates, "")
        .def_readwrite("ee_candidates", &Candidates::ee_candidates, "")
        .def_readwrite("fv_candidates", &Candidates::fv_candidates, "");
}
