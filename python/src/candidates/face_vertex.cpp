#include <common.hpp>

#include <ipc/candidates/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_candidate(py::module_& m)
{
    py::class_<
        FaceVertexCandidate, CollisionStencil, ContinuousCollisionCandidate>(
        m, "FaceVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("face_id"),
            py::arg("vertex_id"))
        .def(
            "__str__",
            [](const FaceVertexCandidate& ev) {
                return fmt::format("[{:d}, {:d}]", ev.face_id, ev.vertex_id);
            })
        .def(
            "__repr__",
            [](const FaceVertexCandidate& ev) {
                return fmt::format(
                    "FaceVertexCandidate({:d}, {:d})", ev.face_id,
                    ev.vertex_id);
            })
        .def("num_vertices", &FaceVertexCandidate::num_vertices, "")
        .def(
            "vertex_ids", &FaceVertexCandidate::vertex_ids, "",
            py::arg("edges"), py::arg("faces"))
        .def(
            "ccd",
            [](FaceVertexCandidate& self, const Eigen::MatrixXd& vertices_t0,
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
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("faces"), py::arg("min_distance") = 0.0,
            py::arg("tmax") = 1.0, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &FaceVertexCandidate::print_ccd_query, "",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"))
        .def("__eq__", &FaceVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &FaceVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &FaceVertexCandidate::operator<,
            "Compare FaceVertexCandidate for sorting.", py::arg("other"))
        .def_readwrite(
            "face_id", &FaceVertexCandidate::face_id, "ID of the face")
        .def_readwrite(
            "vertex_id", &FaceVertexCandidate::vertex_id, "ID of the vertex");
}
