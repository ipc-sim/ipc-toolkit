#include "../common.hpp"

#include <ipc/candidates/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_candidate(py::module_& m)
{
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
        .def("num_vertices", &FaceVertexCandidate::num_vertices, "")
        .def(
            "vertex_indices", &FaceVertexCandidate::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &FaceVertexCandidate::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("dtype") = PointTriangleDistanceType::AUTO)
        .def(
            "compute_distance_gradient",
            &FaceVertexCandidate::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"),
            py::arg("dtype") = PointTriangleDistanceType::AUTO)
        .def(
            "compute_distance_hessian",
            &FaceVertexCandidate::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"),
            py::arg("dtype") = PointTriangleDistanceType::AUTO)
        .def(
            "ccd",
            [](FaceVertexCandidate& self, const Eigen::MatrixXd& V0,
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
            "", py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"),
            py::arg("min_distance") = 0.0, py::arg("tmax") = 1.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
            py::arg("conservative_rescaling") =
                DEFAULT_CCD_CONSERVATIVE_RESCALING)
        .def(
            "print_ccd_query", &FaceVertexCandidate::print_ccd_query, "",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"))
        .def("__eq__", &FaceVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &FaceVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &FaceVertexCandidate::operator<,
            "Compare FaceVertexCandidate for sorting.", py::arg("other"))
        .def_readwrite("face_index", &FaceVertexCandidate::face_index, "")
        .def_readwrite("vertex_index", &FaceVertexCandidate::vertex_index, "");
}
