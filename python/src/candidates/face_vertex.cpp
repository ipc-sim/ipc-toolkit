#include <common.hpp>

#include <ipc/candidates/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_candidate(py::module_& m)
{
    py::class_<
        FaceVertexCandidate, CollisionStencil, ContinuousCollisionCandidate>(
        m, "FaceVertexCandidate")
        .def(py::init<long, long>(), py::arg("face_id"), py::arg("vertex_id"))
        .def("known_dtype", &FaceVertexCandidate::known_dtype)
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
        .def("__eq__", &FaceVertexCandidate::operator==, py::arg("other"))
        .def("__ne__", &FaceVertexCandidate::operator!=, py::arg("other"))
        .def(
            "__lt__", &FaceVertexCandidate::operator<,
            "Compare FaceVertexCandidate for sorting.", py::arg("other"))
        .def_readwrite(
            "face_id", &FaceVertexCandidate::face_id, "ID of the face")
        .def_readwrite(
            "vertex_id", &FaceVertexCandidate::vertex_id, "ID of the vertex");
}
