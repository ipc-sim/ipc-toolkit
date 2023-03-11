#include "../common.hpp"

#include <ipc/candidates/edge_face.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_face_candidate(py::module_& m)
{
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
}
