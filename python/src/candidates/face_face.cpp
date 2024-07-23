#include <common.hpp>

#include <ipc/candidates/face_face.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_face_candidate(py::module_& m)
{
    py::class_<FaceFaceCandidate>(m, "FaceFaceCandidate")
        .def(py::init<long, long>(), py::arg("face0_id"), py::arg("face1_id"))
        .def("__eq__", &FaceFaceCandidate::operator==, py::arg("other"))
        .def("__ne__", &FaceFaceCandidate::operator!=, py::arg("other"))
        .def(
            "__lt__", &FaceFaceCandidate::operator<,
            "Compare FaceFaceCandidate for sorting.", py::arg("other"))
        .def_readwrite(
            "face0_id", &FaceFaceCandidate::face0_id, "ID of the first face.")
        .def_readwrite(
            "face1_id", &FaceFaceCandidate::face1_id, "ID of the second face.");
}
