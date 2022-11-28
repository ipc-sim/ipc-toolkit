#include "../common.hpp"

#include <ipc/candidates/candidates.hpp>

namespace py = pybind11;
using namespace ipc;

void define_candidates(py::module_& m)
{
    py::class_<Candidates>(m, "Candidates")
        .def(py::init(), "")
        .def("__len__", &Candidates::size, "")
        .def("empty", &Candidates::empty, "")
        .def("clear", &Candidates::clear, "")
        .def(
            "__getitem__",
            [](Candidates& self, size_t idx) -> ContinuousCollisionCandidate& {
                return self[idx];
            }, py::return_value_policy::reference)
        .def(
            "save_obj", &Candidates::save_obj, "", py::arg("filename"),
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def_readwrite("ev_candidates", &Candidates::ev_candidates, "")
        .def_readwrite("ee_candidates", &Candidates::ee_candidates, "")
        .def_readwrite("fv_candidates", &Candidates::fv_candidates, "");
}
