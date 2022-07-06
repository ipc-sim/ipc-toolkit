#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/utils/intersection.hpp>

namespace py = pybind11;
using namespace ipc;

void define_intersection_members(py::module_& m)
{
    m.def(
        "is_edge_intersecting_triangle", &is_edge_intersecting_triangle, "",
        py::arg("e0"), py::arg("e1"), py::arg("t0"), py::arg("t1"),
        py::arg("t2"));
}
