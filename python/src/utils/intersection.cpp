#include <common.hpp>

#include <ipc/utils/intersection.hpp>

namespace py = pybind11;
using namespace ipc;

void define_intersection(py::module_& m)
{
    m.def(
        "is_edge_intersecting_triangle", &is_edge_intersecting_triangle, "",
        py::arg("e0"), py::arg("e1"), py::arg("t0"), py::arg("t1"),
        py::arg("t2"));
}
