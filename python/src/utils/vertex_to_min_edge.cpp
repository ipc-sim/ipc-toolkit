#include <common.hpp>

#include <ipc/utils/vertex_to_min_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_to_min_edge(py::module_& m)
{
    m.def(
        "vertex_to_min_edge", &vertex_to_min_edge, py::arg("num_vertices"),
        py::arg("edges"));
}
