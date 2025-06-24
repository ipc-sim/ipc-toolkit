#include <common.hpp>

#include <ipc/utils/vertex_to_min_edge.hpp>

using namespace ipc;

void define_vertex_to_min_edge(py::module_& m)
{
    m.def(
        "vertex_to_min_edge", &vertex_to_min_edge, "num_vertices"_a, "edges"_a);
}
