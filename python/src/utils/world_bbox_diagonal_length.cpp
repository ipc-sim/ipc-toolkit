#include "../common.hpp"

#include <ipc/utils/world_bbox_diagonal_length.hpp>

namespace py = pybind11;
using namespace ipc;

void define_world_bbox_diagonal_length(py::module_& m)
{
    m.def(
        "world_bbox_diagonal_length", &world_bbox_diagonal_length, "",
        py::arg("V"));
}
