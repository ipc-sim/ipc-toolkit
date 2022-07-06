#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/utils/world_bbox_diagonal_length.hpp>

namespace py = pybind11;
using namespace ipc;

void define_world_bbox_diagonal_length_members(py::module_& m)
{
    m.def(
        "world_bbox_diagonal_length", &world_bbox_diagonal_length, "",
        py::arg("V"));
}
