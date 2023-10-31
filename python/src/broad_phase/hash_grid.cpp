#include <common.hpp>

#include <ipc/broad_phase/hash_grid.hpp>

namespace py = pybind11;
using namespace ipc;

void define_hash_grid(py::module_& m)
{
    py::class_<HashItem>(m, "HashItem")
        .def(
            py::init<int, int>(),
            "Construct a hash item as a (key, value) pair.", py::arg("key"),
            py::arg("id"))
        .def(
            "__lt__", &HashItem::operator<,
            "Compare HashItems by their keys for sorting.", py::arg("other"))
        .def_readwrite("key", &HashItem::key, "The key of the item.")
        .def_readwrite("id", &HashItem::id, "The value of the item.");

    py::class_<HashGrid, BroadPhase>(m, "HashGrid")
        .def("cellSize", &HashGrid::cellSize)
        .def(
            "gridSize", &HashGrid::gridSize, py::return_value_policy::reference)
        .def(
            "domainMin", &HashGrid::domainMin,
            py::return_value_policy::reference)
        .def(
            "domainMax", &HashGrid::domainMax,
            py::return_value_policy::reference);
}
