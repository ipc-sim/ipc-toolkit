#include <common.hpp>

#include <ipc/broad_phase/hash_grid.hpp>

using namespace ipc;

void define_hash_grid(py::module_& m)
{
    py::class_<HashItem>(m, "HashItem")
        .def(
            py::init<int, int>(),
            "Construct a hash item as a (key, value) pair.", "key"_a, "id"_a)
        .def(
            "__lt__", &HashItem::operator<,
            "Compare HashItems by their keys for sorting.", "other"_a)
        .def_readwrite("key", &HashItem::key, "The key of the item.")
        .def_readwrite("id", &HashItem::id, "The value of the item.");

    py::class_<HashGrid, BroadPhase, std::shared_ptr<HashGrid>>(m, "HashGrid")
        .def(py::init())
        .def_property_readonly("cell_size", &HashGrid::cell_size)
        .def_property_readonly(
            "grid_size", &HashGrid::grid_size,
            py::return_value_policy::reference)
        .def_property_readonly(
            "domain_min", &HashGrid::domain_min,
            py::return_value_policy::reference)
        .def_property_readonly(
            "domain_max", &HashGrid::domain_max,
            py::return_value_policy::reference);
}
