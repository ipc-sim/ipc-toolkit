#include <pybind11/pybind11.h>

#include <ipc/broad_phase/broad_phase.hpp>

namespace py = pybind11;
using namespace ipc;

void define_broad_phase_functions(py::module_& m)
{
    py::enum_<BroadPhaseMethod>(m, "BroadPhaseMethod")
        .value("BRUTE_FORCE", BroadPhaseMethod::BRUTE_FORCE)
        .value("HASH_GRID", BroadPhaseMethod::HASH_GRID)
        .value("SPATIAL_HASH", BroadPhaseMethod::SPATIAL_HASH)
        .export_values();
}
