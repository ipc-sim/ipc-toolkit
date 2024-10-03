#include <common.hpp>

#include <ipc/broad_phase/sweep_and_prune.hpp>

namespace py = pybind11;
using namespace ipc;

void define_sweep_and_prune(py::module_& m)
{
    py::class_<SweepAndPrune, BroadPhase>(m, "SweepAndPrune").def(py::init());
}
