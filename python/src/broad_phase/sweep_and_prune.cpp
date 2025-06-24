#include <common.hpp>

#include <ipc/broad_phase/sweep_and_prune.hpp>

using namespace ipc;

void define_sweep_and_prune(py::module_& m)
{
    py::class_<SweepAndPrune, BroadPhase, std::shared_ptr<SweepAndPrune>>(
        m, "SweepAndPrune")
        .def(py::init());
}
