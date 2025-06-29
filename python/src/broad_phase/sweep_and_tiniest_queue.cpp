#include <common.hpp>

#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
using namespace ipc; // not defined if IPC_TOOLKIT_WITH_CUDA is not defined
#endif

void define_sweep_and_tiniest_queue(py::module_& m)
{
#ifdef IPC_TOOLKIT_WITH_CUDA
    py::class_<
        SweepAndTiniestQueue, BroadPhase,
        std::shared_ptr<SweepAndTiniestQueue>>(m, "SweepAndTiniestQueue")
        .def(py::init());
#endif
}
