#include <common.hpp>

#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>

namespace py = pybind11;
using namespace ipc;

void define_sweep_and_tiniest_queue(py::module_& m)
{
    py::class_<CopyMeshBroadPhase, BroadPhase>(m, "CopyMeshBroadPhase");

    py::class_<SweepAndTiniestQueue, CopyMeshBroadPhase>(
        m, "SweepAndTiniestQueue");

#ifdef IPC_TOOLKIT_WITH_CUDA
    py::class_<SweepAndTiniestQueueGPU, CopyMeshBroadPhase>(
        m, "SweepAndTiniestQueueGPU");
#endif
}
