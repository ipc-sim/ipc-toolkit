#include <common.hpp>

#include <ipc/utils/cuda/is_available.hpp>

void define_cuda_is_available(py::module_& m) // m is the ipctk.cuda submodule
{
    m.def(
        "is_available", &ipc::cuda::is_available,
        R"ipc_Qu8mg5v7(
        Check whether the ipctk.cuda implementations can run here.

        A CUDA build still imports on machines without a GPU or an NVIDIA
        driver; its ipctk.cuda calls then raise on first use. Call this first
        to choose between the GPU and CPU implementations.

        Returns:
            True if this is a CUDA build and the CUDA runtime finds at least
            one usable device.
        )ipc_Qu8mg5v7");
}
