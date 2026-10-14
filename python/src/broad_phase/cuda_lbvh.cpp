#include <common.hpp>

#include <ipc/broad_phase/cuda/lbvh.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
using namespace ipc; // not defined if IPC_TOOLKIT_WITH_CUDA is not defined
#endif

void define_cuda_lbvh(py::module_& m) // m is the ipctk.cuda submodule
{
#ifdef IPC_TOOLKIT_WITH_CUDA
    py::class_<cuda::LBVH, BroadPhase, std::shared_ptr<cuda::LBVH>>(
        m, "LBVH",
        R"ipc_Qu8mg5v7(
        GPU Linear Bounding Volume Hierarchy (LBVH) broad phase (ipc::cuda::LBVH).

        Builds the vertex/edge/face AABBs and BVHs, and runs the traversal and
        mesh-connectivity filtering, on the device. Produces the same candidates
        as ipctk.LBVH for any vertex filter. Available only in CUDA builds.
        )ipc_Qu8mg5v7")
        .def(py::init())
        .def_property_readonly(
            "num_vertex_nodes", &cuda::LBVH::num_vertex_nodes,
            "Number of nodes in the vertex BVH (2 * n_leaves - 1, or 0).")
        .def_property_readonly(
            "num_edge_nodes", &cuda::LBVH::num_edge_nodes,
            "Number of nodes in the edge BVH (2 * n_leaves - 1, or 0).")
        .def_property_readonly(
            "num_face_nodes", &cuda::LBVH::num_face_nodes,
            "Number of nodes in the face BVH (2 * n_leaves - 1, or 0).");
#endif
}
