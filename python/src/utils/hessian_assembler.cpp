#include <common.hpp>

#include <ipc/utils/hessian_assembler.hpp>
#include <ipc/utils/meshfem_hessian_assembler.hpp>

using namespace ipc;

void define_hessian_assembler(py::module_& m)
{
    // NOTE: HessianAssembler is intentionally not extendable from Python. The
    // driver calls add_local_hessian() once per collision, so a Python-defined
    // assembler would acquire the GIL hundreds of thousands of times per
    // assembly. Only the built-in backends are exposed.
    py::class_<HessianAssembler>(
        m, "HessianAssembler", py::is_final(), R"ipc_Qu8mg5v7(
        Abstract sink for assembling local (per-collision) Hessians into a global matrix.

        Cannot be constructed or subclassed from Python; use TripletHessianAssembler or MeshFEMHessianAssembler.
        )ipc_Qu8mg5v7");

    py::class_<TripletHessianAssembler, HessianAssembler>(
        m, "TripletHessianAssembler", py::is_final(), R"ipc_Qu8mg5v7(
        Assembles through thread-local triplet caches and Eigen's setFromTriplets.
        )ipc_Qu8mg5v7")
        .def(py::init())
        .def(
            "get_matrix", &TripletHessianAssembler::get_matrix,
            R"ipc_Qu8mg5v7(
            Merge the thread-local caches and build the global matrix.

            Call once, after assembly; the internal caches are consumed.

            Returns:
                The assembled Hessian.
            )ipc_Qu8mg5v7");

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE
    py::class_<MeshFEMHessianAssembler, HessianAssembler>(
        m, "MeshFEMHessianAssembler", py::is_final(), R"ipc_Qu8mg5v7(
        Assembles into MeshFEMSparse's block-CSC data structures.

        Reuse one instance across assemblies (e.g., one per Newton solve) to also reuse the sparsity pattern, which is where most of the speedup over the triplet assembler comes from.
        )ipc_Qu8mg5v7")
        .def(py::init())
        .def(
            "get_matrix", &MeshFEMHessianAssembler::get_matrix,
            R"ipc_Qu8mg5v7(
            Convert the assembled matrix to a full symmetric sparse matrix.

            Returns:
                The assembled Hessian.
            )ipc_Qu8mg5v7")
        .def_property(
            "stale_block_tolerance",
            &MeshFEMHessianAssembler::stale_block_tolerance,
            &MeshFEMHessianAssembler::set_stale_block_tolerance,
            R"ipc_Qu8mg5v7(
            Number of vanished blocks tolerated before the pattern is rebuilt.

            A stencil that introduces a new vertex pair always forces a rebuild. When blocks merely disappear (e.g., a contact separates), the pattern is reused as long as at most this many blocks vanished. Stale blocks assemble to explicit zeros, so values stay correct. Defaults to 0.
            )ipc_Qu8mg5v7")
        .def_property(
            "assume_unchanged_stencils",
            &MeshFEMHessianAssembler::assume_unchanged_stencils,
            &MeshFEMHessianAssembler::set_assume_unchanged_stencils,
            R"ipc_Qu8mg5v7(
            Whether to skip change detection entirely.

            When enabled, the cached pattern is reused without comparing the stencils against it, as long as the stencil count is unchanged. Use this only when you know the collision set is identical to the previous assembly (e.g., reassembling with a different PSD projection or stiffness); on large scenes change detection costs about as much as rebuilding the pattern.

            Warning:
                If the stencils did change while the count stayed equal, assembly reads out of bounds.
            )ipc_Qu8mg5v7")
        .def_property_readonly(
            "reused_pattern", &MeshFEMHessianAssembler::reused_pattern,
            R"ipc_Qu8mg5v7(
            Whether the last assembly reused the cached sparsity pattern.
            )ipc_Qu8mg5v7");
#endif
}
